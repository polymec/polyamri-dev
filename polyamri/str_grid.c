// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/unordered_set.h"
#include "core/unordered_map.h"
#include "polyamri/str_grid.h"
#include "polyamri/str_grid_patch_filler.h"

#if POLYMEC_HAVE_OPENMP
#include <omp.h>
#endif

struct str_grid_t
{
  // Intrinsic metadata.
  int nx, ny, nz, px, py, pz;
  bool x_periodic, y_periodic, z_periodic;
  
  // Information about which patches are present.
  int_unordered_set_t* patches;
  int* patch_indices;

  // Patch fillers.
  int_ptr_unordered_map_t* patch_fillers;

  // This flag is set by str_grid_finalize() after a grid has been assembled.
  bool finalized;

  // A mapping from a grid's ghost fill tokens to lists of tokens for the 
  // underlying patch fill operations.
  int_ptr_unordered_map_t* tokens;
};

str_grid_t* str_grid_new(int nx, int ny, int nz, 
                         int px, int py, int pz, 
                         bool periodic_in_x, 
                         bool periodic_in_y, 
                         bool periodic_in_z)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  ASSERT(px > 0);
  ASSERT(py > 0);
  ASSERT(pz > 0);

  str_grid_t* grid = polymec_malloc(sizeof(str_grid_t));
  grid->nx = nx;
  grid->ny = ny;
  grid->nz = nz;
  grid->px = px;
  grid->py = py;
  grid->pz = pz;
  grid->x_periodic = periodic_in_x;
  grid->y_periodic = periodic_in_y;
  grid->z_periodic = periodic_in_z;
  grid->patches = int_unordered_set_new();
  grid->patch_indices = NULL;
  grid->patch_fillers = int_ptr_unordered_map_new();
  grid->finalized = false;
  grid->tokens = int_ptr_unordered_map_new();
  return grid;
}

static inline int patch_index(str_grid_t* grid, int i, int j, int k)
{
  return grid->ny*grid->nz*i + grid->nz*j + k;
}

static inline void get_patch_indices(str_grid_t* grid, int index, int* i, int* j, int *k)
{
  *i = index / (grid->ny*grid->nz);
  *j = (index - grid->ny*grid->nz*(*i) ) / grid->nz;
  *k = index - grid->ny*grid->nz*(*i) - grid->nz*(*j);
}

void str_grid_free(str_grid_t* grid)
{
  int_ptr_unordered_map_free(grid->tokens);
  int_ptr_unordered_map_free(grid->patch_fillers);
  int_unordered_set_free(grid->patches);
  if (grid->patch_indices != NULL)
    polymec_free(grid->patch_indices);
  polymec_free(grid);
}

void str_grid_insert_patch(str_grid_t* grid, int i, int j, int k)
{
  ASSERT(!grid->finalized);
  int index = patch_index(grid, i, j, k);
  ASSERT(!int_unordered_set_contains(grid->patches, index));
  int_unordered_set_insert(grid->patches, index);
}

void str_grid_append_patch_filler(str_grid_t* grid, int i, int j, int k,
                                  str_grid_patch_filler_t* patch_filler)
{
  int index = patch_index(grid, i, j, k);
  ASSERT(int_unordered_set_contains(grid->patches, index));
  ptr_array_t** fillers_p = (ptr_array_t**)int_ptr_unordered_map_get(grid->patch_fillers, index);
  ptr_array_t* fillers;
  if (fillers_p != NULL)
    fillers = *fillers_p;
  else
  {
    fillers = ptr_array_new();
    int_ptr_unordered_map_insert_with_v_dtor(grid->patch_fillers, index, fillers, DTOR(ptr_array_free));
  }
  ptr_array_append(fillers, patch_filler);
}

void str_grid_finalize(str_grid_t* grid)
{
  ASSERT(!grid->finalized);

  // Make an array of indices for locally-present patches.
  grid->patch_indices = polymec_malloc(sizeof(int) * 3 * grid->patches->size);
  int l = 0;
  for (int i = 0; i < grid->nx; ++i)
  {
    for (int j = 0; j < grid->ny; ++j)
    {
      for (int k = 0; k < grid->nz; ++k, ++l)
      {
        grid->patch_indices[3*l]   = i;
        grid->patch_indices[3*l+1] = j;
        grid->patch_indices[3*l+2] = k;
      }
    }
  }

  grid->finalized = true;
}

bool str_grid_next_patch(str_grid_t* grid, int* pos, int* i, int* j, int* k)
{
  ASSERT(grid->finalized);
  ASSERT(*pos >= 0);
#if POLYMEC_HAVE_OPENMP
  int num_threads = omp_get_num_threads();
  int tid = omp_get_thread_num();
#else
  int num_threads = 1;
  int tid = 0;
#endif
  if (*pos == 0) 
    *pos = tid;
  bool result = (*pos < grid->patches->size);
  if (result)
  {
    int l = *pos;
    *i = grid->patch_indices[3*l];
    *j = grid->patch_indices[3*l+1];
    *k = grid->patch_indices[3*l+2];
    *pos += num_threads;
  }
  return result;
}

void str_grid_get_periodicity(str_grid_t* grid, bool* periodicity)
{
  ASSERT(periodicity != NULL);
  periodicity[0] = grid->x_periodic;
  periodicity[1] = grid->y_periodic;
  periodicity[2] = grid->z_periodic;
}

void str_grid_get_extents(str_grid_t* grid, int* nx, int* ny, int* nz)
{
  *nx = grid->nx;
  *ny = grid->ny;
  *nz = grid->nz;
}

void str_grid_get_patch_size(str_grid_t* grid, int* pnx, int* pny, int* pnz)
{
  *pnx = grid->px;
  *pny = grid->py;
  *pnz = grid->pz;
}

int str_grid_num_patches(str_grid_t* grid)
{
  return grid->patches->size;
}

bool str_grid_has_patch(str_grid_t* grid, int i, int j, int k)
{
  int index = patch_index(grid, i, j, k);
  return int_unordered_set_contains(grid->patches, index);
}

int str_grid_start_filling_ghost_cells(str_grid_t* grid, str_grid_cell_data_t* data)
{
  ASSERT(grid->finalized);
  ASSERT(str_grid_cell_data_grid(data) == grid);

  // Dream up a new token for ghost fills.
  int token = 0;
  while (int_ptr_unordered_map_contains(grid->tokens, token))
    ++token;

  // Map our token to a list of underlying operation tokens.`
  int_array_t* op_tokens = int_array_new();
  int_ptr_unordered_map_insert_with_v_dtor(grid->tokens, token, op_tokens, DTOR(int_array_free));

  // Begin the ghost fill operations and gather tokens.
  int pos = 0, patch_index;
  ptr_array_t* fillers;
  while (int_ptr_unordered_map_next(grid->patch_fillers, &pos, &patch_index, (void**)&fillers))
  {
    int i, j, k;
    get_patch_indices(grid, patch_index, &i, &j, &k);

    // Go through all the fillers for this patch.
    for (int l = 0; l < fillers->size; ++l)
    {
      str_grid_patch_filler_t* filler = fillers->data[l];
      int op_token = str_grid_patch_filler_start(filler, i, j, k, data);
      ASSERT((op_token >= 0) || (op_token == -1));

      if (op_token != -1) // -1 <-> no finish operation.
      {
        // We stash both the patch index and the op token.
        int_array_append(op_tokens, patch_index);
        int_array_append(op_tokens, op_token);
      }
    }
  }

  return token;
}

void str_grid_finish_filling_ghost_cells(str_grid_t* grid, int token)
{
  ASSERT(grid->finalized);

  // Retrieve the operation tokens corresponding to this one and 
  // finish out the fill operations.
  int_array_t** op_tokens_p = (int_array_t**)int_ptr_unordered_map_get(grid->tokens, token);
  ASSERT(op_tokens_p != NULL); // Token exists in our list?
  int_array_t* op_tokens = *op_tokens_p;
  ASSERT((op_tokens->size % 2) == 0); // Should be patch index/op token pairs
  int n = op_tokens->size/2;
  for (int l = 0; l < n; ++l)
  {
    int patch_index = op_tokens->data[2*n];
    int op_token = op_tokens->data[2*n+1];
    ASSERT(op_token >= 0);
    str_grid_patch_filler_t* filler = (str_grid_patch_filler_t*)(*int_ptr_unordered_map_get(grid->patch_fillers, patch_index));
    str_grid_patch_filler_finish(filler, op_token);
  }

  // Remove this entry from our token mapping.
  int_ptr_unordered_map_delete(grid->tokens, token);
}

void str_grid_fill_ghost_cells(str_grid_t* grid, str_grid_cell_data_t* data)
{
  int token = str_grid_start_filling_ghost_cells(grid, data);
  str_grid_finish_filling_ghost_cells(grid, token);
}

