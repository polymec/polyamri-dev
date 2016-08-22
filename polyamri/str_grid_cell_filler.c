// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "polyamri/str_grid_patch_filler.h"
#include "polyamri/str_grid_cell_filler.h"

struct str_grid_cell_filler_t 
{
  // Underlying grid and metadata.
  str_grid_t* grid;
  int nx, ny, nz;

  // Patch fillers.
  int_ptr_unordered_map_t* patch_fillers;

  // A mapping from a grid's ghost fill tokens to lists of tokens for the 
  // underlying patch fill operations.
  int_ptr_unordered_map_t* tokens;
};

str_grid_cell_filler_t* str_grid_cell_filler_new(str_grid_t* grid)
{
  str_grid_cell_filler_t* filler = polymec_malloc(sizeof(str_grid_cell_filler_t));
  filler->grid = grid;
  str_grid_get_extents(grid, &filler->nx, &filler->ny, &filler->nz);
  filler->patch_fillers = int_ptr_unordered_map_new();
  filler->tokens = int_ptr_unordered_map_new();
  return filler;
}

void str_grid_cell_filler_free(str_grid_cell_filler_t* cell_filler)
{
  int_ptr_unordered_map_free(cell_filler->tokens);
  int_ptr_unordered_map_free(cell_filler->patch_fillers);
  polymec_free(cell_filler);
}

static inline int patch_index(str_grid_cell_filler_t* cell_filler, int i, int j, int k)
{
  return cell_filler->ny*cell_filler->nz*i + cell_filler->nz*j + k;
}

static inline void get_patch_indices(str_grid_cell_filler_t* cell_filler, int index, int* i, int* j, int *k)
{
  *i = index / (cell_filler->ny*cell_filler->nz);
  *j = (index - cell_filler->ny*cell_filler->nz*(*i) ) / cell_filler->nz;
  *k = index - cell_filler->ny*cell_filler->nz*(*i) - cell_filler->nz*(*j);
}

void str_grid_cell_filler_insert(str_grid_cell_filler_t* cell_filler, 
                                 int i, int j, int k,
                                 str_grid_patch_filler_t* patch_filler)
{
  int index = patch_index(cell_filler, i, j, k);
  ASSERT(str_grid_has_patch(cell_filler->grid, i, j, k));
  ptr_array_t** fillers_p = (ptr_array_t**)int_ptr_unordered_map_get(cell_filler->patch_fillers, index);
  ptr_array_t* fillers;
  if (fillers_p != NULL)
    fillers = *fillers_p;
  else
  {
    fillers = ptr_array_new();
    int_ptr_unordered_map_insert_with_v_dtor(cell_filler->patch_fillers, index, fillers, DTOR(ptr_array_free));
  }
  ptr_array_append(fillers, patch_filler);
}

str_grid_t* str_grid_cell_filler_grid(str_grid_cell_filler_t* cell_filler)
{
  return cell_filler->grid;
}

void str_grid_cell_filler_fill(str_grid_cell_filler_t* cell_filler, str_grid_cell_data_t* data)
{
  int token = str_grid_cell_filler_start(cell_filler, data);
  str_grid_cell_filler_finish(cell_filler, token);
}

int str_grid_cell_filler_start(str_grid_cell_filler_t* cell_filler, str_grid_cell_data_t* data)
{
  ASSERT(str_grid_cell_data_grid(data) == cell_filler->grid);

  // Dream up a new token for ghost fills.
  int token = 0;
  while (int_ptr_unordered_map_contains(cell_filler->tokens, token))
    ++token;

  // Map our token to a list of underlying operation tokens.`
  int_array_t* op_tokens = int_array_new();
  int_ptr_unordered_map_insert_with_v_dtor(cell_filler->tokens, token, op_tokens, DTOR(int_array_free));

  // Begin the ghost fill operations and gather tokens.
  int pos = 0, patch_index;
  ptr_array_t* fillers;
  while (int_ptr_unordered_map_next(cell_filler->patch_fillers, &pos, &patch_index, (void**)&fillers))
  {
    int i, j, k;
    get_patch_indices(cell_filler, patch_index, &i, &j, &k);

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

void str_grid_cell_filler_finish(str_grid_cell_filler_t* cell_filler, int token)
{
  // Retrieve the operation tokens corresponding to this one and 
  // finish out the fill operations.
  int_array_t** op_tokens_p = (int_array_t**)int_ptr_unordered_map_get(cell_filler->tokens, token);
  ASSERT(op_tokens_p != NULL); // Token exists in our list?
  int_array_t* op_tokens = *op_tokens_p;
  ASSERT((op_tokens->size % 2) == 0); // Should be patch index/op token pairs
  int n = (int)op_tokens->size/2;
  for (int l = 0; l < n; ++l)
  {
    int patch_index = op_tokens->data[2*n];
    int op_token = op_tokens->data[2*n+1];
    ASSERT(op_token >= 0);
    str_grid_patch_filler_t* filler = (str_grid_patch_filler_t*)(*int_ptr_unordered_map_get(cell_filler->patch_fillers, patch_index));
    str_grid_patch_filler_finish(filler, op_token);
  }

  // Remove this entry from our token mapping.
  int_ptr_unordered_map_delete(cell_filler->tokens, token);
}

