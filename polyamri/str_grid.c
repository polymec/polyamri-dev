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

  // This flag is set by str_grid_finalize() after a grid has been assembled.
  bool finalized;

  // Properties stored in this grid.
  string_ptr_unordered_map_t* properties;
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
  grid->finalized = false;
  grid->properties = string_ptr_unordered_map_new();
  return grid;
}

static inline int patch_index(str_grid_t* grid, int i, int j, int k)
{
  return grid->ny*grid->nz*i + grid->nz*j + k;
}

void str_grid_free(str_grid_t* grid)
{
  string_ptr_unordered_map_free(grid->properties);
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

// Properties stashed on the grid.
typedef struct
{
  void* data;
  serializer_t* serializer;
} prop_t;

static prop_t* prop_new(void* data, serializer_t* ser)
{
  prop_t* prop = polymec_malloc(sizeof(prop_t));
  prop->data = data;
  prop->serializer = ser;
  return prop;
}

static void prop_dtor(void* context)
{
  prop_t* prop = context;
  if (prop->serializer != NULL)
    serializer_destroy_object(prop->serializer, prop->data);
  prop->data = NULL;
  prop->serializer = NULL;
  polymec_free(prop);
}

void str_grid_set_property(str_grid_t* grid, 
                           const char* property, 
                           void* data, 
                           serializer_t* serializer)
{
  prop_t* prop = prop_new(data, serializer);
  string_ptr_unordered_map_insert_with_kv_dtors(grid->properties, 
                                                string_dup(property),
                                                prop,
                                                string_free,
                                                prop_dtor);
}

void* str_grid_property(str_grid_t* grid, const char* property)
{
  void** prop_p = string_ptr_unordered_map_get(grid->properties, (char*)property);
  if (prop_p != NULL)
  {
    prop_t* prop = *((prop_t**)prop_p);
    return prop->data;
  }
  else
    return NULL;
}

void str_grid_delete_property(str_grid_t* grid, const char* property)
{
  string_ptr_unordered_map_delete(grid->properties, (char*)property);
}

bool str_grid_next_property(str_grid_t* grid, int* pos, 
                            char** prop_name, void** prop_data, 
                            serializer_t** prop_serializer)
{
  void* val;
  bool result = string_ptr_unordered_map_next(grid->properties, pos,
                                              prop_name, &val);
  if (result)
  {
    prop_t* prop = val;
    *prop_data = prop->data;
    *prop_serializer = prop->serializer;
  }
  return result;
}

