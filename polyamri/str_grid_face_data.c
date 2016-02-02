// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "polyamri/str_grid_face_data.h"

struct str_grid_face_data_t 
{
  str_grid_t* grid;
  int nx, ny, nz, nc;
  int_ptr_unordered_map_t* x_patches;
  int_ptr_unordered_map_t* y_patches;
  int_ptr_unordered_map_t* z_patches;
};

static inline int patch_index(str_grid_face_data_t* face_data, int i, int j, int k)
{
  return face_data->ny*face_data->nz*i + face_data->nz*j + k;
}

str_grid_face_data_t* str_grid_face_data_new(str_grid_t* grid, 
                                             int num_components)
{
  ASSERT(num_components > 0);

  // Allocate storage.
  int num_patches = str_grid_num_patches(grid);
  size_t patches_size = sizeof(str_grid_patch_t*) * num_patches;
  size_t face_data_size = sizeof(str_grid_face_data_t) + patches_size;
  str_grid_face_data_t* face_data = polymec_malloc(face_data_size);
  face_data->grid = grid;
  face_data->nc = num_components;
  int nx, ny, nz;
  str_grid_get_extents(grid, &nx, &ny, &nz);
  face_data->x_patches = int_ptr_unordered_map_new();
  face_data->y_patches = int_ptr_unordered_map_new();
  face_data->z_patches = int_ptr_unordered_map_new();

  // Now populate the patches for x-, y-, and z-faces.
  int px, py, pz;
  str_grid_get_patch_size(grid, &px, &py, &pz);
  int pos = 0, i, j, k;
  while (str_grid_next_patch(grid, &pos, &i, &j, &k))
  {
    int index = patch_index(face_data, i, j, k);

    int_ptr_unordered_map_insert_with_v_dtor(face_data->x_patches, index, 
                                             str_grid_patch_new(px+1, py, pz, num_components, 0),
                                             DTOR(str_grid_patch_free));

    int_ptr_unordered_map_insert_with_v_dtor(face_data->y_patches, index, 
                                             str_grid_patch_new(px, py+1, pz, num_components, 0),
                                             DTOR(str_grid_patch_free));

    int_ptr_unordered_map_insert_with_v_dtor(face_data->z_patches, index, 
                                             str_grid_patch_new(px, py, pz+1, num_components, 0),
                                             DTOR(str_grid_patch_free));
  }

  return face_data;
}

void str_grid_face_data_free(str_grid_face_data_t* face_data)
{
  int_ptr_unordered_map_free(face_data->x_patches);
  int_ptr_unordered_map_free(face_data->y_patches);
  int_ptr_unordered_map_free(face_data->z_patches);
  polymec_free(face_data);
}

int str_grid_face_data_num_components(str_grid_face_data_t* face_data)
{
  return face_data->nc;
}

str_grid_t* str_grid_face_data_grid(str_grid_face_data_t* face_data)
{
  return face_data->grid;
}

str_grid_patch_t* str_grid_face_data_x_patch(str_grid_face_data_t* face_data, int i, int j, int k)
{
  int index = patch_index(face_data, i, j, k);
  str_grid_patch_t** pp = (str_grid_patch_t**)int_ptr_unordered_map_get(face_data->x_patches, index);
  if (pp != NULL)
    return *pp;
  else
    return NULL;
}

str_grid_patch_t* str_grid_face_data_y_patch(str_grid_face_data_t* face_data, int i, int j, int k)
{
  int index = patch_index(face_data, i, j, k);
  str_grid_patch_t** pp = (str_grid_patch_t**)int_ptr_unordered_map_get(face_data->y_patches, index);
  if (pp != NULL)
    return *pp;
  else
    return NULL;
}

str_grid_patch_t* str_grid_face_data_z_patch(str_grid_face_data_t* face_data, int i, int j, int k)
{
  int index = patch_index(face_data, i, j, k);
  str_grid_patch_t** pp = (str_grid_patch_t**)int_ptr_unordered_map_get(face_data->z_patches, index);
  if (pp != NULL)
    return *pp;
  else
    return NULL;
}

bool str_grid_face_data_next_x_patch(str_grid_face_data_t* face_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** x_patch)
{
  bool result = str_grid_next_patch(face_data->grid, pos, i, j, k);
  if (result)
    *x_patch = str_grid_face_data_x_patch(face_data, *i, *j, *k);
  return result;
}

bool str_grid_face_data_next_y_patch(str_grid_face_data_t* face_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** y_patch)
{
  bool result = str_grid_next_patch(face_data->grid, pos, i, j, k);
  if (result)
    *y_patch = str_grid_face_data_y_patch(face_data, *i, *j, *k);
  return result;
}

bool str_grid_face_data_next_z_patch(str_grid_face_data_t* face_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** z_patch)
{
  bool result = str_grid_next_patch(face_data->grid, pos, i, j, k);
  if (result)
    *z_patch = str_grid_face_data_z_patch(face_data, *i, *j, *k);
  return result;
}

