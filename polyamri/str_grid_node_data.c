// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "polyamri/str_grid_node_data.h"

struct str_grid_node_data_t 
{
  str_grid_t* grid;
  int nx, ny, nz, nc;
  int_ptr_unordered_map_t* patches;
  int* patch_offsets;
  int num_nodes;

  void* buffer;
  bool owns_buffer;
};

static inline int patch_index(str_grid_node_data_t* node_data, int i, int j, int k)
{
  return node_data->ny*node_data->nz*i + node_data->nz*j + k;
}

static void count_nodes(str_grid_node_data_t* node_data)
{
  int pos = 0, ip, jp, kp, l = 1;
  str_grid_patch_t* patch;
  node_data->num_nodes = 0;
  node_data->patch_offsets[0] = 0;
  while (str_grid_node_data_next_patch(node_data, &pos, &ip, &jp, &kp, &patch))
  {
    int nx = patch->i2 - patch->i1 + 1;
    int ny = patch->j2 - patch->j1 + 1;
    int nz = patch->k2 - patch->k1 + 1;
    node_data->num_nodes += nx * ny * nz;
    node_data->patch_offsets[l] = node_data->num_nodes * patch->nc;
    ++l;
  }
}

str_grid_node_data_t* str_grid_node_data_new(str_grid_t* grid, 
                                             int num_components)
{
  str_grid_node_data_t* node_data = 
    str_grid_node_data_with_buffer(grid, num_components, NULL);
  int data_size = node_data->num_nodes * node_data->nc;
  void* buffer = polymec_malloc(sizeof(real_t) * data_size);
  str_grid_node_data_set_buffer(node_data, buffer, true);
  return node_data;
}

str_grid_node_data_t* str_grid_node_data_with_buffer(str_grid_t* grid, 
                                                     int num_components, 
                                                     void* buffer)
{
  ASSERT(num_components > 0);

  // Allocate storage.
  int num_patches = str_grid_num_patches(grid);
  size_t patches_size = sizeof(str_grid_patch_t*) * num_patches;
  size_t node_data_size = sizeof(str_grid_node_data_t) + patches_size;
  str_grid_node_data_t* node_data = polymec_malloc(node_data_size);
  node_data->grid = grid;
  node_data->nc = num_components;
  str_grid_get_extents(grid, &node_data->nx, &node_data->ny, &node_data->nz);
  node_data->patches = int_ptr_unordered_map_new();
  node_data->patch_offsets = polymec_malloc(sizeof(int) * (num_patches+1));

  // Now populate the patches (with NULL buffers).
  int px, py, pz;
  str_grid_get_patch_size(grid, &px, &py, &pz);
  int pos = 0, i, j, k, l = 0;
  while (str_grid_next_patch(grid, &pos, &i, &j, &k))
  {
    int index = patch_index(node_data, i, j, k);
    int_ptr_unordered_map_insert_with_v_dtor(node_data->patches, index, 
      str_grid_patch_with_buffer(px+1, py+1, pz+1, num_components, 0, NULL),
      DTOR(str_grid_patch_free));
    ++l;
  }

  count_nodes(node_data);

  // Set the buffer.
  str_grid_node_data_set_buffer(node_data, buffer, false);

  return node_data;
}

void str_grid_node_data_free(str_grid_node_data_t* node_data)
{
  int_ptr_unordered_map_free(node_data->patches);
  polymec_free(node_data->patch_offsets);
  if (node_data->owns_buffer)
    polymec_free(node_data->buffer);
  polymec_free(node_data);
}

int str_grid_node_data_num_components(str_grid_node_data_t* node_data)
{
  return node_data->nc;
}

int str_grid_node_data_num_patches(str_grid_node_data_t* node_data)
{
  return str_grid_num_patches(node_data->grid);
}

int str_grid_node_data_num_nodes(str_grid_node_data_t* node_data)
{
  return node_data->num_nodes;
}

str_grid_t* str_grid_node_data_grid(str_grid_node_data_t* node_data)
{
  return node_data->grid;
}

str_grid_patch_t* str_grid_node_data_patch(str_grid_node_data_t* node_data, int i, int j, int k)
{
  int index = patch_index(node_data, i, j, k);
  str_grid_patch_t** pp = (str_grid_patch_t**)int_ptr_unordered_map_get(node_data->patches, index);
  if (pp != NULL)
    return *pp;
  else
    return NULL;
}

bool str_grid_node_data_next_patch(str_grid_node_data_t* node_data, int* pos, 
                                   int* i, int* j, int* k, 
                                   str_grid_patch_t** patch)
{
  bool result = str_grid_next_patch(node_data->grid, pos, i, j, k);
  if (result)
    *patch = str_grid_node_data_patch(node_data, *i, *j, *k);
  return result;
}

void* str_grid_node_data_buffer(str_grid_node_data_t* node_data)
{
  return node_data->buffer;
}

void str_grid_node_data_set_buffer(str_grid_node_data_t* node_data, 
                                   void* buffer, 
                                   bool assume_control)
{
  if ((node_data->buffer != NULL) && node_data->owns_buffer)
    polymec_free(node_data->buffer);
  node_data->buffer = buffer;
  node_data->owns_buffer = assume_control;

  // Point the patches at the buffer.
  int pos = 0, l = 0, ip, jp, kp;
  str_grid_patch_t* patch;
  while (str_grid_node_data_next_patch(node_data, &pos, &ip, &jp, &kp, &patch))
  {
    int patch_offset = node_data->patch_offsets[l];
    patch->data = &(((real_t*)buffer)[patch_offset]);
    ++l;
  }
}

