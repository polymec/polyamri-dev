// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "polyamri/str_grid_edge_data.h"

struct str_grid_edge_data_t 
{
  str_grid_t* grid;
  int nx, ny, nz, nc;
  real_t patch_lx, patch_ly, patch_lz;
  int_ptr_unordered_map_t* x_patches;
  int_ptr_unordered_map_t* y_patches;
  int_ptr_unordered_map_t* z_patches;
  int* patch_offsets;
  int num_edges;

  void* buffer;
  bool owns_buffer;
};

static inline int patch_index(str_grid_edge_data_t* edge_data, int i, int j, int k)
{
  return edge_data->ny*edge_data->nz*i + edge_data->nz*j + k;
}

static void count_edges(str_grid_edge_data_t* edge_data)
{
  int pos = 0, ip, jp, kp, l = 1;
  str_grid_patch_t* patch;
  edge_data->num_edges = 0;
  edge_data->patch_offsets[0] = 0;
  while (str_grid_edge_data_next_x_patch(edge_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int nx = patch->i2 - patch->i1;
    int ny = patch->j2 - patch->j1 + 1;
    int nz = patch->k2 - patch->k1 + 1;
    edge_data->num_edges += nx * ny * nz;
    edge_data->patch_offsets[l] = edge_data->num_edges * patch->nc;
    ++l;
  }

  pos = 0;
  while (str_grid_edge_data_next_y_patch(edge_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int nx = patch->i2 - patch->i1 + 1;
    int ny = patch->j2 - patch->j1;
    int nz = patch->k2 - patch->k1 + 1;
    edge_data->num_edges += nx * ny * nz;
    edge_data->patch_offsets[l] = edge_data->num_edges * patch->nc;
    ++l;
  }

  pos = 0;
  while (str_grid_edge_data_next_z_patch(edge_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int nx = patch->i2 - patch->i1 + 1;
    int ny = patch->j2 - patch->j1 + 1;
    int nz = patch->k2 - patch->k1;
    edge_data->num_edges += nx * ny * nz;
    edge_data->patch_offsets[l] = edge_data->num_edges * patch->nc;
    ++l;
  }
}

str_grid_edge_data_t* str_grid_edge_data_new(str_grid_t* grid, 
                                             int num_components)
{
  str_grid_edge_data_t* edge_data = 
    str_grid_edge_data_with_buffer(grid, num_components, NULL);
  int data_size = edge_data->num_edges * edge_data->nc;
  void* buffer = polymec_malloc(sizeof(real_t) * data_size);
  str_grid_edge_data_set_buffer(edge_data, buffer, true);
  return edge_data;
}

str_grid_edge_data_t* str_grid_edge_data_with_buffer(str_grid_t* grid, 
                                                     int num_components,
                                                     void* buffer)
{
  ASSERT(num_components > 0);

  // Allocate storage.
  int num_patches = 3 * str_grid_num_patches(grid);
  size_t patches_size = sizeof(str_grid_patch_t*) * num_patches;
  size_t edge_data_size = sizeof(str_grid_edge_data_t) + patches_size;
  str_grid_edge_data_t* edge_data = polymec_malloc(edge_data_size);
  edge_data->grid = grid;
  edge_data->nc = num_components;

  str_grid_get_extents(grid, &edge_data->nx, &edge_data->ny, &edge_data->nz);
  edge_data->x_patches = int_ptr_unordered_map_new();
  edge_data->y_patches = int_ptr_unordered_map_new();
  edge_data->z_patches = int_ptr_unordered_map_new();
  edge_data->patch_offsets = polymec_malloc(sizeof(int) * (3*num_patches+1));
  edge_data->buffer = NULL;

  // Now populate the patches for x-, y-, and z-edges.
  int px, py, pz;
  str_grid_get_patch_size(grid, &px, &py, &pz);
  edge_data->patch_lx = 1.0 / px;
  edge_data->patch_ly = 1.0 / py;
  edge_data->patch_lz = 1.0 / pz;
  int pos = 0, i, j, k;
  while (str_grid_next_patch(grid, &pos, &i, &j, &k))
  {
    int index = patch_index(edge_data, i, j, k);

    int_ptr_unordered_map_insert_with_v_dtor(edge_data->x_patches, index, 
                                             str_grid_patch_with_buffer(px+1, py, pz, num_components, 0, NULL),
                                             DTOR(str_grid_patch_free));

    int_ptr_unordered_map_insert_with_v_dtor(edge_data->y_patches, index, 
                                             str_grid_patch_with_buffer(px, py+1, pz, num_components, 0, NULL),
                                             DTOR(str_grid_patch_free));

    int_ptr_unordered_map_insert_with_v_dtor(edge_data->z_patches, index, 
                                             str_grid_patch_with_buffer(px, py, pz+1, num_components, 0, NULL),
                                             DTOR(str_grid_patch_free));
  }

  // Set the buffer.
  str_grid_edge_data_set_buffer(edge_data, buffer, false);

  return edge_data;
}

void str_grid_edge_data_free(str_grid_edge_data_t* edge_data)
{
  int_ptr_unordered_map_free(edge_data->x_patches);
  int_ptr_unordered_map_free(edge_data->y_patches);
  int_ptr_unordered_map_free(edge_data->z_patches);
  if (edge_data->owns_buffer)
    polymec_free(edge_data->buffer);
  polymec_free(edge_data);
}

int str_grid_edge_data_num_components(str_grid_edge_data_t* edge_data)
{
  return edge_data->nc;
}

int str_grid_edge_data_num_edges(str_grid_edge_data_t* edge_data)
{
  return edge_data->num_edges;
}

str_grid_t* str_grid_edge_data_grid(str_grid_edge_data_t* edge_data)
{
  return edge_data->grid;
}

str_grid_patch_t* str_grid_edge_data_x_patch(str_grid_edge_data_t* edge_data, int i, int j, int k)
{
  int index = patch_index(edge_data, i, j, k);
  str_grid_patch_t** pp = (str_grid_patch_t**)int_ptr_unordered_map_get(edge_data->x_patches, index);
  if (pp != NULL)
    return *pp;
  else
    return NULL;
}

str_grid_patch_t* str_grid_edge_data_y_patch(str_grid_edge_data_t* edge_data, int i, int j, int k)
{
  int index = patch_index(edge_data, i, j, k);
  str_grid_patch_t** pp = (str_grid_patch_t**)int_ptr_unordered_map_get(edge_data->y_patches, index);
  if (pp != NULL)
    return *pp;
  else
    return NULL;
}

str_grid_patch_t* str_grid_edge_data_z_patch(str_grid_edge_data_t* edge_data, int i, int j, int k)
{
  int index = patch_index(edge_data, i, j, k);
  str_grid_patch_t** pp = (str_grid_patch_t**)int_ptr_unordered_map_get(edge_data->z_patches, index);
  if (pp != NULL)
    return *pp;
  else
    return NULL;
}

bool str_grid_edge_data_next_x_patch(str_grid_edge_data_t* edge_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** x_patch,
                                     bbox_t* bbox)
{
  bool result = str_grid_next_patch(edge_data->grid, pos, i, j, k);
  if (result)
  {
    *x_patch = str_grid_edge_data_x_patch(edge_data, *i, *j, *k);
    if (bbox != NULL)
    {
      bbox->x1 = (*i) * edge_data->patch_lx;
      bbox->x2 = bbox->x1 + edge_data->patch_lx;
      bbox->y1 = (*j) * edge_data->patch_ly;
      bbox->y2 = bbox->y1 + edge_data->patch_ly;
      bbox->z1 = (*k) * edge_data->patch_lz;
      bbox->z2 = bbox->z1 + edge_data->patch_lz;
    }
  }
  return result;
}

bool str_grid_edge_data_next_y_patch(str_grid_edge_data_t* edge_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** y_patch,
                                     bbox_t* bbox)
{
  bool result = str_grid_next_patch(edge_data->grid, pos, i, j, k);
  if (result)
  {
    *y_patch = str_grid_edge_data_y_patch(edge_data, *i, *j, *k);
    if (bbox != NULL)
    {
      bbox->x1 = (*i) * edge_data->patch_lx;
      bbox->x2 = bbox->x1 + edge_data->patch_lx;
      bbox->y1 = (*j) * edge_data->patch_ly;
      bbox->y2 = bbox->y1 + edge_data->patch_ly;
      bbox->z1 = (*k) * edge_data->patch_lz;
      bbox->z2 = bbox->z1 + edge_data->patch_lz;
    }
  }
  return result;
}

bool str_grid_edge_data_next_z_patch(str_grid_edge_data_t* edge_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** z_patch,
                                     bbox_t* bbox)
{
  bool result = str_grid_next_patch(edge_data->grid, pos, i, j, k);
  if (result)
  {
    *z_patch = str_grid_edge_data_z_patch(edge_data, *i, *j, *k);
    if (bbox != NULL)
    {
      bbox->x1 = (*i) * edge_data->patch_lx;
      bbox->x2 = bbox->x1 + edge_data->patch_lx;
      bbox->y1 = (*j) * edge_data->patch_ly;
      bbox->y2 = bbox->y1 + edge_data->patch_ly;
      bbox->z1 = (*k) * edge_data->patch_lz;
      bbox->z2 = bbox->z1 + edge_data->patch_lz;
    }
  }
  return result;
}

void* str_grid_edge_data_buffer(str_grid_edge_data_t* edge_data)
{
  return edge_data->buffer;
}

void str_grid_edge_data_set_buffer(str_grid_edge_data_t* edge_data, 
                                   void* buffer, 
                                   bool assume_control)
{
  if ((edge_data->buffer != NULL) && edge_data->owns_buffer)
    polymec_free(edge_data->buffer);
  edge_data->buffer = buffer;
  edge_data->owns_buffer = assume_control;

  // Point the patches at the buffer.
  int pos = 0, l = 0, ip, jp, kp;
  str_grid_patch_t* patch;
  while (str_grid_edge_data_next_x_patch(edge_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int patch_offset = edge_data->patch_offsets[l];
    patch->data = &(((real_t*)buffer)[patch_offset]);
    ++l;
  }

  pos = 0;
  while (str_grid_edge_data_next_y_patch(edge_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int patch_offset = edge_data->patch_offsets[l];
    patch->data = &(((real_t*)buffer)[patch_offset]);
    ++l;
  }

  pos = 0;
  while (str_grid_edge_data_next_z_patch(edge_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int patch_offset = edge_data->patch_offsets[l];
    patch->data = &(((real_t*)buffer)[patch_offset]);
    ++l;
  }
}

