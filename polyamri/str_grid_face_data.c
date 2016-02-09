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
  real_t patch_lx, patch_ly, patch_lz;
  int_ptr_unordered_map_t* x_patches;
  int_ptr_unordered_map_t* y_patches;
  int_ptr_unordered_map_t* z_patches;
  int* patch_offsets;
  int num_faces;

  void* buffer;
  bool owns_buffer;
};

static inline int patch_index(str_grid_face_data_t* face_data, int i, int j, int k)
{
  return face_data->ny*face_data->nz*i + face_data->nz*j + k;
}

static void count_faces(str_grid_face_data_t* face_data)
{
  int pos = 0, ip, jp, kp, l = 1;
  str_grid_patch_t* patch;
  face_data->num_faces = 0;
  face_data->patch_offsets[0] = 0;
  while (str_grid_face_data_next_x_patch(face_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int nx = patch->i2 - patch->i1 + 1;
    int ny = patch->j2 - patch->j1;
    int nz = patch->k2 - patch->k1;
    face_data->num_faces += nx * ny * nz;
    face_data->patch_offsets[l] = face_data->num_faces * patch->nc;
    ++l;
  }

  pos = 0;
  while (str_grid_face_data_next_y_patch(face_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int nx = patch->i2 - patch->i1;
    int ny = patch->j2 - patch->j1 + 1;
    int nz = patch->k2 - patch->k1;
    face_data->num_faces += nx * ny * nz;
    face_data->patch_offsets[l] = face_data->num_faces * patch->nc;
    ++l;
  }

  pos = 0;
  while (str_grid_face_data_next_z_patch(face_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int nx = patch->i2 - patch->i1;
    int ny = patch->j2 - patch->j1;
    int nz = patch->k2 - patch->k1 + 1;
    face_data->num_faces += nx * ny * nz;
    face_data->patch_offsets[l] = face_data->num_faces * patch->nc;
    ++l;
  }
}

str_grid_face_data_t* str_grid_face_data_new(str_grid_t* grid, 
                                             int num_components)
{
  str_grid_face_data_t* face_data = 
    str_grid_face_data_with_buffer(grid, num_components, NULL);
  int data_size = face_data->num_faces * face_data->nc;
  void* buffer = polymec_malloc(sizeof(real_t) * data_size);
  str_grid_face_data_set_buffer(face_data, buffer, true);
  return face_data;
}

str_grid_face_data_t* str_grid_face_data_with_buffer(str_grid_t* grid, 
                                                     int num_components,
                                                     void* buffer)
{
  ASSERT(num_components > 0);

  // Allocate storage.
  int num_patches = 3 * str_grid_num_patches(grid);
  size_t patches_size = sizeof(str_grid_patch_t*) * num_patches;
  size_t face_data_size = sizeof(str_grid_face_data_t) + patches_size;
  str_grid_face_data_t* face_data = polymec_malloc(face_data_size);
  face_data->grid = grid;
  face_data->nc = num_components;

  str_grid_get_extents(grid, &face_data->nx, &face_data->ny, &face_data->nz);
  face_data->x_patches = int_ptr_unordered_map_new();
  face_data->y_patches = int_ptr_unordered_map_new();
  face_data->z_patches = int_ptr_unordered_map_new();
  face_data->patch_offsets = polymec_malloc(sizeof(int) * (3*num_patches+1));
  face_data->buffer = NULL;

  // Now populate the patches for x-, y-, and z-faces.
  int px, py, pz;
  str_grid_get_patch_size(grid, &px, &py, &pz);
  face_data->patch_lx = 1.0 / face_data->nx;
  face_data->patch_ly = 1.0 / face_data->ny;
  face_data->patch_lz = 1.0 / face_data->nz;
  int pos = 0, i, j, k;
  while (str_grid_next_patch(grid, &pos, &i, &j, &k))
  {
    int index = patch_index(face_data, i, j, k);

    int_ptr_unordered_map_insert_with_v_dtor(face_data->x_patches, index, 
                                             str_grid_patch_with_buffer(px+1, py, pz, num_components, 0, NULL),
                                             DTOR(str_grid_patch_free));

    int_ptr_unordered_map_insert_with_v_dtor(face_data->y_patches, index, 
                                             str_grid_patch_with_buffer(px, py+1, pz, num_components, 0, NULL),
                                             DTOR(str_grid_patch_free));

    int_ptr_unordered_map_insert_with_v_dtor(face_data->z_patches, index, 
                                             str_grid_patch_with_buffer(px, py, pz+1, num_components, 0, NULL),
                                             DTOR(str_grid_patch_free));
  }
  count_faces(face_data);

  // Set the buffer.
  str_grid_face_data_set_buffer(face_data, buffer, false);

  return face_data;
}

void str_grid_face_data_free(str_grid_face_data_t* face_data)
{
  int_ptr_unordered_map_free(face_data->x_patches);
  int_ptr_unordered_map_free(face_data->y_patches);
  int_ptr_unordered_map_free(face_data->z_patches);
  polymec_free(face_data->patch_offsets);
  if (face_data->owns_buffer)
    polymec_free(face_data->buffer);
  polymec_free(face_data);
}

int str_grid_face_data_num_components(str_grid_face_data_t* face_data)
{
  return face_data->nc;
}

int str_grid_face_data_num_faces(str_grid_face_data_t* face_data)
{
  return face_data->num_faces;
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
                                     str_grid_patch_t** x_patch,
                                     bbox_t* bbox)
{
  bool result = str_grid_next_patch(face_data->grid, pos, i, j, k);
  if (result)
  {
    *x_patch = str_grid_face_data_x_patch(face_data, *i, *j, *k);
    if (bbox != NULL)
    {
      bbox->x1 = (*i) * face_data->patch_lx;
      bbox->x2 = bbox->x1 + face_data->patch_lx;
      bbox->y1 = (*j) * face_data->patch_ly;
      bbox->y2 = bbox->y1 + face_data->patch_ly;
      bbox->z1 = (*k) * face_data->patch_lz;
      bbox->z2 = bbox->z1 + face_data->patch_lz;
    }
  }
  return result;
}

bool str_grid_face_data_next_y_patch(str_grid_face_data_t* face_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** y_patch,
                                     bbox_t* bbox)
{
  bool result = str_grid_next_patch(face_data->grid, pos, i, j, k);
  if (result)
  {
    *y_patch = str_grid_face_data_y_patch(face_data, *i, *j, *k);
    if (bbox != NULL)
    {
      bbox->x1 = (*i) * face_data->patch_lx;
      bbox->x2 = bbox->x1 + face_data->patch_lx;
      bbox->y1 = (*j) * face_data->patch_ly;
      bbox->y2 = bbox->y1 + face_data->patch_ly;
      bbox->z1 = (*k) * face_data->patch_lz;
      bbox->z2 = bbox->z1 + face_data->patch_lz;
    }
  }
  return result;
}

bool str_grid_face_data_next_z_patch(str_grid_face_data_t* face_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** z_patch,
                                     bbox_t* bbox)
{
  bool result = str_grid_next_patch(face_data->grid, pos, i, j, k);
  if (result)
  {
    *z_patch = str_grid_face_data_z_patch(face_data, *i, *j, *k);
    if (bbox != NULL)
    {
      bbox->x1 = (*i) * face_data->patch_lx;
      bbox->x2 = bbox->x1 + face_data->patch_lx;
      bbox->y1 = (*j) * face_data->patch_ly;
      bbox->y2 = bbox->y1 + face_data->patch_ly;
      bbox->z1 = (*k) * face_data->patch_lz;
      bbox->z2 = bbox->z1 + face_data->patch_lz;
    }
  }
  return result;
}

void* str_grid_face_data_buffer(str_grid_face_data_t* face_data)
{
  return face_data->buffer;
}

void str_grid_face_data_set_buffer(str_grid_face_data_t* face_data, 
                                   void* buffer, 
                                   bool assume_control)
{
  if ((face_data->buffer != NULL) && face_data->owns_buffer)
    polymec_free(face_data->buffer);
  face_data->buffer = buffer;
  face_data->owns_buffer = assume_control;

  // Point the patches at the buffer.
  int pos = 0, l = 0, ip, jp, kp;
  str_grid_patch_t* patch;
  while (str_grid_face_data_next_x_patch(face_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int patch_offset = face_data->patch_offsets[l];
    patch->data = &(((real_t*)buffer)[patch_offset]);
    ++l;
  }

  pos = 0;
  while (str_grid_face_data_next_y_patch(face_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int patch_offset = face_data->patch_offsets[l];
    patch->data = &(((real_t*)buffer)[patch_offset]);
    ++l;
  }

  pos = 0;
  while (str_grid_face_data_next_z_patch(face_data, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    int patch_offset = face_data->patch_offsets[l];
    patch->data = &(((real_t*)buffer)[patch_offset]);
    ++l;
  }
}

