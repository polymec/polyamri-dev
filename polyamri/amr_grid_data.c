// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/amr_grid_data.h"

struct amr_grid_data_t
{
  amr_grid_t* grid;
  int nx, ny, nz;
  void* patches;
  void* bboxes;
};

// Multi-dimensional array accessors for patches, bboxes.
#define DECLARE_PATCH_ARRAY(array, grid_data) \
amr_patch_t* (*array)[grid_data->nx][grid_data->ny] = grid_data->patches

#define DECLARE_BBOX_ARRAY(array, grid_data) \
bbox_t* (*array)[grid_data->nx][grid_data->ny] = grid_data->bboxes

amr_grid_data_t* amr_grid_data_new(amr_grid_t* grid, int num_components)
{
  ASSERT(num_components > 0);

  // Allocate storage.
  int num_patches = amr_grid_num_patches(grid);
  size_t patches_size = sizeof(amr_patch_t*) * num_patches;
  size_t bboxes_size = sizeof(bbox_t*) * num_patches;
  size_t grid_data_size = sizeof(amr_grid_data_t) + patches_size + bboxes_size;
  amr_grid_data_t* grid_data = polymec_malloc(grid_data_size);
  grid_data->grid = grid;
  amr_grid_get_extents(grid, &grid_data->nx, &grid_data->ny, &grid_data->nz);
  grid_data->patches = (char*)grid_data + sizeof(amr_grid_data_t);
  memset(grid_data->patches, 0, patches_size);
  grid_data->bboxes = (char*)grid_data + sizeof(amr_grid_data_t) + patches_size;
  memset(grid_data->bboxes, 0, bboxes_size);

  // Now populate the patches/bboxes.
  bbox_t* domain = amr_grid_domain(grid);
  int px, py, pz, ng;
  amr_grid_get_patch_size(grid, &px, &py, &pz, &ng);
  DECLARE_PATCH_ARRAY(patches, grid_data);
  DECLARE_BBOX_ARRAY(bboxes, grid_data);
  real_t dx = (domain->x2 - domain->x1) / grid_data->nx;
  real_t dy = (domain->y2 - domain->y1) / grid_data->ny;
  real_t dz = (domain->z2 - domain->z1) / grid_data->nz;
  for (int i = 0; i < grid_data->nx; ++i)
  {
    for (int j = 0; j < grid_data->ny; ++j)
    {
      for (int k = 0; k < grid_data->nz; ++k)
      {
        if (amr_grid_has_patch(grid, i, j, k))
        {
          patches[i][j][k] = amr_patch_new(px, py, pz, num_components, ng);
          bboxes[i][j][k] = bbox_new(domain->x1 + i*dx, 
                                     domain->x1 + (i+1)*dx,
                                     domain->y1 + j*dx, 
                                     domain->y1 + (j+1)*dy,
                                     domain->z1 + k*dz, 
                                     domain->z1 + (k+1)*dz);
        }
      }
    }
  }


  return grid_data;
}

void amr_grid_data_free(amr_grid_data_t* grid_data)
{
  DECLARE_PATCH_ARRAY(patches, grid_data);
  DECLARE_BBOX_ARRAY(bboxes, grid_data);
  for (int i = 0; i < grid_data->nx; ++i)
  {
    for (int j = 0; j < grid_data->ny; ++j)
    {
      for (int k = 0; i < grid_data->nz; ++k)
      {
        if (patches[i][j][k] != NULL)
          amr_patch_free(patches[i][j][k]);
        if (bboxes[i][j][k] != NULL)
          bboxes[i][j][k] = NULL;
      }
    }
  }
  polymec_free(grid_data);
}

int amr_grid_data_num_patches(amr_grid_data_t* grid_data)
{
  return amr_grid_num_patches(grid_data->grid);
}

amr_grid_t* amr_grid_data_grid(amr_grid_data_t* grid_data)
{
  return grid_data->grid;
}

amr_patch_t* amr_grid_data_patch(amr_grid_data_t* grid_data, int i, int j, int k, bbox_t** bbox)
{
  DECLARE_PATCH_ARRAY(patches, grid_data);
  DECLARE_BBOX_ARRAY(bboxes, grid_data);
  if (bbox != NULL)
    *bbox = bboxes[i][j][k];
  return patches[i][j][k];
}

bool amr_grid_data_next(amr_grid_data_t* grid_data, int* pos, 
                        int* i, int* j, int* k, 
                        amr_patch_t** patch, bbox_t** bbox)
{
  int index = *pos;
  int num_patches = amr_grid_num_patches(grid_data->grid);
  amr_patch_t** patches = grid_data->patches;
  if (index < num_patches)
  {
    *patch = patches[index];
    if (bbox != NULL)
    {
      bbox_t** bboxes = grid_data->bboxes;
      *bbox = bboxes[index];
    }
    while ((patches[index] == NULL) && (index < num_patches))
      ++(*pos);
    return true;
  }
  return false;
}

void amr_grid_data_fill_ghosts(amr_grid_data_t* grid_data)
{
  amr_grid_fill_ghosts(grid_data->grid, grid_data);
}

void amr_grid_data_start_filling_ghosts(amr_grid_data_t* grid_data)
{
  amr_grid_start_filling_ghosts(grid_data->grid, grid_data);
}

void amr_grid_data_finish_filling_ghosts(amr_grid_data_t* grid_data)
{
  amr_grid_finish_filling_ghosts(grid_data->grid, grid_data);
}

