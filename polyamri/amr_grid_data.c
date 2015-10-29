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
  amr_grid_data_centering_t centering;
  int nx, ny, nz, nc, ng;
  void* patches;

  int token;
};

// Multi-dimensional array accessors for patches.
#define DECLARE_PATCH_ARRAY(array, grid_data) \
amr_patch_t* (*array)[grid_data->nx][grid_data->ny] = grid_data->patches

amr_grid_data_t* amr_grid_data_new(amr_grid_t* grid, 
                                   amr_grid_data_centering_t centering,
                                   int num_components,
                                   int num_ghosts)
{
  ASSERT(num_components > 0);
  ASSERT(num_ghosts >= 0);

  // Allocate storage.
  int num_local_patches = amr_grid_num_local_patches(grid);
  size_t patches_size = sizeof(amr_patch_t*) * num_local_patches;
  size_t grid_data_size = sizeof(amr_grid_data_t) + patches_size;
  amr_grid_data_t* grid_data = polymec_malloc(grid_data_size);
  grid_data->grid = grid;
  grid_data->centering = centering;
  grid_data->nc = num_components;
  grid_data->ng = num_ghosts;
  amr_grid_get_extents(grid, &grid_data->nx, &grid_data->ny, &grid_data->nz);
  grid_data->patches = (char*)grid_data + sizeof(amr_grid_data_t);
  memset(grid_data->patches, 0, patches_size);

  // Now populate the patches.
  int px, py, pz;
  amr_grid_get_patch_size(grid, &px, &py, &pz);
  DECLARE_PATCH_ARRAY(patches, grid_data);
  int pos = 0, i, j, k;
  if (centering == AMR_GRID_CELL)
  {
    while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k))
      patches[i][j][k] = amr_patch_new(px, py, pz, num_components, num_ghosts);
  }
  else if (centering == AMR_GRID_X_FACE)
  {
    while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k))
      patches[i][j][k] = amr_patch_new(px+1, py, pz, num_components, num_ghosts);
  }
  else if (centering == AMR_GRID_Y_FACE)
  {
    while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k))
      patches[i][j][k] = amr_patch_new(px, py+1, pz, num_components, num_ghosts);
  }
  else if (centering == AMR_GRID_Z_FACE)
  {
    while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k))
      patches[i][j][k] = amr_patch_new(px, py, pz+1, num_components, num_ghosts);
  }
  else if (centering == AMR_GRID_X_EDGE)
  {
    while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k))
      patches[i][j][k] = amr_patch_new(px, py+1, pz+1, num_components, num_ghosts);
  }
  else if (centering == AMR_GRID_Y_EDGE)
  {
    while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k))
      patches[i][j][k] = amr_patch_new(px+1, py, pz+1, num_components, num_ghosts);
  }
  else if (centering == AMR_GRID_Z_EDGE)
  {
    while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k))
      patches[i][j][k] = amr_patch_new(px+1, py+1, pz, num_components, num_ghosts);
  }
  else
  {
    ASSERT(centering == AMR_GRID_NODE);
    while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k))
      patches[i][j][k] = amr_patch_new(px+1, py+1, pz+1, num_components, num_ghosts);
  }

  grid_data->token = -1; // No data in flight.
  return grid_data;
}

void amr_grid_data_free(amr_grid_data_t* grid_data)
{
  int pos = 0, i, j, k;
  amr_patch_t* patch;
  while (amr_grid_data_next_local_patch(grid_data, &pos, &i, &j, &k, &patch))
    amr_patch_free(patch);
  polymec_free(grid_data);
}

amr_grid_data_centering_t amr_grid_data_centering(amr_grid_data_t* grid_data)
{
  return grid_data->centering;
}

int amr_grid_data_num_ghosts(amr_grid_data_t* grid_data)
{
  return grid_data->ng;
}

int amr_grid_data_num_components(amr_grid_data_t* grid_data)
{
  return grid_data->nc;
}

int amr_grid_data_num_local_patches(amr_grid_data_t* grid_data)
{
  return amr_grid_num_local_patches(grid_data->grid);
}

amr_grid_t* amr_grid_data_grid(amr_grid_data_t* grid_data)
{
  return grid_data->grid;
}

amr_patch_t* amr_grid_data_patch(amr_grid_data_t* grid_data, int i, int j, int k)
{
  DECLARE_PATCH_ARRAY(patches, grid_data);
  return patches[i][j][k];
}

bool amr_grid_data_next_local_patch(amr_grid_data_t* grid_data, int* pos, 
                                    int* i, int* j, int* k, 
                                    amr_patch_t** patch)
{
  bool result = amr_grid_next_local_patch(grid_data->grid, pos, i, j, k);
  if (result)
  {
    DECLARE_PATCH_ARRAY(patches, grid_data);
    *patch = patches[*i][*j][*k];
  }
  return result;
}

void amr_grid_data_fill_ghosts(amr_grid_data_t* grid_data)
{
  amr_grid_fill_ghosts(grid_data->grid, grid_data);
}

void amr_grid_data_start_filling_ghosts(amr_grid_data_t* grid_data)
{
  ASSERT(grid_data->token == -1);
  grid_data->token = amr_grid_start_filling_ghosts(grid_data->grid, grid_data);
}

void amr_grid_data_finish_filling_ghosts(amr_grid_data_t* grid_data)
{
  if (grid_data->token != -1)
  {
    amr_grid_finish_filling_ghosts(grid_data->grid, grid_data->token);
    grid_data->token = -1;
  }
}

