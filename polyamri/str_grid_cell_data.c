// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "polyamri/str_grid_cell_data.h"

struct str_grid_cell_data_t 
{
  str_grid_t* grid;
  int nx, ny, nz, nc, ng;
  int_ptr_unordered_map_t* patches;

  int token;
};

static inline int patch_index(str_grid_cell_data_t* cell_data, int i, int j, int k)
{
  return cell_data->ny*cell_data->nz*i + cell_data->nz*j + k;
}

str_grid_cell_data_t* str_grid_cell_data_new(str_grid_t* grid, 
                                             int num_components, 
                                             int num_ghost_layers)
{
  ASSERT(num_components > 0);
  ASSERT(num_ghost_layers >= 0);

  // Allocate storage.
  int num_patches = str_grid_num_patches(grid);
  size_t patches_size = sizeof(str_grid_patch_t*) * num_patches;
  size_t cell_data_size = sizeof(str_grid_cell_data_t) + patches_size;
  str_grid_cell_data_t* cell_data = polymec_malloc(cell_data_size);
  cell_data->grid = grid;
  cell_data->nc = num_components;
  cell_data->ng = num_ghost_layers;
  str_grid_get_extents(grid, &cell_data->nx, &cell_data->ny, &cell_data->nz);
  cell_data->patches = int_ptr_unordered_map_new();

  // Now populate the patches.
  int px, py, pz;
  str_grid_get_patch_size(grid, &px, &py, &pz);
  int pos = 0, i, j, k;
  while (str_grid_next_patch(grid, &pos, &i, &j, &k))
  {
    int index = patch_index(cell_data, i, j, k);
    int_ptr_unordered_map_insert_with_v_dtor(cell_data->patches, index, 
      str_grid_patch_new(px, py, pz, num_components, num_ghost_layers),
      DTOR(str_grid_patch_free));
  }

  cell_data->token = -1; // No data in flight.
  return cell_data;
}

void str_grid_cell_data_free(str_grid_cell_data_t* cell_data)
{
  int_ptr_unordered_map_free(cell_data->patches);
  polymec_free(cell_data);
}

int str_grid_cell_data_num_components(str_grid_cell_data_t* cell_data)
{
  return cell_data->nc;
}

int str_grid_cell_data_num_ghost_layers(str_grid_cell_data_t* cell_data)
{
  return cell_data->ng;
}

int str_grid_cell_data_num_patches(str_grid_cell_data_t* cell_data)
{
  return str_grid_num_patches(cell_data->grid);
}

str_grid_t* str_grid_cell_data_grid(str_grid_cell_data_t* cell_data)
{
  return cell_data->grid;
}

str_grid_patch_t* str_grid_cell_data_patch(str_grid_cell_data_t* cell_data, int i, int j, int k)
{
  int index = patch_index(cell_data, i, j, k);
  str_grid_patch_t** pp = (str_grid_patch_t**)int_ptr_unordered_map_get(cell_data->patches, index);
  if (pp != NULL)
    return *pp;
  else
    return NULL;
}

bool str_grid_cell_data_next_patch(str_grid_cell_data_t* cell_data, int* pos, 
                                   int* i, int* j, int* k, 
                                   str_grid_patch_t** patch)
{
  bool result = str_grid_next_patch(cell_data->grid, pos, i, j, k);
  if (result)
    *patch = str_grid_cell_data_patch(cell_data, *i, *j, *k);
  return result;
}

void str_grid_cell_data_fill_ghosts(str_grid_cell_data_t* cell_data)
{
  str_grid_fill_ghost_cells(cell_data->grid, cell_data);
}

void str_grid_cell_data_start_filling_ghosts(str_grid_cell_data_t* cell_data)
{
  cell_data->token = str_grid_start_filling_ghost_cells(cell_data->grid, cell_data);
}

void str_grid_cell_data_finish_filling_ghosts(str_grid_cell_data_t* cell_data)
{
  ASSERT(cell_data->token != -1); // No ghost fill in progress?
  str_grid_finish_filling_ghost_cells(cell_data->grid, cell_data->token);
  cell_data->token = -1;
}

