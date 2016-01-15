// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "amr_data_hierarchy.h"

struct amr_data_hierarchy_t 
{
  amr_grid_hierarchy_t* grids;
  int num_components, num_ghosts;
  ptr_array_t* grid_data;
};

amr_data_hierarchy_t* amr_data_hierarchy_new(amr_grid_hierarchy_t* grids,
                                             int num_components,
                                             int num_ghosts)
{
  ASSERT(num_components > 0);
  amr_data_hierarchy_t* data = polymec_malloc(sizeof(amr_data_hierarchy_t));
  data->grids = grids;
  data->num_components = num_components;
  data->num_ghosts = num_ghosts;
  data->grid_data = ptr_array_new();

  int pos = 0;
  amr_grid_t* level;
  while (amr_grid_hierarchy_next_coarsest(grids, &pos, &level))
  {
    amr_grid_data_t* grid_data = amr_grid_data_new(level, AMR_GRID_CELL, data->num_components, data->num_ghosts);
    ptr_array_append_with_dtor(data->grid_data, grid_data, DTOR(amr_grid_data_free));
  }
  return data;
}

void amr_data_hierarchy_free(amr_data_hierarchy_t* data)
{
  ptr_array_free(data->grid_data);
  polymec_free(data);
}

amr_grid_hierarchy_t* amr_data_hierarchy_grids(amr_data_hierarchy_t* data)
{
  return data->grids;
}

int amr_data_hierarchy_num_components(amr_data_hierarchy_t* data)
{
  return data->num_components;
}

int amr_data_hierarchy_num_ghosts(amr_data_hierarchy_t* data)
{
  return data->num_ghosts;
}

int amr_data_hierarchy_num_levels(amr_data_hierarchy_t* data)
{
  return amr_grid_hierarchy_num_levels(data->grids);
}

bool amr_data_hierarchy_next_coarsest(amr_data_hierarchy_t* data, int* pos, amr_grid_data_t** grid_data)
{
  if (*pos >= data->grid_data->size)
    return false;
  *grid_data = data->grid_data->data[*pos];
  ++(*pos);
  return true;
}

bool amr_data_hierarchy_next_finest(amr_data_hierarchy_t* data, int* pos, amr_grid_data_t** grid_data)
{
  if (*pos >= data->grid_data->size)
    return false;
  *grid_data = data->grid_data->data[data->grid_data->size - *pos - 1];
  ++(*pos);
  return true;
}

