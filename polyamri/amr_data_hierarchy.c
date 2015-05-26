// Copyright (c) 2014-2015, Jeffrey N. Johnson
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
  int num_components;
  ptr_array_t* patch_sets;
};

amr_data_hierarchy_t* amr_data_hierarchy_new(amr_grid_hierarchy_t* grids,
                                             int num_components)
{
  ASSERT(num_components > 0);
  amr_data_hierarchy_t* data = polymec_malloc(sizeof(amr_data_hierarchy_t));
  data->grids = grids;
  data->num_components = num_components;
  data->patch_sets = ptr_array_new();

  int pos = 0;
  amr_grid_t* level;
  while (amr_grid_hierarchy_next_coarsest(grids, &pos, &level))
  {
    amr_patch_set_t* patches = amr_grid_patch_set(level, data->num_components);
    ptr_array_append_with_dtor(data->patch_sets, patches, DTOR(amr_patch_set_free));
  }
  return data;
}

void amr_data_hierarchy_free(amr_data_hierarchy_t* data)
{
  ptr_array_free(data->patch_sets);
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

int amr_data_hierarchy_num_levels(amr_data_hierarchy_t* data)
{
  return amr_grid_hierarchy_num_levels(data->grids);
}

bool amr_data_hierarchy_next_coarsest(amr_data_hierarchy_t* data, int* pos, amr_patch_set_t** patch_set, amr_grid_t** level)
{
  bool found = amr_grid_hierarchy_next_coarsest(data->grids, pos, level);
  if (found)
    *patch_set = data->patch_sets->data[*pos-1];
  return found;
}

bool amr_data_hierarchy_next_finest(amr_data_hierarchy_t* data, int* pos, amr_patch_set_t** patch_set, amr_grid_t** level)
{
  bool found = amr_grid_hierarchy_next_finest(data->grids, pos, level);
  if (found)
    *patch_set = data->patch_sets->data[*pos-1];
  return found;
}

