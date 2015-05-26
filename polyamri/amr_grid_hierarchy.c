// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "polyamri/amr_grid_hierarchy.h"

struct amr_grid_hierarchy_t 
{
  bbox_t domain;
  int nx, ny, nz, px, py, pz;
  int num_ghosts, ref_ratio;
  bool x_periodic, y_periodic, z_periodic;
  ptr_array_t* levels;
};

amr_grid_hierarchy_t* amr_grid_hierarchy_new(bbox_t* domain, 
                                             int nx, int ny, int nz, 
                                             int px, int py, int pz,
                                             int num_ghosts, int ref_ratio,
                                             bool periodic_in_x, bool periodic_in_y, bool periodic_in_z)
{
  ASSERT(domain->x2 > domain->x1);
  ASSERT(domain->y2 > domain->y1);
  ASSERT(domain->z2 > domain->z1);
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  ASSERT(px > 0);
  ASSERT(py > 0);
  ASSERT(pz > 0);
  ASSERT(num_ghosts >= 0);
  ASSERT((ref_ratio % 2) == 0); // FIXME: Not good enough!

  amr_grid_hierarchy_t* hierarchy = malloc(sizeof(amr_grid_hierarchy_t));
  hierarchy->domain = *domain;
  hierarchy->nx = nx, hierarchy->ny = ny, hierarchy->nz = nz;
  hierarchy->px = px, hierarchy->py = py, hierarchy->pz = pz;
  hierarchy->num_ghosts = num_ghosts;
  hierarchy->ref_ratio = ref_ratio;
  hierarchy->x_periodic = periodic_in_x;
  hierarchy->y_periodic = periodic_in_y;
  hierarchy->z_periodic = periodic_in_z;
  hierarchy->levels = ptr_array_new();
  return hierarchy;
}

void amr_grid_hierarchy_free(amr_grid_hierarchy_t* hierarchy)
{
  ptr_array_free(hierarchy->levels);
  free(hierarchy);
}

bbox_t* amr_grid_hierarchy_domain(amr_grid_hierarchy_t* hierarchy)
{
  return &hierarchy->domain;
}

int amr_grid_hierarchy_num_levels(amr_grid_hierarchy_t* hierarchy)
{
  return hierarchy->levels->size;
}

int amr_grid_hierarchy_num_ghosts(amr_grid_hierarchy_t* hierarchy)
{
  return hierarchy->num_ghosts;
}

int amr_grid_hierarchy_ref_ratio(amr_grid_hierarchy_t* hierarchy)
{
  return hierarchy->ref_ratio;
}

void amr_grid_hierarchy_get_periodicity(amr_grid_hierarchy_t* hierarchy, bool* periodicity)
{
  ASSERT(periodicity != NULL);
  periodicity[0] = hierarchy->x_periodic;
  periodicity[1] = hierarchy->y_periodic;
  periodicity[2] = hierarchy->z_periodic;
}


amr_grid_t* amr_grid_hierarchy_add_level(amr_grid_hierarchy_t* hierarchy)
{
  int num_levels = hierarchy->levels->size;
  int nx = hierarchy->nx, ny = hierarchy->ny, nz = hierarchy->nz;
  int px = hierarchy->px, py = hierarchy->py, pz = hierarchy->pz;
  int ref_ratio = hierarchy->ref_ratio;
  for (int l = 0; l < num_levels; ++l)
  {
    nx *= ref_ratio;
    ny *= ref_ratio;
    nz *= ref_ratio;
  }
  amr_grid_t* new_level = amr_grid_new(&hierarchy->domain, nx, ny, nz,
                                       px, py, pz, hierarchy->num_ghosts,
                                       hierarchy->x_periodic, 
                                       hierarchy->y_periodic, 
                                       hierarchy->z_periodic);
  ptr_array_append_with_dtor(hierarchy->levels, new_level, DTOR(amr_grid_free));

  // Hook this grid level up to the next coarser one.
  if (num_levels > 0)
  {
    amr_grid_t* coarser = hierarchy->levels->data[num_levels-1];
    amr_grid_associate_finer(coarser, new_level, ref_ratio);
    amr_grid_associate_coarser(new_level, coarser, ref_ratio);
  }

  return new_level;
}

bool amr_grid_hierarchy_next_coarsest(amr_grid_hierarchy_t* hierarchy, int* pos, amr_grid_t** level)
{
  if (*pos >= hierarchy->levels->size)
    return false;
  *level = hierarchy->levels->data[*pos];
  ++(*pos);
  return true;
}

bool amr_grid_hierarchy_next_finest(amr_grid_hierarchy_t* hierarchy, int* pos, amr_grid_t** level)
{
  if (*pos >= hierarchy->levels->size)
    return false;
  *level = hierarchy->levels->data[hierarchy->levels->size - *pos - 1];
  ++(*pos);
  return true;
}

