// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "polyamri/str_grid_cell_filler_factory.h"

struct str_grid_cell_filler_factory_t 
{
  MPI_Comm comm;
};

str_grid_cell_filler_factory_t* str_grid_cell_filler_factory_new(MPI_Comm comm)
{
  str_grid_cell_filler_factory_t* factory = GC_MALLOC(sizeof(str_grid_cell_filler_factory_t));
  factory->comm = comm;
  return factory;
}

str_grid_cell_filler_t* str_grid_cell_filler_factory_ghost_filler(str_grid_cell_filler_factory_t* factory, 
                                                                  str_grid_t* grid,
                                                                  str_grid_patch_filler_t* x1_boundary_filler,
                                                                  str_grid_patch_filler_t* x2_boundary_filler,
                                                                  str_grid_patch_filler_t* y1_boundary_filler,
                                                                  str_grid_patch_filler_t* y2_boundary_filler,
                                                                  str_grid_patch_filler_t* z1_boundary_filler,
                                                                  str_grid_patch_filler_t* z2_boundary_filler)
{
  // FIXME: For now, we assume this is a serial configuration.
  int nprocs;
  MPI_Comm_size(factory->comm, &nprocs);
  ASSERT(nprocs == 1);

  // Get the number of patches in x, y, z.
  int npx, npy, npz;
  str_grid_get_extents(grid, &npx, &npy, &npz);

  // Create a new cell filler and insert patch fillers into it.
  str_grid_cell_filler_t* filler = str_grid_cell_filler_new(grid);
  str_grid_patch_filler_t* fill_from_east = copy_str_grid_patch_filler_new(STR_GRID_PATCH_X1_BOUNDARY, STR_GRID_PATCH_X2_BOUNDARY);
  str_grid_patch_filler_t* fill_from_west = copy_str_grid_patch_filler_new(STR_GRID_PATCH_X2_BOUNDARY, STR_GRID_PATCH_X1_BOUNDARY);
  str_grid_patch_filler_t* fill_from_south = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Y2_BOUNDARY, STR_GRID_PATCH_Y1_BOUNDARY);
  str_grid_patch_filler_t* fill_from_north = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Y1_BOUNDARY, STR_GRID_PATCH_Y2_BOUNDARY);
  str_grid_patch_filler_t* fill_from_above = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Z1_BOUNDARY, STR_GRID_PATCH_Z2_BOUNDARY);
  str_grid_patch_filler_t* fill_from_below = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Z2_BOUNDARY, STR_GRID_PATCH_Z1_BOUNDARY);

  int pos = 0, ip, jp, kp;
  while (str_grid_next_patch(grid, &pos, &ip, &jp, &kp))
  {
    if (ip > 0)
      str_grid_cell_filler_insert(filler, ip, jp, kp, fill_from_west);
    else
      str_grid_cell_filler_insert(filler, ip, jp, kp, x1_boundary_filler);
    if (ip < npx-1)
      str_grid_cell_filler_insert(filler, ip, jp, kp, fill_from_east);
    else
      str_grid_cell_filler_insert(filler, ip, jp, kp, x2_boundary_filler);
    if (jp > 0)
      str_grid_cell_filler_insert(filler, ip, jp, kp, fill_from_south);
    else
      str_grid_cell_filler_insert(filler, ip, jp, kp, y1_boundary_filler);
    if (jp < npy-1)
      str_grid_cell_filler_insert(filler, ip, jp, kp, fill_from_north);
    else
      str_grid_cell_filler_insert(filler, ip, jp, kp, y2_boundary_filler);
    if (kp > 0)
      str_grid_cell_filler_insert(filler, ip, jp, kp, fill_from_below);
    else
      str_grid_cell_filler_insert(filler, ip, jp, kp, z1_boundary_filler);
    if (kp < npz-1)
      str_grid_cell_filler_insert(filler, ip, jp, kp, fill_from_above);
    else
      str_grid_cell_filler_insert(filler, ip, jp, kp, z2_boundary_filler);
  }

  return filler;
}

