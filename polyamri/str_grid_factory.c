// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "polyamri/str_grid_factory.h"
#include "polyamri/str_grid_patch_filler.h"
#include "polyamri/grid_to_bbox_coord_mapping.h"

struct str_grid_factory_t 
{
  // MPI communicator.
  MPI_Comm comm;

  // Patch size.
  int patch_nx, patch_ny, patch_nz;
};

str_grid_factory_t* str_grid_factory_new(MPI_Comm comm,
                                         int patch_nx, int patch_ny, int patch_nz)
{
  ASSERT(patch_nx > 0);
  ASSERT(patch_ny > 0);
  ASSERT(patch_nz > 0);

  str_grid_factory_t* factory = GC_MALLOC(sizeof(str_grid_factory_t));
  factory->comm = comm;
  factory->patch_nx = patch_nx;
  factory->patch_ny = patch_ny;
  factory->patch_nz = patch_nz;
  return factory;
}

str_grid_t* str_grid_factory_brick(str_grid_factory_t* factory,
                                   int npx, int npy, int npz,
                                   bbox_t* domain)
{
  ASSERT(npx > 0);
  ASSERT(npy > 0);
  ASSERT(npz > 0);

  // Here we set up the machinery to fill ghost cells in the grid.
  str_grid_patch_filler_t* fill_from_east = copy_str_grid_patch_filler_new(STR_GRID_PATCH_X1_BOUNDARY, STR_GRID_PATCH_X2_BOUNDARY);
  str_grid_patch_filler_t* fill_from_west = copy_str_grid_patch_filler_new(STR_GRID_PATCH_X2_BOUNDARY, STR_GRID_PATCH_X1_BOUNDARY);
  str_grid_patch_filler_t* fill_from_south = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Y2_BOUNDARY, STR_GRID_PATCH_Y1_BOUNDARY);
  str_grid_patch_filler_t* fill_from_north = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Y1_BOUNDARY, STR_GRID_PATCH_Y2_BOUNDARY);
  str_grid_patch_filler_t* fill_from_above = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Z1_BOUNDARY, STR_GRID_PATCH_Z2_BOUNDARY);
  str_grid_patch_filler_t* fill_from_below = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Z2_BOUNDARY, STR_GRID_PATCH_Z1_BOUNDARY);

  // By default, we impose a zero "flux" on all boundaries, on the assumption that the 
  // solution is zero at the boundary.
  str_grid_patch_filler_t* zero_flux_x1 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_X1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_x2 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_X2_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_y1 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Y1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_y2 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Y2_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_z1 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Z1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_z2 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Z2_BOUNDARY);

  // Create the grid.
  str_grid_t* grid = str_grid_new(npx, npy, npz, 
                                  factory->patch_nx, factory->patch_ny, factory->patch_nz, 
                                  false, false, false);

  // FIXME: For now, we assume this is a serial grid and just assign all the patches to rank 0.
  int nprocs;
  MPI_Comm_size(factory->comm, &nprocs);
  ASSERT(nprocs == 1);
  for (int ip = 0; ip < npx; ++ip)
    for (int jp = 0; jp < npy; ++jp)
      for (int kp = 0; kp < npz; ++kp)
        str_grid_insert_patch(grid, ip, jp, kp);

  // Finalize the grid.
  str_grid_finalize(grid);

  // If we are given a bounding box, store it as a mapping in the grid.
  if (domain != NULL)
  {
    coord_mapping_t* mapping = grid_to_bbox_coord_mapping_new(domain);
    str_grid_set_property(grid, "mapping", mapping, NULL);
  }

  return grid;
}

