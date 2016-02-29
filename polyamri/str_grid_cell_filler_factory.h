// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_CELL_FILLER_FACTORY_H
#define POLYAMRI_STR_GRID_CELL_FILLER_FACTORY_H

#include "polyamri/str_grid_cell_filler.h"

// A structured grid cell filler factory creates grid cell filler objects for 
// use with cell-centered data on grids. Objects of this type are garbage-collected.
typedef struct str_grid_cell_filler_factory_t str_grid_cell_filler_factory_t;

// Creates a new grid cell filler factory for the given MPI communicator.
str_grid_cell_filler_factory_t* str_grid_cell_filler_factory_new(MPI_Comm comm);

// Creates a grid cell filler object for cell-centered fields on the given grid, 
// defining the given boundary conditions for the given boundaries on that grid 
// and performing local/remote copies to fill ghost cells in the patches of the 
// grid based on the underlying topology of those patches.
str_grid_cell_filler_t* str_grid_cell_filler_factory_ghost_filler(str_grid_cell_filler_factory_t* factory, 
                                                                  str_grid_t* grid,
                                                                  str_grid_patch_filler_t* x1_boundary_filler,
                                                                  str_grid_patch_filler_t* x2_boundary_filler,
                                                                  str_grid_patch_filler_t* y1_boundary_filler,
                                                                  str_grid_patch_filler_t* y2_boundary_filler,
                                                                  str_grid_patch_filler_t* z1_boundary_filler,
                                                                  str_grid_patch_filler_t* z2_boundary_filler);

#endif

