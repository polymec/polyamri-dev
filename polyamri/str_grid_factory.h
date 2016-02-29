// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_FACTORY_H
#define POLYAMRI_STR_GRID_FACTORY_H

#include "polyamri/str_grid.h"

// The structured grid factory creates str_grid objects for use in simulations.
// Objects of this type are garbage-collected.
typedef struct str_grid_factory_t str_grid_factory_t;

// Creates a new structured grid factory that creates grids on the given MPI
// communicator, comprising patches with the given number of cells in x, y, and z
str_grid_factory_t* str_grid_factory_new(MPI_Comm comm, 
                                         int patch_nx, int patch_ny, int patch_nz);

//------------------------------------------------------------------------
//                            Factory Methods
//------------------------------------------------------------------------
// The following methods produce single-block structured grids that are 
// understood to be mapped from a "logical" domain [0,1]x[0,1]x[0,1] to 
// some sort of "physical" domain. The coordinate mapping for the grid is 
// stored in its "mapping" property.
//------------------------------------------------------------------------

// Creates a "brick"-shaped structured grid consisting of npx x npy x npz 
// patches, filling the rectangular domain defined by the given bounding box.
str_grid_t* str_grid_factory_brick(str_grid_factory_t* factory, 
                                   int npx, int npy, int npz,
                                   bbox_t* domain);

#endif

