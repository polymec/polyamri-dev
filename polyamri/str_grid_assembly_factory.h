// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_ASSEMBLY_FACTORY_H
#define POLYAMRI_STR_GRID_ASSEMBLY_FACTORY_H

#include "polyamri/str_grid_assembly.h"

// The structured grid assembly factory creates str_grid_assembly objects for 
// use in simulations. Objects of this type are garbage-collected.
typedef struct str_grid_assembly_factory_t str_grid_assembly_factory_t;

// Creates a new structured grid assembly factory that creates grids on the 
// given MPI communicator, comprising patches with the given number of cells 
// in x, y, and z.
str_grid_assembly_factory_t* str_grid_assembly_factory_new(MPI_Comm comm, 
                                                           int patch_nx, 
                                                           int patch_ny, 
                                                           int patch_nz);

//------------------------------------------------------------------------
//                            Factory Methods
//------------------------------------------------------------------------
// The following methods produce assemble of multi-block structured grids
// (blocks), each of which is understood to be mapped from a "logical" domain 
// [0,1]x[0,1]x[0,1] to some sort of "physical" domain. The coordinate 
// mapping for each block is stored in its "mapping" property.
//------------------------------------------------------------------------

// This type identifies different cubed-sphere coordinate mappings.
typedef enum
{
  CUBED_SPHERE_CONFORMAL,
  CUBED_SPHERE_GNOMONIC
} cubed_sphere_mapping_t;

// Creates a cubed sphere assembly consisting of 6 blocks, 4 spanning
// the equator of a sphere, and 2 polar blocks. The cells in the blocks are 
// mapped according to the given mapping type. For each block, the number of 
// patches in logical coordinates are given by npx, npy, and npz. The radius 
// of the sphere in physical space is R.
str_grid_assembly_t* str_grid_assembly_factory_cubed_sphere(str_grid_assembly_factory_t* factory, 
                                                            cubed_sphere_mapping_t mapping_type,
                                                            int npx, int npy, int npz,
                                                            real_t R);

#endif

