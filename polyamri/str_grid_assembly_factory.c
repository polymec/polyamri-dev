// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "polyamri/str_grid_assembly_factory.h"

struct str_grid_assembly_factory_t 
{
  // MPI communicator.
  MPI_Comm comm;

  // Patch size.
  int patch_nx, patch_ny, patch_nz;
};

str_grid_assembly_factory_t* str_grid_assembly_factory_new(MPI_Comm comm, 
                                                           int patch_nx, 
                                                           int patch_ny, 
                                                           int patch_nz)
{
  ASSERT(patch_nx > 0);
  ASSERT(patch_ny > 0);
  ASSERT(patch_nz > 0);

  str_grid_assembly_factory_t* factory = GC_MALLOC(sizeof(str_grid_assembly_factory_t));
  factory->comm = comm;
  factory->patch_nx = patch_nx;
  factory->patch_ny = patch_ny;
  factory->patch_nz = patch_nz;
  return factory;
}

str_grid_assembly_t* str_grid_assembly_factory_cubed_sphere(str_grid_assembly_factory_t* factory, 
                                                            cubed_sphere_mapping_t mapping_type,
                                                            int npx, int npy, int npz,
                                                            real_t R)
{
  str_grid_assembly_t* assembly = str_grid_assembly_new();

  // Insert 6 identical blocks into the assembly.
  const char* block_names[6] = {"equatorial_1", "equatorial_2", 
                                "equatorial_3", "equatorial_4", 
                                "north_pole", "south_pole"};
  for (int b = 0; b < 6; ++b)
  {
    str_grid_t* block = str_grid_new(npx, npy, npz, factory->patch_nx, factory->patch_ny, factory->patch_nz);
    str_grid_assembly_add_block(assembly, block_names[b], block);
  }

  // Connect the equatorial blocks to one another.
  for (int b = 0; b < 4; ++b)
  {
    str_grid_assembly_connect(assembly, 
                              block_names[b], STR_GRID_X2_BOUNDARY, 
                              block_names[(b+1)%4], STR_GRID_X1_BOUNDARY);
  }

  // Connect the equatorial blocks to the north pole.
  str_grid_assembly_connect(assembly, 
                            block_names[0], STR_GRID_Y2_BOUNDARY, 
                            block_names[5], STR_GRID_Y1_BOUNDARY);
  str_grid_assembly_connect(assembly, 
                            block_names[1], STR_GRID_Y2_BOUNDARY, 
                            block_names[5], STR_GRID_X2_BOUNDARY);
  str_grid_assembly_connect(assembly, 
                            block_names[2], STR_GRID_Y2_BOUNDARY, 
                            block_names[5], STR_GRID_Y2_BOUNDARY);
  str_grid_assembly_connect(assembly, 
                            block_names[3], STR_GRID_Y2_BOUNDARY, 
                            block_names[5], STR_GRID_X1_BOUNDARY);

  // Connect the equatorial blocks to the south pole.
  str_grid_assembly_connect(assembly, 
                            block_names[0], STR_GRID_Y1_BOUNDARY, 
                            block_names[6], STR_GRID_Y2_BOUNDARY);
  str_grid_assembly_connect(assembly, 
                            block_names[1], STR_GRID_Y1_BOUNDARY, 
                            block_names[6], STR_GRID_X2_BOUNDARY);
  str_grid_assembly_connect(assembly, 
                            block_names[2], STR_GRID_Y1_BOUNDARY, 
                            block_names[6], STR_GRID_Y1_BOUNDARY);
  str_grid_assembly_connect(assembly, 
                            block_names[3], STR_GRID_Y1_BOUNDARY, 
                            block_names[6], STR_GRID_X1_BOUNDARY);

  // Set up mappings on each block.

  return assembly;
}

