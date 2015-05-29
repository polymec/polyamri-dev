// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_GRID_ASSEMBLY_H
#define POLYAMRI_AMR_GRID_ASSEMBLY_H

#include "polyamri/amr_grid_hierarchy.h"

// A grid assembly is a set of grid hierarchies representing a domain in 3D 
// that can be topologically represented by multiple rectangular blocks.
typedef struct amr_grid_assembly_t amr_grid_assembly_t;

// Creates a new empty grid assembly.
amr_grid_assembly_t* amr_grid_assembly_new();

// Destroys the given grid assembly and all of its blocks.
void amr_grid_assembly_free(amr_grid_assembly_t* assembly);

// Returns the current number of blocks (grid hierarchies) within this assembly.
int amr_grid_assembly_num_blocks(amr_grid_assembly_t* assembly);

// Adds the given block (grid_hierarchy) to this assembly, associating the 
// block with the given name. The assembly assumes responsibility for the 
// block's resources.
void amr_grid_assembly_add_block(amr_grid_assembly_t* assembly,
                                 const char* block_name,       
                                 amr_grid_hierarchy_t* block);

// Returns the block with the given name, or NULL if no such block can be 
// found in the assembly.
amr_grid_hierarchy_t* amr_grid_assembly_block(amr_grid_assembly_t* assembly,
                                              const char* block_name);

// Traverses the blocks in the grid assembly.
// Set *pos to 0 to reset the iteration.
bool amr_grid_assembly_next_block(amr_grid_assembly_t* assembly, 
                                  int* pos, 
                                  char** block_name,
                                  amr_grid_hierarchy_t** block);

#endif

