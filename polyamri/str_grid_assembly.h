// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_ASSEMBLY_H
#define POLYAMRI_STR_GRID_ASSEMBLY_H

#include "polyamri/str_grid.h"

// A structured grid assembly is a set of structured grids representing a 
// domain in 3D that can be topologically represented by multiple rectangular 
// blocks.
typedef struct str_grid_assembly_t str_grid_assembly_t;

// Creates a new empty structured grid assembly.
str_grid_assembly_t* str_grid_assembly_new();

// Destroys the given grid assembly and all of its blocks.
void str_grid_assembly_free(str_grid_assembly_t* assembly);

// Returns the current number of blocks (structured grid) within this assembly.
int str_grid_assembly_num_blocks(str_grid_assembly_t* assembly);

// Adds the given block (structured grid) to this assembly, associating the 
// block with the given name. The assembly assumes responsibility for the 
// block's resources. The block is not connected to other blocks when added.
void str_grid_assembly_add_block(str_grid_assembly_t* assembly,
                                 const char* block_name,       
                                 str_grid_t* block);

// Connects the given boundaries of blocks with the given names. Blocks with
// these names are assumed to be present in the assembly.
void str_grid_assembly_connect(str_grid_assembly_t* assembly,
                               const char* block1_name,
                               str_grid_boundary_t block1_boundary,
                               const char* block2_name,
                               str_grid_boundary_t block2_boundary);

// Returns the block with the given name, or NULL if no such block can be 
// found in the assembly.
str_grid_t* str_grid_assembly_block(str_grid_assembly_t* assembly,
                                    const char* block_name);

// Traverses the blocks in the grid assembly. The neighbors of each block
// are stored in their respective pointers (which are set to NULL if the block
// has no neighbor on a given boundary). Set *pos to 0 to reset the iteration.
bool str_grid_assembly_next_block(str_grid_assembly_t* assembly, 
                                  int* pos, 
                                  char** block_name,
                                  str_grid_t** block,
                                  str_grid_t** block_x1_neighbor,
                                  str_grid_t** block_x2_neighbor,
                                  str_grid_t** block_y1_neighbor,
                                  str_grid_t** block_y2_neighbor,
                                  str_grid_t** block_z1_neighbor,
                                  str_grid_t** block_z2_neighbor);

#endif

