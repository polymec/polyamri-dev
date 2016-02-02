// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_CELL_DATA_H
#define POLYAMRI_STR_GRID_CELL_DATA_H

#include "core/point.h"
#include "polyamri/str_grid.h"

// A structured grid cell data object is a collection of cell-centered patches 
// that are associated with a structured grid in 3D space.
typedef struct str_grid_cell_data_t str_grid_cell_data_t;

// Creates a str_grid_cell_data object associated with the given grid, with 
// the given number of components and ghost layers.
str_grid_cell_data_t* str_grid_cell_data_new(str_grid_t* grid, 
                                             int num_components, 
                                             int num_ghost_layers);

// Frees the given str_grid_cell_data object.
void str_grid_cell_data_free(str_grid_cell_data_t* cell_data);

// Returns the number of components in the str_grid_cell_data object.
int str_grid_cell_data_num_components(str_grid_cell_data_t* cell_data);

// Returns the number of ghost layers in the str_grid_cell_data object.
int str_grid_cell_data_num_ghost_layers(str_grid_cell_data_t* cell_data);

// Returns the number of patches in the str_grid_cell_data object.
int str_grid_cell_data_num_patches(str_grid_cell_data_t* cell_data);

// Returns an internal pointer to the given object's underlying str_grid.
str_grid_t* str_grid_cell_data_grid(str_grid_cell_data_t* cell_data);

// Given a tuple (i, j, k) identifying a patch in the underlying str_grid,
// returns a patch containing data, or NULL if the patch is not present 
// in the str_grid. 
str_grid_patch_t* str_grid_cell_data_patch(str_grid_cell_data_t* cell_data, int i, int j, int k);

// Traverses the grid data, returning true if a patch was found and false if not.
// Set *pos to 0 to reset the traversal. patch is set to the cell patch.
bool str_grid_cell_data_next_patch(str_grid_cell_data_t* cell_data, int* pos, 
                                   int* i, int* j, int* k, 
                                   str_grid_patch_t** patch);

// Fills all ghost values in the str_grid_cell_data's patches.
void str_grid_cell_data_fill_ghosts(str_grid_cell_data_t* cell_data);

// Begins an asynchronous ghost-filling operation.
void str_grid_cell_data_start_filling_ghosts(str_grid_cell_data_t* cell_data);

// Concludes an asynchronous ghost-filling operation initiated by 
// a call to str_grid_cell_data_start_filling_ghosts().
void str_grid_cell_data_finish_filling_ghosts(str_grid_cell_data_t* cell_data);

#endif

