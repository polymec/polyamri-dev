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

// Creates a str_grid_cell_data object associated with the given grid, with 
// patches whose data is aliased to data in the given patch buffer, which is 
// a serialized into a sequential buffer. This object does not own the data 
// in the buffer--it only accesses it. It is up to the caller to ensure that 
// the lifetime of the buffer exceeds that of the resulting str_grid_cell_data
// object.
str_grid_cell_data_t* str_grid_cell_data_alias(str_grid_t* grid, 
                                               int num_components, 
                                               int num_ghost_layers,
                                               real_t* buffer);

// Frees the given str_grid_cell_data object.
void str_grid_cell_data_free(str_grid_cell_data_t* cell_data);

// Returns the number of components in the str_grid_cell_data object.
int str_grid_cell_data_num_components(str_grid_cell_data_t* cell_data);

// Returns the number of ghost layers in the str_grid_cell_data object.
int str_grid_cell_data_num_ghost_layers(str_grid_cell_data_t* cell_data);

// Returns the number of patches in the str_grid_cell_data object.
int str_grid_cell_data_num_patches(str_grid_cell_data_t* cell_data);

// Returns the total number of cells in the str_grid_cell_data object. 
// If include_ghosts is set to true, ghost cells are included -- otherwise 
// only interior cells are included.
int str_grid_cell_data_num_cells(str_grid_cell_data_t* cell_data, bool include_ghosts);

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

// Copies all of the cell data to a contiguous buffer. The buffer is assumed 
// to be large enough to hold all of the cells (including ghost cells).
void str_grid_cell_data_pack(str_grid_cell_data_t* cell_data, real_t* buffer);

// Copies all of the cell data from a contigous buffer into this object.
void str_grid_cell_data_unpack(str_grid_cell_data_t* cell_data, real_t* buffer);

#endif

