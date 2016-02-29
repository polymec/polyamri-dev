// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_CELL_FILLER_H
#define POLYAMRI_STR_GRID_CELL_FILLER_H

#include "polyamri/str_grid.h"
#include "polyamri/str_grid_patch_filler.h"
#include "polyamri/str_grid_cell_data.h"

// A structured grid cell filler is an object that fills cells in a 
// structured_grid according to a given algorithm. Cell fillers can be used 
// to fill ghost cells for cell-centered data between patches on a grid, or 
// to implement boundary conditions on grids, or even to implement source terms
// in patches. The options are many. Because different cell data can use these 
// cell fillers in different ways, many cell fillers can be used in the 
// implementation of a single numerical model.
// 
// Each str_grid_cell_filler manages and coordinates several 
// str_grid_patch_filler objects to fill the cells in a str_grid. These 
// patch fillers must be explicitly added to a cell filler during its 
// construction.
typedef struct str_grid_cell_filler_t str_grid_cell_filler_t;

// Creates a new empty grid cell filler that operates on data for the given 
// structured grid.
str_grid_cell_filler_t* str_grid_cell_filler_new(str_grid_t* grid);

// Destroys the given cell filler and all of its patch fillers.
void str_grid_cell_filler_free(str_grid_cell_filler_t* cell_filler);

// Inserts a patch filler that will fill cells within the the patch at the 
// given patch coordinates. This function produces a fatal error if there is 
// no patch at (i, j, k).
void str_grid_cell_filler_insert(str_grid_cell_filler_t* cell_filler, 
                                 int i, int j, int k,
                                 str_grid_patch_filler_t* patch_filler);

// Returns an internal pointer to the underlying str_grid.
str_grid_t* str_grid_cell_filler_grid(str_grid_cell_filler_t* cell_filler);

// Fills the cells in the patches within the given grid data using the logic
// in this cell filler. This function may involve communication.
void str_grid_cell_filler_fill(str_grid_cell_filler_t* cell_filler, str_grid_cell_data_t* data);

// Begins an asynchronous cell-filling operation in the patches within 
// the given grid data, communicating as needed. Returns an integer token that 
// can be used with str_grid_cell_filler_finish, below.
int str_grid_cell_filler_start(str_grid_cell_filler_t* cell_filler, str_grid_cell_data_t* data);

// Concludes an asynchronous cell-filling operation initiated by a call to 
// str_grid_cell_filler_start(), using the token returned by that call.
void str_grid_cell_filler_finish(str_grid_cell_filler_t* cell_filler, int token);

#endif

