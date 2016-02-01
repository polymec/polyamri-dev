// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_H
#define POLYAMRI_STR_GRID_H

#include "polyamri/str_grid_patch.h"
#include "polyamri/str_grid_patch_filler.h"

// A structured grid is three-dimensional logically rectangular grid. It 
// consists of a set of uniformly-sized patches. The grid manages
// these patches and their connectivity.
typedef struct str_grid_t str_grid_t;

//------------------------------------------------------------------------
//                          Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct structured grids.
// str_grid_finalize() must be called after a grid has been properly
// constructed.
//------------------------------------------------------------------------

// Creates a new empty grid level defined on the region filling [0,1]x[0,1]x[0,1]
// with nx x ny x nz patches of size px x py x pz. This grid is not associated 
// with any other grids.
str_grid_t* str_grid_new(int nx, int ny, int nz, 
                         int px, int py, int pz,
                         bool periodic_in_x, 
                         bool periodic_in_y, 
                         bool periodic_in_z);

// Destroys the given grid and all of its patches.
void str_grid_free(str_grid_t* grid);

// Inserts a new patch at (i, j, k) in the nx x ny x nz array of patches.
void str_grid_insert_patch(str_grid_t* grid, int i, int j, int k);

// Appends a patch filler to the list of fillers that fill the ghost cells of 
// the patch at the given patch coordinates. This function produces a fatal 
// error if there is no patch at (i, j, k).
void str_grid_append_patch_filler(str_grid_t* grid, int i, int j, int k,
                                  str_grid_patch_filler_t* patch_filler);

// Finalizes the construction process for the grid. This must be called 
// before any of the grid's usage methods (below) are invoked. Should only 
// be called once.
void str_grid_finalize(str_grid_t* grid);

//------------------------------------------------------------------------
//                          Usage methods
//------------------------------------------------------------------------
// The following methods can only be used after a structured grid has been 
// fully constructed and finalized.
//------------------------------------------------------------------------

// Fetches the number of patches in this grid in the x, y, and z directions, 
// placing them in nx, ny, nz.
void str_grid_get_extents(str_grid_t* grid, int* nx, int* ny, int* nz);

// Fetches the number of cells in each patch on this grid in the x, y, and z 
// directions, placing them in pnx, pny, pnz.
void str_grid_get_patch_size(str_grid_t* grid, int* pnx, int* pny, int* pnz);

// Returns the number of patches that can be stored on this grid.
int str_grid_num_patches(str_grid_t* grid);

// Traverses the patches in the grid, returning true and the next (i, j, k) 
// triple if the traversal is incomplete, false otherwise. 
// Set *pos to zero to reset the traversal.
bool str_grid_next_patch(str_grid_t* grid, int* pos, int* i, int* j, int* k);

// Returns true if the grid has a patch at (i, j, k), false if not.
bool str_grid_has_patch(str_grid_t* grid, int i, int j, int k);

// Queries the periodicity of the grid, placing booleans for the 
// x, y, and z periodicity into the given (3-wide) periodicity array.
void str_grid_get_periodicity(str_grid_t* grid, bool* periodicity);

//------------------------------------------------------------------------
//                          Service methods
//------------------------------------------------------------------------
// The following methods provide services to other classes. They can be 
// used directly by folks who Really Understand What They're Doing.
//------------------------------------------------------------------------

// Data can be centered on cells, faces, edges, or nodes, but only cell data
// can have ghost data.
typedef struct str_grid_cell_data_t str_grid_cell_data_t;

// Fills all ghost cells in the patches within the given grid data. This 
// function may involve communication.
void str_grid_fill_ghost_cells(str_grid_t* grid, str_grid_cell_data_t* data);

// Begins an asynchronous ghost-cell-filling operation in the patches within 
// the given grid data, communicating as needed. Returns an integer token that 
// can be used with str_grid_finish_filling_ghosts, below.
int str_grid_start_filling_ghost_cells(str_grid_t* grid, str_grid_cell_data_t* data);

// Concludes an asynchronous ghost-cell-filling operation initiated by 
// a call to str_grid_start_filling_ghosts(), using the token returned by
// that call.
void str_grid_finish_filling_ghost_cells(str_grid_t* grid, int token);

#endif

