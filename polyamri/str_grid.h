// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_H
#define POLYAMRI_STR_GRID_H

#include "core/serializer.h"
#include "polyamri/str_grid_patch.h"

// This type is used to fill ghost cells in patches.
typedef struct str_grid_patch_filler_t str_grid_patch_filler_t;

// A structured grid is a three-dimensional logically rectangular grid. It 
// consists of a set of uniformly-sized patches. The grid manages
// these patches and their connectivity.
typedef struct str_grid_t str_grid_t;

// This type identifies the six different logical grid boundaries.
typedef enum
{
  STR_GRID_X1_BOUNDARY,
  STR_GRID_X2_BOUNDARY,
  STR_GRID_Y1_BOUNDARY,
  STR_GRID_Y2_BOUNDARY,
  STR_GRID_Z1_BOUNDARY,
  STR_GRID_Z2_BOUNDARY
} str_grid_boundary_t;

//------------------------------------------------------------------------
//                          Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct structured grids.
// str_grid_finalize() must be called after a grid has been properly
// constructed.
//------------------------------------------------------------------------

// Creates a new empty grid level defined on the region filling [0,1]x[0,1]x[0,1]
// with npx x npy x npz patches of size nx x ny x nz. This grid is not associated 
// with any other grids.
str_grid_t* str_grid_new(int npx, int npy, int npz, 
                         int nx, int ny, int nz);

// Destroys the given grid and all of its patches.
void str_grid_free(str_grid_t* grid);

// Inserts a new patch at (i, j, k) in the nx x ny x nz array of patches.
void str_grid_insert_patch(str_grid_t* grid, int i, int j, int k);

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
// placing them in npx, npy, npz.
void str_grid_get_extents(str_grid_t* grid, int* npx, int* npy, int* npz);

// Fetches the number of cells in each patch on this grid in the x, y, and z 
// directions, placing them in nx, ny, nz.
void str_grid_get_patch_size(str_grid_t* grid, int* nx, int* ny, int* nz);

// Returns the number of patches that can be stored on this grid.
int str_grid_num_patches(str_grid_t* grid);

// Traverses the patches in the grid, returning true and the next (i, j, k) 
// triple if the traversal is incomplete, false otherwise. 
// Set *pos to zero to reset the traversal.
bool str_grid_next_patch(str_grid_t* grid, int* pos, int* i, int* j, int* k);

// Returns true if the grid has a patch at (i, j, k), false if not.
bool str_grid_has_patch(str_grid_t* grid, int i, int j, int k);

// Associates a named piece of metadata (a "property") with the grid itself.
// This can be used to store information about (for example) how the grid 
// was generated, which can sometimes be useful. A serializer can 
// be given so that any partitioning or repartitioning of the grid can 
// preserve this property on subdomains. If the given property exists on the 
// grid, it is overwritten.
void str_grid_set_property(str_grid_t* grid, 
                           const char* property, 
                           void* data, 
                           serializer_t* serializer);

// Retrieves the given property from the grid, if any. If the 
// property is not found, this returns NULL.
void* str_grid_property(str_grid_t* grid, const char* property);

// Deletes the given property from the grid. This has no effect if the 
// property is not found.
void str_grid_delete_property(str_grid_t* grid, const char* property);

// Allows the traversal of grid properties. Set *pos to 0 to reset the 
// iteration.
bool str_grid_next_property(str_grid_t* grid, int* pos, 
                            char** prop_name, void** prop_data, 
                            serializer_t** prop_serializer);

#endif

