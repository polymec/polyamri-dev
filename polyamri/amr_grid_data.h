// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_GRID_DATA_H
#define POLYAMRI_AMR_GRID_DATA_H

#include "core/point.h"
#include "polyamri/amr_grid.h"
#include "polyamri/amr_patch.h"

// An AMR grid data object is a collection of patches that are associated 
// with an amr_grid in 3D space.
typedef struct amr_grid_data_t amr_grid_data_t;

// Creates an amr_grid_data object with the given number of components, 
// associated with the given amr_grid.
amr_grid_data_t* amr_grid_data_new(amr_grid_t* grid, int num_components);

// Frees the given amr_grid_data object.
void amr_grid_data_free(amr_grid_data_t* grid_data);

// Returns the number of (local) patches in the amr_grid_data object.
int amr_grid_data_num_patches(amr_grid_data_t* grid_data);

// Returns an internal pointer to the given object's underlying amr_grid.
amr_grid_t* amr_grid_data_grid(amr_grid_data_t* grid_data);

// Given a tuple (i, j, k) identifying a patch in the underlying amr_grid,
// returns a patch containing data, or NULL if the patch is not present 
// locally on the amr_grid. If bbox is non-NULL, it is set to the bounding 
// box for the patch.
amr_patch_t* amr_grid_data_patch(amr_grid_data_t* grid_data, int i, int j, int k, bbox_t** bbox);

// Traverses the grid data, returning true if a patch was found and false if not.
// Set *pos to 0 to reset the traversal. Pointers to the patch and its 
// bounding box are returned in the given locations.
bool amr_grid_data_next(amr_grid_data_t* grid_data, int* pos, 
                        int* i, int* j, int* k, 
                        amr_patch_t** patch, bbox_t** bbox);


// Fills all ghost cells in the amr_grid_data's patches.
void amr_grid_data_fill_ghosts(amr_grid_data_t* grid_data);

// Begins an asynchronous ghost-cell-filling operation.
// communicating with other grids as needed. 
void amr_grid_data_start_filling_ghosts(amr_grid_data_t* grid_data);

// Concludes an asynchronous ghost-cell-filling operation initiated by 
// a call to amr_grid_data_start_filling_ghosts().
void amr_grid_data_finish_filling_ghosts(amr_grid_data_t* grid_data);

#endif

