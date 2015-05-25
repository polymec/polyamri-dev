// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_GRID_LEVEL_H
#define POLYAMRI_AMR_GRID_LEVEL_H

#include "polyamri/amr_patch_set.h"

// A grid level is a single level in an AMR hierarchy. It consists of a set 
// of uniformly-sized patches (with any associated data). The grid level manages
// these patches, their connectivity with each other, and their connectivity 
// with patches on other processes and other grid levels.
typedef struct amr_grid_level_t amr_grid_level_t;

// Creates a new empty grid level defined on the region filling the given 
// bounding box, with nx x ny x nz patches of size tx x ty x tz. Each patch has 
// the given number of ghost cells.
// This grid level is not associated with any other grid levels.
amr_grid_level_t* amr_grid_level_new(bbox_t* domain, 
                                     int nx, int ny, int nz, 
                                     int tx, int ty, int tz,
                                     int num_ghosts,
                                     bool periodic_in_x, bool periodic_in_y, bool periodic_in_z);

// Destroys the given grid level and all of its patches.
void amr_grid_level_free(amr_grid_level_t* level);

// Associates a finer level with this one, with the given refinement ratio 
// (which must be a power of 2).
void amr_grid_level_associate_finer_level(amr_grid_level_t* level, amr_grid_level_t* finer_level, int ref_ratio);

// Associates a coarser level with this one, with the given refinement ratio 
// (which must be a power of 2).
void amr_grid_level_associate_coarser_level(amr_grid_level_t* level, amr_grid_level_t* coarser_level, int ref_ratio);

// Returns the bounding box describing the region represented by this grid level.
bbox_t* amr_grid_level_domain(amr_grid_level_t* level);

// Queries the periodicity of the grid level, placing booleans for the 
// x, y, and z periodicity into the given periodicity array.
void amr_grid_level_get_periodicity(amr_grid_level_t* level, bool* periodicity);

// Inserts a new patch location at (i, j, k) in the nx x ny x nz array of patches.
// No patch may exist (locally, remotely, or at another level) at (i, j, k).
void amr_grid_level_add_patch(amr_grid_level_t* level, int i, int j, int k);

// Returns the number of patches in this grid level.
int amr_grid_level_num_patches(amr_grid_level_t* level);

// Creates a new patch set associated with the given grid level.
amr_patch_set_t* amr_grid_level_patch_set(amr_grid_level_t* level, int num_components);

// Fills all ghost cells in the patches within the given patch set, 
// communicating with other grid levels as needed.
void amr_grid_level_fill_ghosts(amr_grid_level_t* level, amr_patch_set_t* patches);

// Begins an asynchronous ghost-cell-filling operation in the patches within 
// the given patch set, communicating with other grid levels as needed. 
void amr_grid_level_start_filling_ghosts(amr_grid_level_t* level, amr_patch_set_t* patches);

// Concludes an asynchronous ghost-cell-filling operation initiated by 
// a call to amr_grid_level_start_filling_ghosts().
void amr_grid_level_finish_filling_ghosts(amr_grid_level_t* level, amr_patch_set_t* patches);

#endif

