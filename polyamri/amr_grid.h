// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_GRID_H
#define POLYAMRI_AMR_GRID_H

#include "polyamri/amr_patch_set.h"

// An AMR grid is a single level in an AMR hierarchy. It consists of a set 
// of uniformly-sized patches (with any associated data). The grid manages
// these patches, their connectivity with each other, and their connectivity 
// with patches on other processes and other grids.
//
// An AMR grid may have up to 6 other AMR grids as its neighbors: one each 
// for each of the +/- x, y, and z directions.
typedef struct amr_grid_t amr_grid_t;

// Neighbor slots.
typedef enum
{
  AMR_GRID_X1 = 0,
  AMR_GRID_X2 = 1,
  AMR_GRID_Y1 = 2,
  AMR_GRID_Y2 = 3,
  AMR_GRID_Z1 = 4,
  AMR_GRID_Z2 = 5
} amr_grid_neighbor_slot_t;

// Creates a new empty grid level defined on the region filling the given 
// bounding box, with nx x ny x nz patches of size px x py x pz. Each patch has 
// the given number of ghost cells.
// This grid is not associated with any other grids.
amr_grid_t* amr_grid_new(bbox_t* domain, 
                         int nx, int ny, int nz, 
                         int px, int py, int pz,
                         int num_ghosts,
                         bool periodic_in_x, 
                         bool periodic_in_y, 
                         bool periodic_in_z);

// Destroys the given grid and all of its patches.
void amr_grid_free(amr_grid_t* grid);

// Sets the given AMR grid as a neighbor of this one. "The neighbor slot" 
// identifies which of the neighbor "slots" will be occupied by the neighbor
// grid. This AMR grid borrows the reference to the neighbor grid.
void amr_grid_set_neighbor(amr_grid_t* grid, 
                           amr_grid_neighbor_slot_t neighbor_slot,
                           amr_grid_t* neighbor);

// Associates a piece of data with this grid. The association is made using a 
// data index that uniquely identifies the datum. An optional destructor 
// function specifies how the datum should be destroyed when the grid is 
// destroyed.
void amr_grid_set_data(amr_grid_t* grid,
                       int data_index,
                       void* data,
                       void (*data_dtor)(void* data));

// Retrieves a piece of data from this grid given its index, or returns 
// NULL if no data is associated using that index.
void* amr_grid_data(amr_grid_t* grid,
                    int data_index);

// Associates a finer grid with this one, with the given refinement ratio 
// (which must be a power of 2).
void amr_grid_associate_finer(amr_grid_t* grid, amr_grid_t* finer_grid, int ref_ratio);

// Associates a coarser grid with this one, with the given refinement ratio 
// (which must be a power of 2).
void amr_grid_associate_coarser(amr_grid_t* grid, amr_grid_t* coarser_grid, int ref_ratio);

// Returns the bounding box describing the region represented by this grid.
bbox_t* amr_grid_domain(amr_grid_t* grid);

// Queries the periodicity of the grid, placing booleans for the 
// x, y, and z periodicity into the given periodicity array.
void amr_grid_get_periodicity(amr_grid_t* grid, bool* periodicity);

// Inserts a new local patch at (i, j, k) in the nx x ny x nz array of patches.
// No patch may exist (locally, remotely, or at another grid) at (i, j, k).
void amr_grid_add_local_patch(amr_grid_t* grid, int i, int j, int k);

// Inserts a new remote patch at (i, j, k) in the nx x ny x nz array of patches.
// No patch may exist (locally, remotely, or at another grid) at (i, j, k).
// The MPI rank of the remote process is given by remote_owner.
void amr_grid_add_remote_patch(amr_grid_t* grid, int i, int j, int k, int remote_owner);

// Returns the number of patches in this grid.
int amr_grid_num_patches(amr_grid_t* grid);

// Creates a new set of AMR patches associated with the given grid. This 
// patch set must be deallocated with amr_patch_set_free.
amr_patch_set_t* amr_grid_create_patches(amr_grid_t* grid, int num_components);

// Fills all ghost cells in the patches within the given patch set, 
// communicating with other grids as needed.
void amr_grid_fill_ghosts(amr_grid_t* grid, amr_patch_set_t* patches);

// Begins an asynchronous ghost-cell-filling operation in the patches within 
// the given patch set, communicating with other grids as needed. 
void amr_grid_start_filling_ghosts(amr_grid_t* grid, amr_patch_set_t* patches);

// Concludes an asynchronous ghost-cell-filling operation initiated by 
// a call to amr_grid_start_filling_ghosts().
void amr_grid_finish_filling_ghosts(amr_grid_t* grid, amr_patch_set_t* patches);

#endif

