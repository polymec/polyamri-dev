// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_GRID_H
#define POLYAMRI_AMR_GRID_H

#include "polyamri/amr_grid_interpolator.h"

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
  AMR_GRID_X1_NEIGHBOR = 0,
  AMR_GRID_X2_NEIGHBOR = 1,
  AMR_GRID_Y1_NEIGHBOR = 2,
  AMR_GRID_Y2_NEIGHBOR = 3,
  AMR_GRID_Z1_NEIGHBOR = 4,
  AMR_GRID_Z2_NEIGHBOR = 5
} amr_grid_neighbor_slot_t;

//------------------------------------------------------------------------
//                          Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct AMR grids and associate them 
// with their neighbors and with their coarser/finer counterparts.
// amr_grid_finalize() must be called after an AMR grid has been properly
// constructed.
//------------------------------------------------------------------------

// Creates a new empty grid level defined on the region filling [0,1]x[0,1]x[0,1]
// with nx x ny x nz patches of size px x py x pz. Each patch has the given 
// number of ghost cells. This grid is not associated with any other grids.
amr_grid_t* amr_grid_new(int nx, int ny, int nz, 
                         int px, int py, int pz,
                         int num_ghosts,
                         bool periodic_in_x, 
                         bool periodic_in_y, 
                         bool periodic_in_z);

// Destroys the given grid and all of its patches.
void amr_grid_free(amr_grid_t* grid);

// Sets the given AMR grid as a neighbor of this one. "The neighbor slot" 
// identifies which of the neighbor "slots" will be occupied by the neighbor
// grid. The given interpolator is used to transfer data from the other grid 
// to this one. The AMR grid borrows the reference to the neighbor grid.
void amr_grid_set_neighbor(amr_grid_t* grid, 
                           amr_grid_neighbor_slot_t neighbor_slot,
                           amr_grid_t* neighbor,
                           amr_grid_interpolator_t* interpolator);

// Associates a named property (datum) with this grid. An optional 
// destructor function specifies how the datum should be destroyed when the 
// grid is destroyed.
void amr_grid_set_property(amr_grid_t* grid,
                           const char* name,
                           void* property,
                           void (*property_dtor)(void* property));

// Associates a finer grid with this one, with the given refinement ratio 
// (which must be a power of 2). The given interpolator is used to transfer 
// data from the finer grid to this one.
void amr_grid_associate_finer(amr_grid_t* grid, 
                              amr_grid_t* finer_grid, 
                              int ref_ratio,
                              amr_grid_interpolator_t* interpolator);

// Associates a coarser grid with this one, with the given refinement ratio 
// (which must be a power of 2). The given interpolator is used to transfer
// data from the coarser grid to this one.
void amr_grid_associate_coarser(amr_grid_t* grid, 
                                amr_grid_t* coarser_grid, 
                                int ref_ratio,
                                amr_grid_interpolator_t* interpolator);

// Inserts a new local patch at (i, j, k) in the nx x ny x nz array of patches.
// No patch may exist (locally, remotely, or at another grid) at (i, j, k).
void amr_grid_add_local_patch(amr_grid_t* grid, int i, int j, int k);

// Inserts a new remote patch at (i, j, k) in the nx x ny x nz array of patches.
// No patch may exist (locally, remotely, or at another grid) at (i, j, k).
// The MPI rank of the remote process is given by remote_owner.
void amr_grid_add_remote_patch(amr_grid_t* grid, int i, int j, int k, int remote_owner);

// Finalizes the construction process for an AMR grid. This must be called 
// before any of the grid's usage methods (below) are invoked. Should only 
// be called once.
void amr_grid_finalize(amr_grid_t* grid);

//------------------------------------------------------------------------
//                          Usage methods
//------------------------------------------------------------------------
// The following methods can only be used after an AMR grid has been 
// fully constructed and finalized.
//------------------------------------------------------------------------

// Fetches the number of patches in this grid in the x, y, and z directions, 
// placing them in nx, ny, nz.
void amr_grid_get_extents(amr_grid_t* grid, int* nx, int* ny, int* nz);

// Fetches the number of cells in each patch on this grid in the x, y, and z 
// directions, placing them in pnx, pny, pnz. Additionally, the number of ghost 
// cells is placed in png.
void amr_grid_get_patch_size(amr_grid_t* grid, int* pnx, int* pny, int* pnz, int* png);

// Returns the number of patches that can be stored locally on this grid.
int amr_grid_num_local_patches(amr_grid_t* grid);

// Traverses the locally-present patches in the grid, returning true and the 
// next (i, j, k) triple if the traversal is incomplete, false otherwise. 
// Set *pos to zero to reset the traversal.
bool amr_grid_next_local_patch(amr_grid_t* grid, int* pos, int* i, int* j, int* k);

// Returns true if the grid has a local patch at (i, j, k), false if not.
bool amr_grid_has_patch(amr_grid_t* grid, int i, int j, int k);

// Queries the periodicity of the grid, placing booleans for the 
// x, y, and z periodicity into the given (3-wide) periodicity array.
void amr_grid_get_periodicity(amr_grid_t* grid, bool* periodicity);

// Retrieves the property with the given name from this grid, or returns 
// NULL if no data is associated using that index.
void* amr_grid_property(amr_grid_t* grid, const char* name);

//------------------------------------------------------------------------
//                          Service methods
//------------------------------------------------------------------------
// The following methods provide services to other classes. They can be 
// used directly by folks who Really Understand What They're Doing.
//------------------------------------------------------------------------

typedef struct amr_grid_data_t amr_grid_data_t;

// Fills all ghost cells in the patches within the given grid data, 
// communicating with other grids as needed.
void amr_grid_fill_ghosts(amr_grid_t* grid, amr_grid_data_t* data);

// Begins an asynchronous ghost-cell-filling operation in the patches within 
// the given grid data, communicating with other grids as needed. 
void amr_grid_start_filling_ghosts(amr_grid_t* grid, amr_grid_data_t* data);

// Concludes an asynchronous ghost-cell-filling operation initiated by 
// a call to amr_grid_start_filling_ghosts().
void amr_grid_finish_filling_ghosts(amr_grid_t* grid, amr_grid_data_t* data);

#endif

