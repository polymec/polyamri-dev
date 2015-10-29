// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_GRID_HIERARCHY_H
#define POLYAMRI_AMR_GRID_HIERARCHY_H

#include "polyamri/amr_grid_data.h"

// A grid hierarchy is an object that manages a set of grid levels representing 
// a composite grid for a quantity whose resolution is adaptive in space. The 
// grid levels can be traversed from coarsest to finest or from finest to 
// coarsest.
typedef struct amr_grid_hierarchy_t amr_grid_hierarchy_t;

// Creates a new empty grid hierarchy on [0,1]x[0,1]x[0,1] with the given 
// number of ghost cells in each tile, the (integer power of 2) refinement 
// ratio between levels, and the given periodicity. The coarsest grid level, 
// when added, will have nx x ny x nz tiles of dimension px x py x pz, 
amr_grid_hierarchy_t* amr_grid_hierarchy_new(MPI_Comm comm,
                                             int nx, int ny, int nz, 
                                             int px, int py, int pz,
                                             int ref_ratio,
                                             bool periodic_in_x, bool periodic_in_y, bool periodic_in_z,
                                             amr_grid_interpolator_t* interpolator);

// Destroys the given grid hierarchy and all of its levels.
void amr_grid_hierarchy_free(amr_grid_hierarchy_t* hierarchy);

// Returns the MPI communicator associated with this AMR grid hierarchy.
MPI_Comm amr_grid_hierarchy_comm(amr_grid_hierarchy_t* hierarchy);

// Associates the given MPI communicator with this AMR grid hierarchy.
void amr_grid_hierarchy_set_comm(amr_grid_hierarchy_t* hierarchy, MPI_Comm comm);


// Returns the current number of grid levels.
int amr_grid_hierarchy_num_levels(amr_grid_hierarchy_t* hierarchy);

// Returns the refinement ratio between grid levels.
int amr_grid_hierarchy_ref_ratio(amr_grid_hierarchy_t* hierarchy);

// Queries the periodicity of the grid hierarchy, placing booleans for the 
// x, y, and z periodicity into the given periodicity array.
void amr_grid_hierarchy_get_periodicity(amr_grid_hierarchy_t* hierarchy, bool* periodicity);

// Adds a new (finer) level of refinement to this hierarchy, associating it with the 
// next-coarsest level. Returns the grid just added.
amr_grid_t* amr_grid_hierarchy_add_level(amr_grid_hierarchy_t* hierarchy);

// Traverses the hierarchy of grids from coarsest to finest.
// Set *pos to 0 to reset the iteration.
bool amr_grid_hierarchy_next_coarsest(amr_grid_hierarchy_t* hierarchy, int* pos, amr_grid_t** level);

// Traverses the hierarchy of grids from finest to coarsest.
// Set *pos to 0 to reset the iteration.
bool amr_grid_hierarchy_next_finest(amr_grid_hierarchy_t* hierarchy, int* pos, amr_grid_t** level);

#endif

