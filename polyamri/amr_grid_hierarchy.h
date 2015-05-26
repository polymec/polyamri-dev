// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_GRID_HIERARCHY_H
#define POLYAMRI_AMR_GRID_HIERARCHY_H

#include "polyamri/amr_grid.h"

// A grid hierarchy is an object that manages a set of grid levels representing 
// a composite grid for a quantity whose resolution is adaptive in space. The 
// grid levels can be traversed from coarsest to finest or from finest to 
// coarsest.
typedef struct amr_grid_hierarchy_t amr_grid_hierarchy_t;

// Creates a new empty grid hierarchy on the domain defined by the given 
// bounding box, with the given number of ghost cells in each tile, the 
// (integer power of 2) refinement ratio between levels, and the given 
// periodicity. The coarsest grid level, when added, will have nx x ny x nz 
// tiles of dimension px x py x pz, 
amr_grid_hierarchy_t* amr_grid_hierarchy_new(bbox_t* domain, 
                                             int nx, int ny, int nz, 
                                             int px, int py, int pz,
                                             int num_ghosts, int ref_ratio,
                                             bool periodic_in_x, bool periodic_in_y, bool periodic_in_z);

// Destroys the given grid hierarchy and all of its levels.
void amr_grid_hierarchy_free(amr_grid_hierarchy_t* hierarchy);

// Returns the bounding box defining the domain of this hierarchy.
bbox_t* amr_grid_hierarchy_domain(amr_grid_hierarchy_t* hierarchy);

// Returns the current number of grid levels.
int amr_grid_hierarchy_num_levels(amr_grid_hierarchy_t* hierarchy);

// Returns the number of ghost cells in the tiles within this grid hierarchy.
int amr_grid_hierarchy_num_ghosts(amr_grid_hierarchy_t* hierarchy);

// Returns the refinement ratio between grid levels.
int amr_grid_hierarchy_ref_ratio(amr_grid_hierarchy_t* hierarchy);

// Queries the periodicity of the grid hierarchy, placing booleans for the 
// x, y, and z periodicity into the given periodicity array.
void amr_grid_hierarchy_get_periodicity(amr_grid_hierarchy_t* hierarchy, bool* periodicity);

// Adds another level of refinement to this hierarchy, associating it with 
// the next-coarsest level. Returns the grid just added.
amr_grid_t* amr_grid_hierarchy_add_level(amr_grid_hierarchy_t* hierarchy);

// Traverses the grid hierarchy from the coarsest level to the finest.
// Set *pos to 0 to reset the iteration.
bool amr_grid_hierarchy_next_coarsest(amr_grid_hierarchy_t* hierarchy, int* pos, amr_grid_t** level);

// Traverses the grid hierarchy from the finest level to the coarsest.
// Set *pos to 0 to reset the iteration.
bool amr_grid_hierarchy_next_finest(amr_grid_hierarchy_t* hierarchy, int* pos, amr_grid_t** level);

#endif

