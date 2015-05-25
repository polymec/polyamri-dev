// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_DATA_HIERARCHY_H
#define POLYAMRI_AMR_DATA_HIERARCHY_H

#include "polyamri/amr_grid_hierarchy.h"

// A data hierarchy is an object that manages a set of patch sets representing 
// a quantity on a composite grid. The patch sets for the various refinement 
// levels of the data can be traversed from coarsest to finest or from finest 
// to coarsest.
typedef struct amr_data_hierarchy_t amr_data_hierarchy_t;

// Creates a new data hierarchy object defined on the given grid hierarchy.
// The lifetime of the data hierarchy (and the validity of its data) is tied 
// to that of the grid hierarchy, though the two objects are managed separately.
amr_data_hierarchy_t* amr_data_hierarchy_new(amr_grid_hierarchy_t* grids,
                                             int num_components);

// Destroys the given data hierarchy and all of its levels.
void amr_data_hierarchy_free(amr_data_hierarchy_t* data);

// Returns the grid hierarchy associated with this data hierarchy.
amr_grid_hierarchy_t* amr_data_hierarchy_grids(amr_data_hierarchy_t* data);

// Returns the number of components for the data within this hierarchy.
int amr_data_hierarchy_num_components(amr_data_hierarchy_t* data);

// Returns the current number of grid levels.
int amr_data_hierarchy_num_levels(amr_data_hierarchy_t* data);

// Traverses the data hierarchy from the coarsest level to the finest.
// Allows access to the patch set and the grid level on each level.
// Set *pos to 0 to reset the iteration.
bool amr_data_hierarchy_next_coarsest(amr_data_hierarchy_t* data, int* pos, amr_patch_set_t** patch_set, amr_grid_level_t** level);

// Traverses the grid hierarchy from the finest level to the coarsest.
// Allows access to the patch set and the grid level on each level.
// Set *pos to 0 to reset the iteration.
bool amr_data_hierarchy_next_finest(amr_data_hierarchy_t* data, int* pos, amr_patch_set_t** patch_set, amr_grid_level_t** level);

#endif

