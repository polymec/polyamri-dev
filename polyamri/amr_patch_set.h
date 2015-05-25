// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_PATCH_SET_H
#define POLYAMRI_AMR_PATCH_SET_H

#include "core/point.h"
#include "polyamri/amr_patch.h"

// An AMR patch set is a collection of patches that are associated with 
// regions in 3D space.
typedef struct amr_patch_set_t amr_patch_set_t;

// Creates a new empty patch set.
amr_patch_set_t* amr_patch_set_new();

// Adds the given patch to the patch set, associating it with the rectangular 
// region in 3D space defined by the given bounding box. The patch set assumes 
// control of the patch, managing its resources.
void amr_patch_set_add(amr_patch_set_t* patches, amr_patch_t* patch, bbox_t* bbox);

// Adds the given patch to the patch set, associating it with the rectangular 
// region in 3D space defined by the given bounding box AND the given data. 
// The patch set assumes control of the patch, managing its resources, and calls 
// data_dtor to destroy the data when it is freed.
void amr_patch_set_add_with_data(amr_patch_set_t* patches, amr_patch_t* patch, bbox_t* bbox, void* data, void (*data_dtor)(void*));

// Frees the given patch set.
void amr_patch_set_free(amr_patch_set_t* patches);

// Returns the number of patches in the set.
int amr_patch_set_size(amr_patch_set_t* patches);

// Traverses the patch set, returning true if a patch was found and false if not.
// Set *pos to 0 to reset the traversal. Pointers to the patch and its 
// bounding box are returned in the given locations.
bool amr_patch_set_next(amr_patch_set_t* patches, int* pos, amr_patch_t** patch, bbox_t** bbox);

// Traverses the patch set, returning true if a patch was found and false if not.
// Set *pos to 0 to reset the traversal. Pointers to the patch, its bounding 
// box, and its data are returned in the given locations.
bool amr_patch_set_next_with_data(amr_patch_set_t* patches, int* pos, amr_patch_t** patch, bbox_t** bbox, void** data);

#endif

