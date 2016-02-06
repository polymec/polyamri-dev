// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_FACE_DATA_H
#define POLYAMRI_STR_GRID_FACE_DATA_H

#include "core/point.h"
#include "polyamri/str_grid.h"

// A structured grid face data object is a collection of face-centered patches 
// that are associated with a structured grid in 3D space. Face-centered patches
// do not possess ghost values, and faces on neighboring patches overlap.
typedef struct str_grid_face_data_t str_grid_face_data_t;

// Creates a str_grid_face_data object associated with the given grid, with 
// the given number of components.
str_grid_face_data_t* str_grid_face_data_new(str_grid_t* grid, 
                                             int num_components);

// Creates a str_grid_face_data object associated with the given grid, with 
// patches whose data is aliased to data in the given patch buffer, which is 
// a serialized into a sequential buffer. This object does not own the data 
// in the buffer--it only accesses it. It is up to the caller to ensure that 
// the lifetime of the buffer exceeds that of the resulting str_grid_face_data
// object. NOTE: the buffer can be NULL as long as no patch data is accessed, 
// and can be set using str_grid_face_data_set_buffer below.
str_grid_face_data_t* str_grid_face_data_with_buffer(str_grid_t* grid, 
                                                     int num_components, 
                                                     void* buffer);

// Frees the given str_grid_face_data object.
void str_grid_face_data_free(str_grid_face_data_t* face_data);

// Returns the number of components in the str_grid_face_data object.
int str_grid_face_data_num_components(str_grid_face_data_t* face_data);

// Returns the total number of faces in the str_grid_face_data object. 
int str_grid_face_data_num_faces(str_grid_face_data_t* face_data);

// Returns an internal pointer to the given object's underlying str_grid.
str_grid_t* str_grid_face_data_grid(str_grid_face_data_t* face_data);

// Given a tuple (i, j, k) identifying a patch in the underlying str_grid,
// returns a corresponding patch containing data for x-faces, or NULL if the 
// patch is not present in the str_grid. 
str_grid_patch_t* str_grid_face_data_x_patch(str_grid_face_data_t* face_data, int i, int j, int k);

// Given a tuple (i, j, k) identifying a patch in the underlying str_grid,
// returns a corresponding patch containing data for y-faces, or NULL if the 
// patch is not present in the str_grid. 
str_grid_patch_t* str_grid_face_data_y_patch(str_grid_face_data_t* face_data, int i, int j, int k);

// Given a tuple (i, j, k) identifying a patch in the underlying str_grid,
// returns a corresponding patch containing data for z-faces, or NULL if the 
// patch is not present in the str_grid. 
str_grid_patch_t* str_grid_face_data_z_patch(str_grid_face_data_t* face_data, int i, int j, int k);

// Traverses the grid data, returning true if an x-face patch was found and false if not.
// Set *pos to 0 to reset the traversal. x_patch is set to the x-face patch.
bool str_grid_face_data_next_x_patch(str_grid_face_data_t* face_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** x_patch);

// Traverses the grid data, returning true if a y-face patch was found and false if not.
// Set *pos to 0 to reset the traversal. y_patch is set to the y-face patch.
bool str_grid_face_data_next_y_patch(str_grid_face_data_t* face_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** y_patch);

// Traverses the grid data, returning true if an z-face patch was found and false if not.
// Set *pos to 0 to reset the traversal. z_patch is set to the z-face patch.
bool str_grid_face_data_next_z_patch(str_grid_face_data_t* face_data, int* pos, 
                                     int* i, int* j, int* k, 
                                     str_grid_patch_t** z_patch);

// Returns the pointer to the underlying patch data buffer.
void* str_grid_face_data_buffer(str_grid_face_data_t* face_data);

// Resets the pointer to the underlying patch data buffer, destroying or 
// releasing all existing patch data. If assume_control is true, the 
// str_grid_face_data object will assume control over the buffer and free it 
// upon destruction--otherwise it is assumed to be managed elsewhere.
void str_grid_face_data_set_buffer(str_grid_face_data_t* face_data, 
                                   void* buffer, 
                                   bool assume_control);

#endif

