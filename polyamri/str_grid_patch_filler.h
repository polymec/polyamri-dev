// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_PATCH_FILLER_H
#define POLYAMRI_STR_GRID_PATCH_FILLER_H

#include "polyamri/str_grid_cell_data.h"

// A structured grid patch filler fills ghost cells in patches, either by
// applying some boundary condition or by copying values from somewhere else.
// This is a base class that defines an interface for the filling operation.
// Objects of this type are garbage-collected.
typedef struct str_grid_patch_filler_t str_grid_patch_filler_t;

// Every implementation of a structured grid patch filler must fill this
// virtual table.
typedef struct
{
  // Begin the process of filling ghost cells in the patch (which exists at patch 
  // coordinates (i, j, k) in its underlying grid), returning an integer-valued token 
  // with which the process can be completed, or -1 if the process is completed synchronously.
  int (*start_filling_cells)(void* context, 
                             int i, int j, int k, 
                             str_grid_cell_data_t* cell_data);

  // Finish filling ghost cells using the given token.
  void (*finish_filling_cells)(void* context, int token);

  // Destructor.
  void (*dtor)(void* context);
} str_grid_patch_filler_vtable;

// Creates a new structured grid patch filler with the given name, context,
// and virtual table.
str_grid_patch_filler_t* str_grid_patch_filler_new(const char* name,
                                                   void* context,
                                                   str_grid_patch_filler_vtable vtable);

// Start filling the ghost cells in the given patch (situated at patch coordinates (i, j, k)
// in its underlying grid) asynchronously, returning a non-negative token that can be used 
// later to finish the operation, or -1 if the process is completed synchronously.
int str_grid_patch_filler_start(str_grid_patch_filler_t* filler,
                                int i, int j, int k,
                                str_grid_cell_data_t* cell_data);

// Finish filling the ghost cells in the given patch by concluding the fill
// operation that corresponds to the given token.
void str_grid_patch_filler_finish(str_grid_patch_filler_t* filler,
                                  int token);

//------------------------------------------------------------------------
//                  Commonly-used grid patch fillers
//------------------------------------------------------------------------

// Creates a patch filler that copies local values from one patch on a grid to 
// another.
str_grid_patch_filler_t* copy_str_grid_patch_filler_new(str_grid_t* grid,
                                                        str_grid_patch_boundary_t src_boundary,
                                                        str_grid_patch_boundary_t dest_boundary);

#endif

