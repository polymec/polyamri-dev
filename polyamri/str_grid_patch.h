// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_PATCH_H
#define POLYAMRI_STR_GRID_PATCH_H

#include "core/declare_nd_array.h"
#include "polyamri/polyamri.h"

// A structured grid patch (str_grid_patch) is a (3D) ni x nj x nk interval in 
// index space with multi-component data (nc components) at every (i, j, k) 
// triple. It allocates space for ng ghost indices in every dimension. Data are 
// accessed within the data array, e.g. p->data[i][j][k][c].
//
// The index space is bounded by [i1, i2-1] in the "i" direction, [j1, j2-1]
// in the "j" direction, and [k1, k2-1] in the "k" direction. Ghost values pad 
// these indices on either side.
typedef struct 
{
  // Data storage for the patch.
  void* data;
  
  // The number of components.
  int nc;

  // The start and end indices in i, j, k for the patch.
  int i1, i2, j1, j2, k1, k2;

  // The number of ghost values.
  int ng;
} str_grid_patch_t;

// This macro generates a multidimensional array that can access the given
// patch's data using C99 variable-length arrays.
// Note that patch->i1 + tile->i2 == (patch->i2 - tile->i1) + 2 * num_ghosts.
#define DECLARE_STR_GRID_PATCH_ARRAY(array, patch) \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->i1 + patch->i2, patch->j1 + patch->j2, patch->k1 + patch->k2, patch->nc)

// Creates a new structured grid patch of the given size with nc components 
// and ng ghost indices.
str_grid_patch_t* str_grid_patch_new(int ni, int nj, int nk, int nc, int ng);

// Creates a deep copy of the structured grid patch.
str_grid_patch_t* str_grid_patch_clone(str_grid_patch_t* patch);

// Frees the given structured grid patch.
void str_grid_patch_free(str_grid_patch_t* patch);

#endif

