// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_GRID_INTERPOLATOR_H
#define POLYAMRI_AMR_GRID_INTERPOLATOR_H

#include "polyamri/amr_patch.h"

// A grid interpolator interpolates quantities from one or more source patches 
// (to which it is bound for its lifetime) to various destination patches given 
// at the time of interpolation requestions. The patches usually contain cells 
// at different levels of refinement.
typedef struct amr_grid_interpolator_t amr_grid_interpolator_t;

// This function interpolates values from one patch to another.
typedef void (*amr_grid_interpolator_interpolate_func)(void* context, int si1, int si2, int sj1, int sj2, int sk1, int sk2, int di1, int di2, int dj1, int dj2, int dk1, int dk2, amr_patch_t* dest_patch);

// This function destroys an interpolator context.
typedef void (*amr_grid_interpolator_dtor)(void* context);

// This virtual table is used to define the behavior of an integrator.
typedef struct
{
  amr_grid_interpolator_interpolate_func interpolate;
  amr_grid_interpolator_dtor dtor;
} amr_grid_interpolator_vtable;

// Creates a new grid interpolator with the given name, behavior, and 
// source patch.
amr_grid_interpolator_t* amr_grid_interpolator_new(const char* name,
                                                   void* context,
                                                   amr_grid_interpolator_vtable vtable);

// Destroys the given grid interpolator.
void amr_grid_interpolator_free(amr_grid_interpolator_t* interp);

// Allows access to the context pointer within the interpolator.
void* amr_grid_interpolator_context(amr_grid_interpolator_t* interp);

// Interpolates values from the interpolator's source patch to the given 
// destination patch using the given source and destination indices. The 
// indices signify rectangular regions [src_i1, src_i2) x [src_j1, src_j2) x [src_k1, src_k2) 
// and so on.
void amr_grid_interpolator_interpolate(amr_grid_interpolator_t* interp, 
                                       int src_i1, int src_i2, int src_j1, int src_j2, int src_k1, int src_k2,
                                       int dest_i1, int dest_i2, int dest_j1, int dest_j2, int dest_k1, int dest_k2,
                                       amr_patch_t* dest_patch);

#endif

