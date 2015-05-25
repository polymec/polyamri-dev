// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/amr_grid_interpolator.h"

struct amr_grid_interpolator_t 
{
  char* name;
  void* context;
  amr_grid_interpolator_vtable vtable;
};

amr_grid_interpolator_t* amr_grid_interpolator_new(const char* name,
                                                   void* context,
                                                   amr_grid_interpolator_vtable vtable)
{
  ASSERT(vtable.interpolate != NULL);

  amr_grid_interpolator_t* interp = polymec_malloc(sizeof(amr_grid_interpolator_t));
  interp->name = string_dup(name);
  interp->context = context;
  interp->vtable = vtable;
  return interp;
}

void amr_grid_interpolator_free(amr_grid_interpolator_t* interp)
{
  if ((interp->context != NULL) && (interp->vtable.dtor != NULL))
    interp->vtable.dtor(interp->context);
  polymec_free(interp->name);
  polymec_free(interp);
}

void* amr_grid_interpolator_context(amr_grid_interpolator_t* interp)
{
  return interp->context;
}

void amr_grid_interpolator_interpolate(amr_grid_interpolator_t* interp, 
                                       int src_i1, int src_i2, int src_j1, int src_j2, int src_k1, int src_k2,
                                       int dest_i1, int dest_i2, int dest_j1, int dest_j2, int dest_k1, int dest_k2,
                                       amr_patch_t* dest_patch)
{
  interp->vtable.interpolate(interp->context, 
                             src_i1, src_i2, src_j1, src_j2, src_k1, src_k2, 
                             dest_i1, dest_i2, dest_j1, dest_j2, dest_k1, dest_k2, dest_patch);
}

