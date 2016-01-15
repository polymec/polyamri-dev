// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "linear_amr_grid_interpolator.h"

static void static_linear_interpolate(void* context, int si1, int si2, int sj1, int sj2, int sk1, int sk2, 
                                      int di1, int di2, int dj1, int dj2, int dk1, int dk2, amr_patch_t* dest_patch)
{
  // FIXME
}

amr_grid_interpolator_t* static_linear_amr_grid_interpolator_new()
{
  amr_grid_interpolator_vtable vtable = {.interpolate = static_linear_interpolate};
  return amr_grid_interpolator_new("Static linear interpolator", NULL, vtable);
}

typedef struct 
{
  real_t alpha;
} dyn_interp_t;

static void dynamic_linear_interpolate(void* context, int si1, int si2, int sj1, int sj2, int sk1, int sk2, 
                                       int di1, int di2, int dj1, int dj2, int dk1, int dk2, amr_patch_t* dest_tile)
{
  // FIXME
}

static void* dynamic_clone(void* context)
{
  dyn_interp_t* dyn = context;
  dyn_interp_t* new_dyn = polymec_malloc(sizeof(dyn_interp_t));
  memcpy(new_dyn, dyn, sizeof(dyn_interp_t));
  return new_dyn;
}

static void dynamic_free(void* context)
{
  dyn_interp_t* dyn = context;
  polymec_free(dyn);
}

amr_grid_interpolator_t* dynamic_linear_amr_grid_interpolator_new(real_t alpha)
{
  dyn_interp_t* dyn = polymec_malloc(sizeof(dyn_interp_t));
  dyn->alpha = alpha;
  amr_grid_interpolator_vtable vtable = {.interpolate = dynamic_linear_interpolate, 
                                         .clone = dynamic_clone,
                                         .dtor = dynamic_free};
  return amr_grid_interpolator_new("Dynamic linear interpolator", dyn, vtable);
}

void dynamic_linear_amr_grid_interpolator_set_alpha(amr_grid_interpolator_t* interp, real_t alpha)
{
  ASSERT(alpha >= 0.0);
  ASSERT(alpha <= 1.0);

  dyn_interp_t* dyn = amr_grid_interpolator_context(interp);
  ASSERT((dyn->alpha >= 0.0) && (dyn->alpha <= 1.0)); // shouldn't be garbage here.
  dyn->alpha = alpha;
}

