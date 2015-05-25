// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "linear_amr_grid_interpolator.h"

static void static_linear_interpolate(void* context, int si1, int si2, int sj1, int sj2, int sk1, int sk2, 
                                      int di1, int di2, int dj1, int dj2, int dk1, int dk2, amr_patch_t* dest_patch)
{
  // FIXME
}

amr_grid_interpolator_t* static_linear_amr_grid_interpolator_new(amr_patch_t* source_tile)
{
  amr_grid_interpolator_vtable vtable = {.interpolate = static_linear_interpolate};
  return amr_grid_interpolator_new("Static linear interpolator", source_tile, vtable);
}

typedef struct 
{
  real_t alpha;
  amr_patch_t *patch1, *patch2;
} dyn_interp_t;

static void dynamic_linear_interpolate(void* context, int si1, int si2, int sj1, int sj2, int sk1, int sk2, 
                                       int di1, int di2, int dj1, int dj2, int dk1, int dk2, amr_patch_t* dest_tile)
{
  // FIXME
}

static void dynamic_free(void* context)
{
  dyn_interp_t* dyn = context;
  polymec_free(dyn);
}

amr_grid_interpolator_t* dynamic_linear_amr_grid_interpolator_new(amr_patch_t* patch1, amr_patch_t* patch2)
{
  ASSERT(patch1 != patch2);

  dyn_interp_t* dyn = polymec_malloc(sizeof(dyn_interp_t));
  dyn->alpha = 0.0;
  dyn->patch1 = patch1;
  dyn->patch2 = patch2;
  amr_grid_interpolator_vtable vtable = {.interpolate = dynamic_linear_interpolate, 
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

