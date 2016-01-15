// Copyright (c) 2014-2016, Jeffrey N. Johnson
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

#ifndef POLYAMRI_LINEAR_AMR_GRID_INTERPOLATOR_H
#define POLYAMRI_LINEAR_AMR_GRID_INTERPOLATOR_H

#include "polyamri/amr_grid_interpolator.h"

// Creates a "static" interpolator that linearly interpolates data from the
// given source patch.
amr_grid_interpolator_t* static_linear_amr_grid_interpolator_new();

// Creates a "dynamic" interpolator that linearly interpolates data from a 
// linear combination of source patches (1 - alpha) * patch1 + alpha * patch2.
amr_grid_interpolator_t* dynamic_linear_amr_grid_interpolator_new(real_t alpha);

// Sets the interpolation coefficient alpha for the dynamic linear grid 
// interpolator.
void dynamic_linear_amr_grid_interpolator_set_alpha(amr_grid_interpolator_t* interp, real_t alpha);

#endif

