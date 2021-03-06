// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/polyamri.h"
#include "polyamri/str_grid_patch.h"

str_grid_patch_t* str_grid_patch_new(int ni, int nj, int nk, 
                                     int nc, int ng)
{
  ASSERT(ni > 0);
  ASSERT(nj > 0);
  ASSERT(nk > 0);
  ASSERT(nc > 0);
  ASSERT(ng >= 0);

  // We allocate one big slab of memory for storage and lean on C99's 
  // VLA semantics.
  size_t storage_size = sizeof(real_t) * (ni+2*ng) * (nj+2*ng) * (nk+2*ng) * nc;
  str_grid_patch_t* p = polymec_malloc(sizeof(str_grid_patch_t) + sizeof(real_t) * storage_size);
  p->data = (char*)p + sizeof(str_grid_patch_t);
  memset(p->data, 0, sizeof(real_t) * storage_size);
  p->nc = nc;
  p->ng = ng;
  p->i1 = ng;
  p->i2 = ni + ng;
  p->j1 = ng;
  p->j2 = nj + ng;
  p->k1 = ng;
  p->k2 = nk + ng;
  return p;
}

str_grid_patch_t* str_grid_patch_with_buffer(int ni, int nj, int nk, 
                                             int nc, int ng, void* buffer)
{
  ASSERT(ni > 0);
  ASSERT(nj > 0);
  ASSERT(nk > 0);
  ASSERT(nc > 0);
  ASSERT(ng >= 0);

  str_grid_patch_t* p = polymec_malloc(sizeof(str_grid_patch_t));
  p->data = buffer;
  p->nc = nc;
  p->ng = ng;
  p->i1 = ng;
  p->i2 = ni + ng;
  p->j1 = ng;
  p->j2 = nj + ng;
  p->k1 = ng;
  p->k2 = nk + ng;
  return p;
}

str_grid_patch_t* str_grid_patch_clone(str_grid_patch_t* patch)
{
  int nc = patch->nc;
  int ng = patch->i1;
  int ni = patch->i2 - ng;
  int nj = patch->j2 - ng;
  int nk = patch->k2 - ng;
  return str_grid_patch_new(ni, nj, nk, nc, ng);
}

void str_grid_patch_free(str_grid_patch_t* patch)
{
  polymec_free(patch);
}

