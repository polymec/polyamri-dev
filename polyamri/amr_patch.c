// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/polyamri.h"
#include "polyamri/amr_patch.h"

amr_patch_t* amr_patch_new(int ni, int nj, int nk, int nc, int ng)
{
  ASSERT(ni > 0);
  ASSERT(nj > 0);
  ASSERT(nk > 0);
  ASSERT(nc > 0);
  ASSERT(ng >= 0);

  // We allocate one big slab of memory for storage and lean on C99's 
  // VLA semantics.
  size_t storage_size = sizeof(real_t) * (ni+2*ng) * (nj+2*ng) * (nk+2*ng) * nc;
  amr_patch_t* t = polymec_malloc(sizeof(amr_patch_t) + storage_size);
  t->data = (char*)t + sizeof(amr_patch_t);
  t->nc = nc;
  t->i1 = ng;
  t->i2 = ni + ng;
  t->j1 = ng;
  t->j2 = nj + ng;
  t->k1 = ng;
  t->k2 = nk + ng;
  return t;
}

void amr_patch_free(amr_patch_t* t)
{
  polymec_free(t);
}

