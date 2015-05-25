// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "polyamri/amr_patch_set.h"

void test_ctor(void** state) 
{
  amr_patch_set_t* patches = amr_patch_set_new();
  assert_int_equal(0, amr_patch_set_size(patches));
  int pos = 0;
  amr_patch_t* patch;
  bbox_t* bbox;
  bool result = amr_patch_set_next(patches, &pos, &patch, &bbox);
  assert_false(result);
  amr_patch_set_free(patches);
}

static void test_uniform_amr_patch_set(void** state, int resolution, int nt, void* data) 
{ 
  // We set up a set of patches that span a region of space from [0, 1] x [0, 1] x [0, 1]. 
  // The region has the given resolution on a side, which we break into 
  // into ntxntxnt patches.
  amr_patch_set_t* patches = amr_patch_set_new();
  for (int i = 0; i < nt; ++i)
  {
    for (int j = 0; j < nt; ++j)
    {
      for (int k = 0; k < nt; ++k)
      {
        // Create the patch (single component, no ghosts).
        amr_patch_t* patch = amr_patch_new(resolution/nt, resolution/nt, resolution/nt, 1, 0);
        // Create the region of space for this patch.
        real_t x1 = 1.0*i/nt, x2 = 1.0*(i+1)/nt,
               y1 = 1.0*j/nt, y2 = 1.0*(j+1)/nt,
               z1 = 1.0*k/nt, z2 = 1.0*(k+1)/nt;
        bbox_t* bbox = bbox_new(x1, x2, y1, y2, z1, z2);
        amr_patch_set_add_with_data(patches, patch, bbox, data, NULL);
      }
    }
  }
  assert_int_equal(nt*nt*nt, amr_patch_set_size(patches));

  // Now traverse the patches.
  int pos = 0, num_patches = 0;
  amr_patch_t* patch;
  bbox_t* bbox;
  void* my_data;
  while (amr_patch_set_next_with_data(patches, &pos, &patch, &bbox, &my_data))
  {
    assert_int_equal(patch->i2 - patch->i1, resolution/nt);
    assert_int_equal(patch->j2 - patch->j1, resolution/nt);
    assert_int_equal(patch->k2 - patch->k1, resolution/nt);
    assert_true(bbox->x2 > bbox->x1);
    assert_true(bbox->y2 > bbox->y1);
    assert_true(bbox->z2 > bbox->z1);
    assert_true(my_data == data);
    ++num_patches;
  }
  assert_int_equal(nt*nt*nt, num_patches);

  // Clean up.
  amr_patch_set_free(patches); 
} 

void test_128x128x128_set(void** state)
{
  test_uniform_amr_patch_set(state, 128, 16, NULL);
}

void test_128x128x128_set_with_data(void** state)
{
  int data = 67;
  test_uniform_amr_patch_set(state, 128, 16, &data);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ctor),
    unit_test(test_128x128x128_set),
    unit_test(test_128x128x128_set_with_data)
  };
  return run_tests(tests);
}
