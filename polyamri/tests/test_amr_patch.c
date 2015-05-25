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
#include "polyamri/amr_patch.h"

void test_amr_patch_without_ghosts(void** state) 
{ 
  amr_patch_t* patch = amr_patch_new(10, 10, 10, 4, 0); 
  assert_int_equal(10, patch->i2 - patch->i1);
  assert_int_equal(10, patch->j2 - patch->j1);
  assert_int_equal(10, patch->k2 - patch->k1);
  assert_int_equal(4, patch->nc);

  DECLARE_AMR_PATCH_ARRAY(a, patch);
  for (int i = patch->i1; i < patch->i2; ++i)
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int l = 0; l < 4; ++l)
          a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  for (int i = 0; i < 10*10*10*4; ++i)
    assert_true(a[0][0][0][i] == (real_t)i);
  amr_patch_free(patch); 
} 

void test_amr_patch_with_ghosts(void** state) 
{ 
  amr_patch_t* patch = amr_patch_new(10, 10, 10, 4, 1); 
  assert_int_equal(10, patch->i2 - patch->i1);
  assert_int_equal(10, patch->j2 - patch->j1);
  assert_int_equal(10, patch->k2 - patch->k1);
  assert_int_equal(4, patch->nc);
  assert_int_equal(1, patch->i1);
  assert_int_equal(1, patch->j1);
  assert_int_equal(1, patch->k1);

  DECLARE_AMR_PATCH_ARRAY(a, patch);
  for (int i = patch->i1-1; i < patch->i2+1; ++i)
    for (int j = patch->j1-1; j < patch->j2+1; ++j)
      for (int k = patch->k1-1; k < patch->k2+1; ++k)
        for (int l = 0; l < 4; ++l)
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
  amr_patch_free(patch); 
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_amr_patch_without_ghosts),
    unit_test(test_amr_patch_with_ghosts)
  };
  return run_tests(tests);
}
