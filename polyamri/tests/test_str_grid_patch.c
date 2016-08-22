// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "polyamri/str_grid_patch.h"

static void test_str_grid_patch_without_ghosts(void** state) 
{ 
  str_grid_patch_t* patch = str_grid_patch_new(10, 10, 10, 4, 0); 
  assert_int_equal(10, patch->i2 - patch->i1);
  assert_int_equal(10, patch->j2 - patch->j1);
  assert_int_equal(10, patch->k2 - patch->k1);
  assert_int_equal(4, patch->nc);

  DECLARE_STR_GRID_PATCH_ARRAY(a, patch);
  for (int i = patch->i1; i < patch->i2; ++i)
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int l = 0; l < 4; ++l)
          a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  for (int i = 0; i < 10*10*10*4; ++i)
    assert_true(reals_equal(a[0][0][0][i], (real_t)i));
  str_grid_patch_free(patch); 
} 

static void test_str_grid_patch_with_ghosts(void** state) 
{ 
  str_grid_patch_t* patch = str_grid_patch_new(10, 10, 10, 4, 1); 
  assert_int_equal(10, patch->i2 - patch->i1);
  assert_int_equal(10, patch->j2 - patch->j1);
  assert_int_equal(10, patch->k2 - patch->k1);
  assert_int_equal(4, patch->nc);
  assert_int_equal(1, patch->i1);
  assert_int_equal(1, patch->j1);
  assert_int_equal(1, patch->k1);

  DECLARE_STR_GRID_PATCH_ARRAY(a, patch);
  for (int i = patch->i1-1; i < patch->i2+1; ++i)
    for (int j = patch->j1-1; j < patch->j2+1; ++j)
      for (int k = patch->k1-1; k < patch->k2+1; ++k)
        for (int l = 0; l < 4; ++l)
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
  str_grid_patch_free(patch); 
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_str_grid_patch_without_ghosts),
    cmocka_unit_test(test_str_grid_patch_with_ghosts)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
