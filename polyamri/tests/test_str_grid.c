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
#include "polyamri/str_grid.h"

void test_ctor(void** state) 
{
  str_grid_t* grid = str_grid_new(4, 4, 4, 16, 16, 16, false, false, false);
  assert_int_equal(0, str_grid_num_patches(grid));
  str_grid_free(grid);
}

static void test_fill_ghosts(void** state)
{ 
  // Set up a 4 x 4 x 4 array of patches in a grid and fill ghost values.
  str_grid_t* grid = str_grid_new(4, 4, 4, 16, 16, 16, true, true, true);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
        str_grid_insert_patch(grid, i, j, k);
  assert_int_equal(4*4*4, str_grid_num_patches(grid));

#if 0
  // Make a patch set on this grid.
  amr_patch_set_t* patches = str_grid_create_patches(grid, 1);
  assert_int_equal(4*4*4, amr_patch_set_size(patches));

  // Fill each patch with a unique identifier.
  int pos = 0;
  amr_patch_t* patch;
  bbox_t* bb;
  while (amr_patch_set_next(patches, &pos, &patch, &bb))
  {
    DECLARE_AMR_PATCH_ARRAY(a, patch);

    for (int i = 1; i < 17; ++i)
      for (int j = 1; j < 17; ++j)
        for (int k = 1; k < 17; ++k)
          a[i][j][k][0] = 16.0*16.0*i + 16.0*j + 1.0*k + 1.0;

    // Ghost cells should be zero, by the way.
    for (int i = 0; i < 16; ++i)
    {
      for (int j = 0; j < 16; ++j)
      {
        assert_true(a[0][i][j][0] == 0.0);
        assert_true(a[17][i][j][0] == 0.0);
        assert_true(a[i][0][j][0] == 0.0);
        assert_true(a[i][17][j][0] == 0.0);
        assert_true(a[i][j][0][0] == 0.0);
        assert_true(a[i][j][17][0] == 0.0);
      }
    }
  }

  // Fill ghosts.
  str_grid_fill_ghosts(grid, patches); 

  // Now check the ghost cell values.
  pos = 0;
  while (amr_patch_set_next(patches, &pos, &patch, &bb))
  {
    DECLARE_AMR_PATCH_ARRAY(a, patch);

    for (int i = 1; i < 17; ++i)
    {
      for (int j = 1; j < 17; ++j)
      {
        assert_true(a[0][i][j][0] != 0.0);
        assert_true(a[0][i][j][0] != 16.0*16.0*i + 16.0*j + 1.0);
        assert_true(a[17][i][j][0] != 0.0);
        assert_true(a[17][i][j][0] != 16.0*16.0*i + 16.0*j + 16.0 + 1.0);
        assert_true(a[i][0][j][0] != 0.0);
        assert_true(a[i][0][j][0] != 16.0*16.0*i + 16.0*j + 1.0);
        assert_true(a[i][17][j][0] != 0.0);
        assert_true(a[i][17][j][0] != 16.0*16.0*i + 16.0*j + 16.0 + 1.0);
        assert_true(a[i][j][0][0] != 0.0);
        assert_true(a[i][j][0][0] != 16.0*16.0*i + 16.0*j + 1.0);
        assert_true(a[i][j][17][0] != 0.0);
        assert_true(a[i][j][17][0] != 16.0*16.0*i + 16.0*j + 16.0 + 1.0);
      }
    }

    // Interior values should be intact, by the way.
    for (int i = 1; i < 17; ++i)
      for (int j = 1; j < 17; ++j)
        for (int k = 1; k < 17; ++k)
          assert_true(a[i][j][k][0] == 16.0*16.0*i + 16.0*j + 1.0*k + 1.0);
  }

  // Clean up.
  amr_patch_set_free(patches);
#endif
  str_grid_free(grid);
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctor),
    cmocka_unit_test(test_fill_ghosts)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
