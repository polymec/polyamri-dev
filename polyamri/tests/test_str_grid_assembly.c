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
#include "polyamri/str_grid_assembly.h"

void test_ctor(void** state) 
{
  str_grid_assembly_t* assembly = str_grid_assembly_new();
  assert_int_equal(0, str_grid_assembly_num_blocks(assembly));
  assert_true(str_grid_assembly_block(assembly, "NOT HERE!") == NULL);
  int pos = 0;
  char* block_name;
  str_grid_t* block;
  str_grid_t *x1, *x2, *y1, *y2, *z1, *z2;
  assert_false(str_grid_assembly_next_block(assembly, &pos, &block_name, &block,
                                            &x1, &x2, &y1, &y2, &z1, &z2));
  str_grid_assembly_free(assembly);
}

void test_add_block(void** state)
{
  str_grid_assembly_t* assembly = str_grid_assembly_new();
  str_grid_t* block1 = str_grid_new(4, 4, 4, 16, 16, 16, false, false, false);
  str_grid_assembly_add_block(assembly, "block1", block1);
  assert_int_equal(1, str_grid_assembly_num_blocks(assembly));
  assert_true(str_grid_assembly_block(assembly, "NOT HERE!") == NULL);
  assert_true(str_grid_assembly_block(assembly, "block1") != NULL);
  int pos = 0;
  char* block_name;
  str_grid_t *x1, *x2, *y1, *y2, *z1, *z2;
  str_grid_t* block;
  assert_true(str_grid_assembly_next_block(assembly, &pos, &block_name, &block,
                                           &x1, &x2, &y1, &y2, &z1, &z2));
  assert_true(block == block1);
  assert_true(x1 == NULL);
  assert_true(x2 == NULL);
  assert_true(y1 == NULL);
  assert_true(y2 == NULL);
  assert_true(z1 == NULL);
  assert_true(z2 == NULL);
  assert_false(str_grid_assembly_next_block(assembly, &pos, &block_name, &block,
                                            &x1, &x2, &y1, &y2, &z1, &z2));

  str_grid_assembly_free(assembly);
}

void test_connect_blocks(void** state)
{
  str_grid_assembly_t* assembly = str_grid_assembly_new();
  str_grid_t* block1 = str_grid_new(4, 4, 4, 16, 16, 16, false, false, false);
  str_grid_assembly_add_block(assembly, "block1", block1);
  str_grid_t* block2 = str_grid_new(4, 4, 4, 16, 16, 16, false, false, false);
  str_grid_assembly_add_block(assembly, "block2", block2);
  assert_int_equal(2, str_grid_assembly_num_blocks(assembly));
  assert_true(str_grid_assembly_block(assembly, "NOT HERE!") == NULL);
  assert_true(str_grid_assembly_block(assembly, "block1") != NULL);
  assert_true(str_grid_assembly_block(assembly, "block2") != NULL);

  str_grid_assembly_connect(assembly, 
                            "block1", STR_GRID_X2_BOUNDARY, 
                            "block2", STR_GRID_X1_BOUNDARY);

  int pos = 0;
  char* block_name;
  str_grid_t *x1, *x2, *y1, *y2, *z1, *z2;
  str_grid_t* block;
  bool got_block1 = false, got_block2 = false;
  while (str_grid_assembly_next_block(assembly, &pos, &block_name, &block,
                                      &x1, &x2, &y1, &y2, &z1, &z2))
  {
    if (block == block1)
    {
      got_block1 = true;
      assert_true(x1 == NULL);
      assert_true(x2 == block2);
      assert_true(y1 == NULL);
      assert_true(y2 == NULL);
      assert_true(z1 == NULL);
      assert_true(z2 == NULL);
    }
    else
    {
      assert_true(block == block2);
      got_block2 = true;
      assert_true(x1 == block1);
      assert_true(x2 == NULL);
      assert_true(y1 == NULL);
      assert_true(y2 == NULL);
      assert_true(z1 == NULL);
      assert_true(z2 == NULL);
    }
  }

  str_grid_assembly_free(assembly);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctor),
    cmocka_unit_test(test_add_block),
    cmocka_unit_test(test_connect_blocks)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
