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
#include "polyamri/amr_grid_hierarchy.h"
#include "polyamri/linear_amr_grid_interpolator.h"

void test_ctor(void** state) 
{
  bbox_t domain = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  amr_grid_interpolator_t* I = static_linear_amr_grid_interpolator_new();
  amr_grid_hierarchy_t* h = amr_grid_hierarchy_new(&domain, 4, 4, 4, 16, 16, 16, 1, 2, false, false, false, I);
  assert_int_equal(0, amr_grid_hierarchy_num_levels(h));
  assert_int_equal(1, amr_grid_hierarchy_num_ghosts(h));
  assert_int_equal(2, amr_grid_hierarchy_ref_ratio(h));
  bool periodicity[3];
  amr_grid_hierarchy_get_periodicity(h, periodicity);
  assert_false(periodicity[0]);
  assert_false(periodicity[1]);
  assert_false(periodicity[2]);

  amr_grid_t* level = amr_grid_hierarchy_add_level(h);
  assert_int_equal(1, amr_grid_hierarchy_num_levels(h));
  amr_grid_get_periodicity(level, periodicity);
  assert_false(periodicity[0]);
  assert_false(periodicity[1]);
  assert_false(periodicity[2]);

  // Test basic coarse->fine traversal.
  int pos = 0;
  amr_grid_t* l;
  bool result = amr_grid_hierarchy_next_coarsest(h, &pos, &l);
  assert_true(result);
  assert_true(l == level);
  assert_false(amr_grid_hierarchy_next_coarsest(h, &pos, &l));

  // Test basic fine->coarse traversal.
  pos = 0;
  result = amr_grid_hierarchy_next_finest(h, &pos, &l);
  assert_true(result);
  assert_true(l == level);
  assert_false(amr_grid_hierarchy_next_finest(h, &pos, &l));

  amr_grid_hierarchy_free(h);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ctor)
  };
  return run_tests(tests);
}
