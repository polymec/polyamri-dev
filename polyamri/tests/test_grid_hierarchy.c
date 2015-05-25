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

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "polyamri/grid_hierarchy.h"

void test_ctor(void** state) 
{
  bbox_t domain = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  grid_hierarchy_t* h = grid_hierarchy_new(&domain, 4, 4, 4, 16, 16, 16, 1, 2, false, false, false);
  assert_int_equal(0, grid_hierarchy_num_levels(h));
  assert_int_equal(1, grid_hierarchy_num_ghosts(h));
  assert_int_equal(2, grid_hierarchy_ref_ratio(h));
  bool periodicity[3];
  grid_hierarchy_get_periodicity(h, periodicity);
  assert_false(periodicity[0]);
  assert_false(periodicity[1]);
  assert_false(periodicity[2]);

  grid_level_t* level = grid_hierarchy_add_level(h);
  assert_int_equal(1, grid_hierarchy_num_levels(h));
  grid_level_get_periodicity(level, periodicity);
  assert_false(periodicity[0]);
  assert_false(periodicity[1]);
  assert_false(periodicity[2]);

  // Test basic coarse->fine traversal.
  int pos = 0;
  grid_level_t* l;
  bool result = grid_hierarchy_next_coarsest(h, &pos, &l);
  assert_true(result);
  assert_true(l == level);
  assert_false(grid_hierarchy_next_coarsest(h, &pos, &l));

  // Test basic fine->coarse traversal.
  pos = 0;
  result = grid_hierarchy_next_finest(h, &pos, &l);
  assert_true(result);
  assert_true(l == level);
  assert_false(grid_hierarchy_next_finest(h, &pos, &l));

  grid_hierarchy_free(h);
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
