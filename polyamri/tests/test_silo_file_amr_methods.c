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
#include "polyamri/silo_file_amr_methods.h"

void test_write_amr_patch(void** state) 
{ 
  amr_patch_t* patch = amr_patch_new(10, 10, 10, 4, 0); 
  DECLARE_AMR_PATCH_ARRAY(a, patch);
  for (int i = patch->i1; i < patch->i2; ++i)
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int l = 0; l < 4; ++l)
          a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_amr_methods", "test_write_amr_patch", 1, 0, 0, 0.0);
  silo_file_write_amr_patch(silo, "patch_without_bbox", patch, NULL);
  bbox_t bbox = {.x1 = -50.0, .x2 = 50.0,
                 .y1 = -25.0, .y2 = 25.0,
                 .z1 = -12.5, .z2 = 12.5};
  silo_file_write_amr_patch(silo, "patch_with_bbox", patch, &bbox);
//  sp_func_t* mapping = 
//  silo_file_write_mapped_amr_patch(silo, "mapped_patch", patch, mapping);
  silo_file_close(silo);
  amr_patch_free(patch); 
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_write_amr_patch)
  };
  return run_tests(tests);
}
