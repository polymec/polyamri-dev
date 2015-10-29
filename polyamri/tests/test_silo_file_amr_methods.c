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
#include "polyamri/grid_to_bbox_coord_mapping.h"

void test_write_amr_patch(void** state) 
{ 
  amr_patch_t* patch = amr_patch_new(10, 10, 10, 1, 0); 
  DECLARE_AMR_PATCH_ARRAY(a, patch);
  for (int i = patch->i1; i < patch->i2; ++i)
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int l = 0; l < 4; ++l)
          a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  const char* field = "solution";
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_amr_methods", "test_write_amr_patch", 1, 0, 0, 0.0);
  silo_file_write_amr_patch(silo, &field, "patch_without_bbox", patch, NULL, NULL);
  bbox_t bbox = {.x1 = -50.0, .x2 = 50.0,
                 .y1 = -25.0, .y2 = 25.0,
                 .z1 = -12.5, .z2 = 12.5};
  coord_mapping_t* mapping = grid_to_bbox_coord_mapping_new(&bbox);
  silo_file_write_amr_patch(silo, &field, "patch_with_bbox", patch, mapping, NULL);
  silo_file_close(silo);
  amr_patch_free(patch); 
} 

void test_write_amr_grid(void** state) 
{ 
  // Make a grid of 4x4x4 patches, each with 10x10x10 cells, 
  // and no periodicity.
  amr_grid_t* grid = amr_grid_new(MPI_COMM_SELF, 4, 4, 4, 10, 10, 10, false, false, false); 
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
        amr_grid_add_local_patch(grid, i, j, k);
  amr_grid_finalize(grid);

  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_amr_methods", "test_write_amr_grid", 1, 0, 0, 0.0);
  silo_file_write_amr_grid(silo, "grid", grid, NULL);
  silo_file_close(silo);

  amr_grid_free(grid); 
} 

void test_write_amr_grid_data(void** state) 
{ 
  // Make a grid of 4x4x4 patches, each with 10x10x10 cells, 
  // and no periodicity.
  amr_grid_t* grid = amr_grid_new(MPI_COMM_SELF, 4, 4, 4, 10, 10, 10, false, false, false); 
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
        amr_grid_add_local_patch(grid, i, j, k);
  amr_grid_finalize(grid);

  // Make a 4-component solution on this grid (no ghosts).
  amr_grid_data_t* solution = amr_grid_data_new(grid, AMR_GRID_CELL, 4, 0);
  
  // Fill it with goodness.
  int pos = 0, I, J, K;
  amr_patch_t* patch;
  while (amr_grid_data_next_local_patch(solution, &pos, &I, &J, &K, &patch))
  {
    DECLARE_AMR_PATCH_ARRAY(a, patch);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int j = patch->j1; j < patch->j2; ++j)
        for (int k = patch->k1; k < patch->k2; ++k)
          for (int l = 0; l < 4; ++l)
            a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  }

  // Plot the thing.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_amr_methods", "test_write_amr_grid_data", 1, 0, 0, 0.0);
  silo_file_write_amr_grid(silo, "grid", grid, NULL);
  const char* field_names[4] = {"sol1", "sol2", "sol3", "sol4"};
  silo_file_write_amr_grid_data(silo, field_names, "grid", solution, NULL, NULL);
  silo_file_close(silo);

  amr_grid_data_free(solution); 
  amr_grid_free(grid); 
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_write_amr_patch),
    unit_test(test_write_amr_grid),
    unit_test(test_write_amr_grid_data)
  };
  return run_tests(tests);
}
