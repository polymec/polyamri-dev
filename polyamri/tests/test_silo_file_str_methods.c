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
#include "polyamri/silo_file_str_methods.h"
#include "polyamri/grid_to_bbox_coord_mapping.h"

void test_write_str_grid_patch(void** state) 
{ 
  str_grid_patch_t* patch = str_grid_patch_new(10, 10, 10, 1, 0); 
  DECLARE_STR_GRID_PATCH_ARRAY(a, patch);
  for (int i = patch->i1; i < patch->i2; ++i)
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int l = 0; l < 4; ++l)
          a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  const char* field = "solution";
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_str_methods", "test_write_str_grid_patch", 1, 0, 0, 0.0);
  silo_file_write_str_grid_patch(silo, &field, "patch_without_bbox", patch, NULL, NULL);
  bbox_t bbox = {.x1 = -50.0, .x2 = 50.0,
                 .y1 = -25.0, .y2 = 25.0,
                 .z1 = -12.5, .z2 = 12.5};
  coord_mapping_t* mapping = grid_to_bbox_coord_mapping_new(&bbox);
  silo_file_write_str_grid_patch(silo, &field, "patch_with_bbox", patch, mapping, NULL);
  silo_file_close(silo);
  str_grid_patch_free(patch); 
} 

void test_write_str_grid(void** state) 
{ 
  // Make a grid of 4x4x4 patches, each with 10x10x10 cells, 
  // and no periodicity.
  str_grid_t* grid = str_grid_new(4, 4, 4, 10, 10, 10, false, false, false); 
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
        str_grid_insert_patch(grid, i, j, k);
  str_grid_finalize(grid);

  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_str_methods", "test_write_str_grid", 1, 0, 0, 0.0);
  silo_file_write_str_grid(silo, "grid", grid, NULL);
  silo_file_close(silo);

  str_grid_free(grid); 
} 

void test_write_str_grid_cell_data(void** state) 
{ 
  // Make a grid of 4x4x4 patches, each with 10x10x10 cells, 
  // and no periodicity.
  str_grid_t* grid = str_grid_new(4, 4, 4, 10, 10, 10, false, false, false); 
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
        str_grid_insert_patch(grid, i, j, k);
  str_grid_finalize(grid);

  // Make a 4-component solution on this grid (no ghosts).
  str_grid_cell_data_t* solution = str_grid_cell_data_new(grid, 4, 0);
  
  // Fill it with goodness.
  int pos = 0, I, J, K;
  str_grid_patch_t* patch;
  while (str_grid_cell_data_next_patch(solution, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_STR_GRID_PATCH_ARRAY(a, patch);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int j = patch->j1; j < patch->j2; ++j)
        for (int k = patch->k1; k < patch->k2; ++k)
          for (int l = 0; l < 4; ++l)
            a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  }

  // Plot the thing.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_str_methods", "test_write_str_grid_cell_data", 1, 0, 0, 0.0);
  silo_file_write_str_grid(silo, "grid", grid, NULL);
  const char* field_names[4] = {"sol1", "sol2", "sol3", "sol4"};
  silo_file_write_str_grid_cell_data(silo, field_names, "grid", solution, NULL, NULL);
  silo_file_close(silo);

  str_grid_cell_data_free(solution); 
  str_grid_free(grid); 
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_write_str_grid_patch),
    cmocka_unit_test(test_write_str_grid),
    cmocka_unit_test(test_write_str_grid_cell_data)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
