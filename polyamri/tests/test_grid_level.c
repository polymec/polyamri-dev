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
#include "polyamri/grid_level.h"

void test_ctor(void** state) 
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  grid_level_t* level = grid_level_new(&bbox, 4, 4, 4, 16, 16, 16, 1, false, false, false);
  assert_int_equal(0, grid_level_num_tiles(level));
  tile_set_t* tiles = grid_level_tile_set(level, 1);
  assert_int_equal(0, tile_set_size(tiles));
  grid_level_fill_ghosts(level, tiles); // Should do nothing.
  tile_set_free(tiles);
  grid_level_free(level);
}

static void test_fill_ghosts(void** state)
{ 
  // Set up a 4 x 4 x 4 array of tiles in a grid level and fill ghost values.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  grid_level_t* level = grid_level_new(&bbox, 4, 4, 4, 16, 16, 16, 1, true, true, true);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
        grid_level_add_tile(level, i, j, k);
  assert_int_equal(4*4*4, grid_level_num_tiles(level));

  // Make a tile set on this grid level.
  tile_set_t* tiles = grid_level_tile_set(level, 1);
  assert_int_equal(4*4*4, tile_set_size(tiles));

  // Fill each tile with a unique identifier.
  int pos = 0;
  tile_t* tile;
  bbox_t* bb;
  while (tile_set_next(tiles, &pos, &tile, &bb))
  {
    DECLARE_TILE_ARRAY(a, tile);

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
  grid_level_fill_ghosts(level, tiles); 

  // Now check the ghost cell values.
  pos = 0;
  while (tile_set_next(tiles, &pos, &tile, &bb))
  {
    DECLARE_TILE_ARRAY(a, tile);

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
  tile_set_free(tiles);
  grid_level_free(level);
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ctor),
    unit_test(test_fill_ghosts)
  };
  return run_tests(tests);
}
