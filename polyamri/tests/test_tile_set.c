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
#include "polyamri/tile_set.h"

void test_ctor(void** state) 
{
  tile_set_t* tiles = tile_set_new();
  assert_int_equal(0, tile_set_size(tiles));
  int pos = 0;
  tile_t* tile;
  bbox_t* bbox;
  bool result = tile_set_next(tiles, &pos, &tile, &bbox);
  assert_false(result);
  tile_set_free(tiles);
}

static void test_uniform_tile_set(void** state, int resolution, int nt, void* data) 
{ 
  // We set up a set of tiles that span a region of space from [0, 1] x [0, 1] x [0, 1]. 
  // The region has the given resolution on a side, which we break into 
  // into ntxntxnt tiles.
  tile_set_t* tiles = tile_set_new();
  for (int i = 0; i < nt; ++i)
  {
    for (int j = 0; j < nt; ++j)
    {
      for (int k = 0; k < nt; ++k)
      {
        // Create the tile (single component, no ghosts).
        tile_t* tile = tile_new(resolution/nt, resolution/nt, resolution/nt, 1, 0);
        // Create the region of space for this tile.
        real_t x1 = 1.0*i/nt, x2 = 1.0*(i+1)/nt,
               y1 = 1.0*j/nt, y2 = 1.0*(j+1)/nt,
               z1 = 1.0*k/nt, z2 = 1.0*(k+1)/nt;
        bbox_t* bbox = bbox_new(x1, x2, y1, y2, z1, z2);
        tile_set_add_with_data(tiles, tile, bbox, data, NULL);
      }
    }
  }
  assert_int_equal(nt*nt*nt, tile_set_size(tiles));

  // Now traverse the tiles.
  int pos = 0, num_tiles = 0;
  tile_t* tile;
  bbox_t* bbox;
  void* my_data;
  while (tile_set_next_with_data(tiles, &pos, &tile, &bbox, &my_data))
  {
    assert_int_equal(tile->i2 - tile->i1, resolution/nt);
    assert_int_equal(tile->j2 - tile->j1, resolution/nt);
    assert_int_equal(tile->k2 - tile->k1, resolution/nt);
    assert_true(bbox->x2 > bbox->x1);
    assert_true(bbox->y2 > bbox->y1);
    assert_true(bbox->z2 > bbox->z1);
    assert_true(my_data == data);
    ++num_tiles;
  }
  assert_int_equal(nt*nt*nt, num_tiles);

  // Clean up.
  tile_set_free(tiles); 
} 

void test_128x128x128_set(void** state)
{
  test_uniform_tile_set(state, 128, 16, NULL);
}

void test_128x128x128_set_with_data(void** state)
{
  int data = 67;
  test_uniform_tile_set(state, 128, 16, &data);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ctor),
    unit_test(test_128x128x128_set),
    unit_test(test_128x128x128_set_with_data)
  };
  return run_tests(tests);
}
