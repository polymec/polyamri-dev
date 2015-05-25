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
#include "polyamri/tile.h"

void test_tile_without_ghosts(void** state) 
{ 
  tile_t* tile = tile_new(10, 10, 10, 4, 0); 
  assert_int_equal(10, tile->i2 - tile->i1);
  assert_int_equal(10, tile->j2 - tile->j1);
  assert_int_equal(10, tile->k2 - tile->k1);
  assert_int_equal(4, tile->nc);

  DECLARE_TILE_ARRAY(a, tile);
  for (int i = tile->i1; i < tile->i2; ++i)
    for (int j = tile->j1; j < tile->j2; ++j)
      for (int k = tile->k1; k < tile->k2; ++k)
        for (int l = 0; l < 4; ++l)
          a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  for (int i = 0; i < 10*10*10*4; ++i)
    assert_true(a[0][0][0][i] == (real_t)i);
  tile_free(tile); 
} 

void test_tile_with_ghosts(void** state) 
{ 
  tile_t* tile = tile_new(10, 10, 10, 4, 1); 
  assert_int_equal(10, tile->i2 - tile->i1);
  assert_int_equal(10, tile->j2 - tile->j1);
  assert_int_equal(10, tile->k2 - tile->k1);
  assert_int_equal(4, tile->nc);
  assert_int_equal(1, tile->i1);
  assert_int_equal(1, tile->j1);
  assert_int_equal(1, tile->k1);

  DECLARE_TILE_ARRAY(a, tile);
  for (int i = tile->i1-1; i < tile->i2+1; ++i)
    for (int j = tile->j1-1; j < tile->j2+1; ++j)
      for (int k = tile->k1-1; k < tile->k2+1; ++k)
        for (int l = 0; l < 4; ++l)
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
  tile_free(tile); 
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_tile_without_ghosts),
    unit_test(test_tile_with_ghosts)
  };
  return run_tests(tests);
}
