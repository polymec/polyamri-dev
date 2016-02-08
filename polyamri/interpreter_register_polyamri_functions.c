// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/str_grid.h"
#include "model/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static const char* structured_grid_usage = 
  "grid = structured_grid{num_cells = {nx, ny, nz},\n"
  "                       patch_size = {px, py, pz}}\n"
  "  Creates a structured grid spanning the given domain with the given\n"
  "  numbers of cells in the x, y, z directions, comprising patches of the\n"
  "  given dimensions (in cells).";

static int structured_grid(lua_State* lua)
{
  // Check the number of arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 1)
    return luaL_error(lua, structured_grid_usage);

  // Parse the numbers of cells.
  int num_cells[3];
  lua_pushstring(lua, "num_cells");
  lua_gettable(lua, 1); // Reads name from top, replaces with bounds[name].
  if (!lua_issequence(lua, -1))
    return luaL_error(lua, structured_grid_usage);
  int len;
  real_t* seq = lua_tosequence(lua, -1, &len);
  if (len != 3)
    return luaL_error(lua, structured_grid_usage);
  for (int i = 0; i < 3; ++i)
    num_cells[i] = (int)seq[i];
  lua_pop(lua, 1);

  // Parse the patch dimensions.
  int patch_size[3];
  lua_pushstring(lua, "patch_size");
  lua_gettable(lua, 1); // Reads name from top, replaces with bounds[name].
  if (!lua_issequence(lua, -1))
    return luaL_error(lua, structured_grid_usage);
  seq = lua_tosequence(lua, -1, &len);
  if (len != 3)
    return luaL_error(lua, structured_grid_usage);
  for (int i = 0; i < 3; ++i)
    patch_size[i] = (int)seq[i];
  lua_pop(lua, 1);

  // Make sure the patch dimensions square with the grid dimensions.
  const char* axes[3] = {"x", "y", "z"};
  for (int i = 0; i < 3; ++i)
  {
    if ((num_cells[i] % patch_size[i]) != 0)
    {
      luaL_error(lua, "structured_grid: %s patch size (%d) does not evenly divide %s grid size (%d).", 
        axes[i], patch_size[i], axes[i], num_cells[i]);
    }           
  }
  
  // Create the structured grid.
  int nx = num_cells[0]/patch_size[0];
  int ny = num_cells[1]/patch_size[1];
  int nz = num_cells[2]/patch_size[2];
  str_grid_t* grid = str_grid_new(nx, ny, nz, patch_size[0], patch_size[1], patch_size[2], false, false, false);

  // Fill it with patches.
  // Ghost cells are filled with local copies.
  for (int ip = 0; ip < nx; ++ip)
    for (int jp = 0; jp < ny; ++jp)
      for (int kp = 0; kp < nz; ++kp)
        str_grid_insert_patch(grid, ip, jp, kp);

  // Finalize the grid.
  str_grid_finalize(grid);

  lua_pushuserdefined(lua, grid, DTOR(str_grid_free));
  return 1;
}

void interpreter_register_polyamri_functions(interpreter_t* interp)
{
  interpreter_register_function(interp, "structured_grid", structured_grid, 
    docstring_from_string(structured_grid_usage));
}

str_grid_t* interpreter_get_str_grid(interpreter_t* interp, const char* name)
{
  // For now, we simply store structured grids as user-defined objects.
  // Not "watertight," but we'll cross the validation bridge when we come to it.
  return interpreter_get_user_defined(interp, name);
}

