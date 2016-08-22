// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/interpreter.h"
#include "polyamri/str_grid_factory.h"
#include "polyamri/str_grid_patch_filler.h"
#include "polyamri/grid_to_bbox_coord_mapping.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// Type codes for structured grids and their friends.
static int structured_grid_type_code = -1;
static int patch_filler_type_code = -1;

// This helper returns the grid spacing in the normal direction of the 
// given boundary.
static real_t grid_spacing(str_grid_t* grid, str_grid_boundary_t boundary)
{
  int npx, npy, npz, nx, ny, nz;
  str_grid_get_extents(grid, &npx, &npy, &npz);
  str_grid_get_patch_size(grid, &nx, &ny, &nz);
  real_t h;
  switch(boundary)
  {
    case STR_GRID_X1_BOUNDARY:
    case STR_GRID_X2_BOUNDARY: h = 1.0 / (npx*nx); break;
    case STR_GRID_Y1_BOUNDARY:
    case STR_GRID_Y2_BOUNDARY: h = 1.0 / (npy*ny); break;
    case STR_GRID_Z1_BOUNDARY:
    case STR_GRID_Z2_BOUNDARY: h = 1.0 / (npz*nz);
  }
  return h;
}

static const char* block_usage = 
  "grid = structured_grids.block{num_cells = {nx, ny, nz},\n"
  "                              patch_size = {px, py, pz},\n"
  "                              domain = D,\n"
  "                              regions = {name1 = indicator1,\n"
  "                                         name2 = indicator2,\n"
  "                                         ...\n"
  "                                         nameN = indicatorN}}\n"
  "  Creates a structured grid spanning the given domain with the given\n"
  "  numbers of cells in the x, y, z directions, comprising patches of the\n"
  "  given dimensions (in cells). Optionally, mapping may be a bounding box\n"
  "  or a coordinate mapping object that maps the logical regin [0,1]x[0,1]x[0,1]\n"
  "  to a physical region. Additionally, regions may be specified by\n"
  "  a table mapping region names to indicator functions. An indicator\n"
  "  function I(x, y, z, t) for a region R maps a cell C to the value 1 if\n"
  "  the centroid of C falls within R and 0 otherwise.\n";

static int block(lua_State* lua)
{
  // Check the number of arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_istable(lua, 1))
    return luaL_error(lua, block_usage);

  // Parse the numbers of cells.
  int num_cells[3];
  lua_pushstring(lua, "num_cells"); // pushes key onto stack
  lua_gettable(lua, 1); // replaces key with value
  if (!lua_issequence(lua, -1))
    return luaL_error(lua, block_usage);
  int len;
  real_t* seq = lua_tosequence(lua, -1, &len);
  if (len != 3)
    return luaL_error(lua, block_usage);
  for (int i = 0; i < 3; ++i)
    num_cells[i] = (int)seq[i];
  lua_pop(lua, 1);

  // Parse the patch dimensions.
  int patch_size[3];
  lua_pushstring(lua, "patch_size"); // pushes key onto stack
  lua_gettable(lua, 1); // replaces key with value
  if (!lua_issequence(lua, -1)) 
    return luaL_error(lua, block_usage);
  seq = lua_tosequence(lua, -1, &len);
  if (len != 3)
    return luaL_error(lua, block_usage);
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

  // Parse domain if such an argument is found.
  coord_mapping_t* mapping = NULL;
  lua_pushstring(lua, "mapping"); // pushes key onto stack
  lua_gettable(lua, 1); // replaces key with value
  if (!lua_isnil(lua, -1))
  {
    mapping = lua_tocoordmapping(lua, -1);
    if (mapping == NULL)
    {
      bbox_t* bbox = lua_toboundingbox(lua, -1);
      if (bbox == NULL)
        luaL_error(lua, "structured_grid: mapping must be a coordinate mapping or a bounding box.");
      else
        mapping = grid_to_bbox_coord_mapping_new(bbox);
    }
  }

  // Parse regions if such an argument is found.
  string_ptr_unordered_map_t* regions = NULL;
  lua_pushstring(lua, "regions"); // pushes key onto stack
  lua_gettable(lua, 1); // replaces key with value
  if (!lua_isnil(lua, -1))
  {
    regions = string_ptr_unordered_map_new(); 
    int t = lua_absindex(lua, -1);

    // Make sure regions is a table.
    if (!lua_istable(lua, t))
      return luaL_error(lua, block_usage);

    // Traverse the regions table.
    lua_pushnil(lua);
    while (lua_next(lua, t))
    {
      int key = -2, val = -1;
      if (!lua_isstring(lua, key) || !lua_isscalarfunction(lua, val))
        return luaL_error(lua, "regions table must map region names to indicator functions.");

      // Add this entry to our regions table.
      const char* name = lua_tostring(lua, key);
      st_func_t* indicator = lua_toscalarfunction(lua, val);
      string_ptr_unordered_map_insert_with_k_dtor(regions, 
                                                  string_dup(name), 
                                                  indicator, 
                                                  string_free);

      lua_pop(lua, 1); // removes value from stack.
    }
    lua_pop(lua, 1); // removes 'regions' argument from stack.
  }
  
  // Create the structured grid using a grid factory.
  int nx = num_cells[0]/patch_size[0];
  int ny = num_cells[1]/patch_size[1];
  int nz = num_cells[2]/patch_size[2];
  str_grid_factory_t* factory = str_grid_factory_new(MPI_COMM_WORLD, patch_size[0], patch_size[1], patch_size[2]);
  str_grid_t* grid = str_grid_factory_brick(factory, nx, ny, nz, NULL);
  factory = NULL;

  // Stash the mapping into a property in the grid if we have it.
  if (mapping != NULL)
    str_grid_set_property(grid, "mapping", mapping, NULL);

  // Adorn the grid with our regions table if we have it.
  if (regions != NULL)
    str_grid_set_property(grid, "regions", regions, NULL);

  lua_pushuserdefined(lua, grid, structured_grid_type_code, DTOR(str_grid_free));
  return 1;
}

static const char* robin_bc_usage = 
  "bc = structured_grids.robin_bc(grid, boundary, A, B, C, component = 0)\n"
  "  Creates an object that enforces a Robin boundary condition of the form\n"
  "    A * U + B * dU/dn = C\n"
  "  on a component of the solution U at the given boundary on the given grid.\n"
  "  Boundary should be one of 'x1', 'x2', 'y1', 'y2', 'z1', or 'z2'.";

static int robin_bc(lua_State* lua)
{
  // Check the number of arguments.
  int num_args = lua_gettop(lua);
  if (((num_args != 5) && (num_args != 6)) || 
      !lua_isuserdefined(lua, 1, structured_grid_type_code) || 
      !lua_isstring(lua, 2) || !lua_isnumber(lua, 3) || 
      !lua_isnumber(lua, 4) || !lua_isnumber(lua, 5))
    return luaL_error(lua, robin_bc_usage);

  // Get the arguments.
  str_grid_t* grid = lua_touserdefined(lua, 1, structured_grid_type_code);
  if (grid == NULL)
    return luaL_error(lua, robin_bc_usage);
  const char* boundary = lua_tostring(lua, 2);
  const char* boundaries[] = {"x1", "x2", "y1", "y2", "z1", "z2"};
  int b_index = string_find_in_list(boundary, boundaries, false);
  if (b_index == -1)
    return luaL_error(lua, robin_bc_usage);
  real_t A = (real_t)lua_tonumber(lua, 3);
  real_t B = (real_t)lua_tonumber(lua, 4);
  real_t C = (real_t)lua_tonumber(lua, 5);
  int component = 0;
  if (num_args == 6)
  {
    if (!lua_isinteger(lua, 6))
      luaL_error(lua, robin_bc_usage);
    component = (int)lua_tointeger(lua, 6);
  }

  // Extract the grid spacing parameter h in the normal direction of the boundary.
  real_t h = grid_spacing(grid, b_index);

  // Construct the str_grid_patch_filler and push it into the interpreter.
  str_grid_patch_filler_t* bc = robin_bc_str_grid_patch_filler_new(A, B, C, h, component, b_index);
  lua_pushuserdefined(lua, bc, patch_filler_type_code, NULL);

  return 1;
}

static const char* dirichlet_bc_usage = 
  "bc = structured_grids.dirichlet_bc(grid, boundary, F, component = 0)\n"
  "  Creates an object that enforces a Dirichlet boundary condition of the form\n"
  "    U = F\n"
  "  on a component of the solution U at the given boundary on the given grid.\n"
  "  Boundary should be one of 'x1', 'x2', 'y1', 'y2', 'z1', or 'z2'.";

static int dirichlet_bc(lua_State* lua)
{
  // Check the number of arguments.
  int num_args = lua_gettop(lua);
  if (((num_args != 3) && (num_args != 4)) || 
      !lua_isuserdefined(lua, 1, structured_grid_type_code) || 
      !lua_isstring(lua, 2) || !lua_isnumber(lua, 3))
    return luaL_error(lua, dirichlet_bc_usage);

  // Get the arguments.
  str_grid_t* grid = lua_touserdefined(lua, 1, structured_grid_type_code);
  if (grid == NULL)
    return luaL_error(lua, dirichlet_bc_usage);
  const char* boundary = lua_tostring(lua, 2);
  const char* boundaries[] = {"x1", "x2", "y1", "y2", "z1", "z2"};
  int b_index = string_find_in_list(boundary, boundaries, false);
  if (b_index == -1)
    return luaL_error(lua, dirichlet_bc_usage);
  real_t F = (real_t)lua_tonumber(lua, 3);
  int component = 0;
  if (num_args == 4)
  {
    if (!lua_isinteger(lua, 4))
      luaL_error(lua, dirichlet_bc_usage);
    component = (int)lua_tointeger(lua, 4);
  }

  // Construct the str_grid_patch_filler and push it into the interpreter.
  str_grid_patch_filler_t* bc = dirichlet_bc_str_grid_patch_filler_new(F, component, b_index);
  lua_pushuserdefined(lua, bc, patch_filler_type_code, NULL);

  return 1;
}

static const char* neumann_bc_usage = 
  "bc = structured_grids.neumann_bc(grid, boundary, A, B, component = 0)\n"
  "  Creates an object that enforces a Neumann boundary condition of the form\n"
  "    A * dU/dn = B\n"
  "  on a component of the solution U at the given boundary on the given grid.\n"
  "  Boundary should be one of 'x1', 'x2', 'y1', 'y2', 'z1', or 'z2'.";

static int neumann_bc(lua_State* lua)
{
  // Check the number of arguments.
  int num_args = lua_gettop(lua);
  if (((num_args != 4) && (num_args != 5)) || 
      !lua_isuserdefined(lua, 1, structured_grid_type_code) || 
      !lua_isstring(lua, 2) || !lua_isnumber(lua, 3) || !lua_isnumber(lua, 4))
    return luaL_error(lua, neumann_bc_usage);

  // Get the arguments.
  str_grid_t* grid = lua_touserdefined(lua, 1, structured_grid_type_code);
  if (grid == NULL)
    return luaL_error(lua, neumann_bc_usage);
  const char* boundary = lua_tostring(lua, 2);
  const char* boundaries[] = {"x1", "x2", "y1", "y2", "z1", "z2"};
  int b_index = string_find_in_list(boundary, boundaries, false);
  if (b_index == -1)
    return luaL_error(lua, neumann_bc_usage);
  real_t A = (real_t)lua_tonumber(lua, 3);
  real_t B = (real_t)lua_tonumber(lua, 4);
  int component = 0;
  if (num_args == 5)
  {
    if (!lua_isinteger(lua, 5))
      luaL_error(lua, neumann_bc_usage);
    component = (int)lua_tointeger(lua, 5);
  }

  // Extract the grid spacing parameter h in the normal direction of the boundary.
  real_t h = grid_spacing(grid, b_index);

  // Construct the str_grid_patch_filler and push it into the interpreter.
  str_grid_patch_filler_t* bc = neumann_bc_str_grid_patch_filler_new(A, B, h, component, b_index);
  lua_pushuserdefined(lua, bc, patch_filler_type_code, NULL);

  return 1;
}

void interpreter_register_polyamri_functions(interpreter_t* interp);
void interpreter_register_polyamri_functions(interpreter_t* interp)
{
  structured_grid_type_code = interpreter_new_user_defined_type_code(interp);
  patch_filler_type_code = interpreter_new_user_defined_type_code(interp);

  interpreter_register_global_table(interp, "structured_grids", NULL);
  interpreter_register_global_method(interp, "structured_grids", "block", block, docstring_from_string(block_usage));
  interpreter_register_global_method(interp, "structured_grids", "robin_bc", robin_bc, docstring_from_string(robin_bc_usage));
  interpreter_register_global_method(interp, "structured_grids", "dirichlet_bc", dirichlet_bc, docstring_from_string(dirichlet_bc_usage));
  interpreter_register_global_method(interp, "structured_grids", "neumann_bc", neumann_bc, docstring_from_string(neumann_bc_usage));
}

str_grid_t* interpreter_get_str_grid(interpreter_t* interp, const char* name);
str_grid_t* interpreter_get_str_grid(interpreter_t* interp, const char* name)
{
  return interpreter_get_user_defined(interp, name, structured_grid_type_code);
}

str_grid_patch_filler_t* interpreter_get_str_grid_patch_filler(interpreter_t* interp, const char* name);
str_grid_patch_filler_t* interpreter_get_str_grid_patch_filler(interpreter_t* interp, const char* name)
{
  return interpreter_get_user_defined(interp, name, patch_filler_type_code);
}
