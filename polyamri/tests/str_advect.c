// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/polyamri.h"
#include "polyamri/str_grid.h"
#include "polyamri/str_grid_patch_filler.h"
#include "polyamri/str_grid_cell_data.h"
#include "polyamri/str_grid_face_data.h"
#include "polyamri/str_ode_integrator.h"
#include "polyamri/interpreter_register_polyamri_functions.h"
#include "polyamri/grid_to_bbox_coord_mapping.h"
#include "model/model.h"
#include "integrators/ark_ode_integrator.h"

//------------------------------------------------------------------------
//                  Lua functions specific to this model
//------------------------------------------------------------------------
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static const char* rigid_body_rotation_usage = 
  "V = rigid_body_rotation{origin = (x0, y0, z0), angular_velocity = (wx, wy, wz)}\n"
  "  Creates a velocity field representing a rigid body rotating in 3D space.\n"
  "  Arguments are:\n"
  "    origin           - Point about which the body rotates.\n"
  "    angular_velocity - Angular velocity (pseudo)vector describing the rotation.";

typedef struct 
{
  point_t x0; // center of rotation
  vector_t omega; // Angular velocity pseudovector.
} sbr_t;

static void sbr_eval(void* context, point_t* x, real_t t, real_t* val)
{
  sbr_t* sbr = context;

  // Displacement vector r.
  vector_t r;
  point_displacement(&sbr->x0, x, &r);

  // V = dr/dt = omega x r.
  vector_cross(&sbr->omega, &r, (vector_t*)val);
}

static int rigid_body_rotation(lua_State* lua)
{
  // Check the number of arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_istable(lua, 1))
    return luaL_error(lua, rigid_body_rotation_usage);

  // Get the origin if it's given.
  lua_pushstring(lua, "origin"); // pushes key onto stack
  lua_gettable(lua, 1); // replaces key with value
  if (!lua_isnil(lua, -1) && !lua_ispoint(lua, -1))
    return luaL_error(lua, rigid_body_rotation_usage);
  point_t* x0 = lua_topoint(lua, -1);
  lua_pop(lua, 1);

  // Get the angular velocity vector.
  lua_pushstring(lua, "angular_velocity"); // pushes key onto stack
  lua_gettable(lua, 1); // replaces key with value
  if (lua_isnil(lua, -1) || !lua_isvector(lua, -1))
    return luaL_error(lua, rigid_body_rotation_usage);
  vector_t* omega = lua_tovector(lua, -1);
  lua_pop(lua, 1);

  // Set up the rigid body rotation thingy.
  sbr_t* sbr = polymec_malloc(sizeof(sbr_t));
  if (x0 != NULL)
    sbr->x0 = *x0;
  else
    sbr->x0.x = sbr->x0.y = sbr->x0.z = 0.0;
  sbr->omega = *omega;

  // Create and return the function.
  char name[1025];
  snprintf(name, 1024, "Rigid body rotation(x0 = (%g, %g, %g), omega = (%g, %g, %g))",
           sbr->x0.x, sbr->x0.y, sbr->x0.z, omega->x, omega->y, omega->z);
  st_func_vtable vtable = {.eval = sbr_eval, .dtor = polymec_free};
  st_func_t* func = st_func_new(name, sbr, vtable, ST_FUNC_HETEROGENEOUS, ST_FUNC_CONSTANT, 3);

  lua_pushvectorfunction(lua, func);
  return 1;
}

static void interpreter_register_advect_functions(interpreter_t* interp)
{
  interpreter_register_function(interp, "rigid_body_rotation", rigid_body_rotation, 
    docstring_from_string(rigid_body_rotation_usage));
}

//------------------------------------------------------------------------
//                    End str_advect-specific Lua functions
//------------------------------------------------------------------------

typedef struct
{
  MPI_Comm comm;

  // Grid and computational domain.
  str_grid_t* grid;
  coord_mapping_t* domain;
  real_t dx, dy, dz;

  // State information.
  str_grid_cell_data_t* U;

  // Initial state.
  st_func_t* U0;

  // Velocity field.
  st_func_t* V;
  real_t V_max;

  // Time integrator and work vectors.
  str_ode_integrator_t* integ;
  str_grid_cell_data_t* U_work;
  str_grid_cell_data_t* dUdt_work;
  str_grid_face_data_t* F_work;
} advect_t;

static const char advect_desc[] = "PolyAMRI advection (str_advect)\n"
  "This model demonstrates a parallel implementation of a finite volume\n"
  "discretization of the linear advection equation on a uniform grid,\n"
  "integrated in time using the method of lines.\n";

static void advect_read_input(void* context, 
                              interpreter_t* interp, 
                              options_t* options)
{
  advect_t* adv = context;

  adv->U0 = interpreter_get_scalar_function(interp, "initial_state");
  adv->V = interpreter_get_vector_function(interp, "velocity");
  adv->grid = interpreter_get_str_grid(interp, "grid");

  // Computational domain can either be a bounding box or a coordinate 
  // mapping.
  bbox_t* bbox = interpreter_get_bbox(interp, "domain");
  if (bbox != NULL)
    adv->domain = grid_to_bbox_coord_mapping_new(bbox);
  else
  {
    coord_mapping_t* mapping = interpreter_get_coord_mapping(interp, "domain");
    if (mapping != NULL)
      adv->domain = mapping;
    else
      polymec_error("domain must be a bounding box or a coordinate mapping.");
  }
}

static real_t advect_max_dt(void* context, real_t t, char* reason)
{
  advect_t* adv = context;
  real_t dx = MIN(adv->dx, MIN(adv->dy, adv->dz));

  if ((adv->V_max == -FLT_MAX) || !st_func_is_constant(adv->V))
  {
    // We compute the maximum allowable timestep using the CFL condition, 
    // measuring the velocity at each of the face centers at the given time and 
    // taking the max.
    int pos = 0, ip, jp, kp;
    str_grid_patch_t* patch;
    bbox_t bbox;
    while (str_grid_cell_data_next_patch(adv->U, &pos, &ip, &jp, &kp, &patch, &bbox))
    {
      point_t x;
      vector_t v;

      // x face velocities.
      for (int i = patch->i1; i <= patch->i2; ++i)
      {
        x.x = bbox.x1 + i * adv->dx;
        for (int j = patch->j1; j < patch->j2; ++j)
        {
          x.y = bbox.y1 + (j + 0.5) * adv->dy;
          for (int k = patch->k1; k < patch->k2; ++k)
          {
            x.z = bbox.z1 + (k + 0.5) * adv->dz;
            st_func_eval(adv->V, &x, t, (real_t*)&v);
            adv->V_max = MAX(adv->V_max, vector_mag(&v));
          }
        }
      }

      // y face velocities.
      for (int i = patch->i1; i < patch->i2; ++i)
      {
        x.x = bbox.x1 + (i + 0.5) * adv->dx;
        for (int j = patch->j1; j <= patch->j2; ++j)
        {
          x.y = bbox.y1 + j * adv->dy;
          for (int k = patch->k1; k < patch->k2; ++k)
          {
            x.z = bbox.z1 + (k + 0.5) * adv->dz;
            st_func_eval(adv->V, &x, t, (real_t*)&v);
            adv->V_max = MAX(adv->V_max, vector_mag(&v));
          }
        }
      }

      // z face velocities.
      for (int i = patch->i1; i < patch->i2; ++i)
      {
        x.x = bbox.x1 + (i + 0.5) * adv->dx;
        for (int j = patch->j1; j < patch->j2; ++j)
        {
          x.y = bbox.y1 + (j + 0.5) * adv->dy;
          for (int k = patch->k1; k <= patch->k2; ++k)
          {
            x.z = bbox.z1 + k * adv->dz;
            st_func_eval(adv->V, &x, t, (real_t*)&v);
            adv->V_max = MAX(adv->V_max, vector_mag(&v));
          }
        }
      }
    }
  }
  snprintf(reason, POLYMEC_MODEL_MAXDT_REASON_SIZE, 
           "CFL constraint (V_max = %g, dx = %g)", adv->V_max, dx);

  if (adv->V_max > 0.0)
    return dx / adv->V_max;
  else 
    return FLT_MAX;
}

static real_t ark_max_dt_wrapper(void* context, real_t t, real_t* x)
{
  char reason[POLYMEC_MODEL_MAXDT_REASON_SIZE];
  return advect_max_dt(context, t, reason);
}

static void advect_clear(advect_t* adv)
{
  if (adv->integ != NULL)
    str_ode_integrator_free(adv->integ);
  if (adv->U != NULL)
    str_grid_cell_data_free(adv->U);
  if (adv->F_work != NULL)
    str_grid_face_data_free(adv->F_work);
  if (adv->U_work != NULL)
    str_grid_cell_data_free(adv->U_work);
  if (adv->dUdt_work != NULL)
    str_grid_cell_data_free(adv->dUdt_work);
}

// Here's the right-hand side for the ARK integrator.
static int advect_rhs(void* context, real_t t, real_t* X, real_t* dXdt)
{
  advect_t* adv = context;

  // Point our "work vector" at the given buffer.
  str_grid_cell_data_set_buffer(adv->U_work, X, false);

  // Make sure ghost cells are filled.
  str_grid_fill_ghost_cells(adv->grid, adv->U_work);

  // Now go over each cell and compute fluxes.
  int pos, ip, jp, kp;
  str_grid_patch_t* face_patch;
  bbox_t bbox;
  
  // x fluxes.
  pos = 0;
  real_t dx = adv->dx, dy = adv->dy, dz = adv->dz;
  while (str_grid_face_data_next_x_patch(adv->F_work, &pos, &ip, &jp, &kp, &face_patch, &bbox))
  {
    str_grid_patch_t* cell_patch = str_grid_cell_data_patch(adv->U_work, ip, jp, kp);
    DECLARE_STR_GRID_PATCH_ARRAY(F, face_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(U, cell_patch);
    point_t x_low, x_high;
    for (int i = face_patch->i1; i < face_patch->i2; ++i)
    {
      x_low.x = bbox.x1 + i*dx;
      x_high.x = x_low.x + dx;
      for (int j = face_patch->j1; j < face_patch->j2; ++j)
      {
        x_low.y = x_high.y = bbox.y1 + (j+0.5) * dy;
        for (int k = face_patch->k1; k < face_patch->k2; ++k)
        {
          x_low.z = x_high.z = bbox.z1 + (k+0.5) * dz;

          // Get the velocity at the low and high faces.
          vector_t v_low, v_high;
          st_func_eval(adv->V, &x_low, t, (real_t*)&v_low);
          st_func_eval(adv->V, &x_high, t, (real_t*)&v_high);

          // Get the upwinded values of U at the low and high faces.
          // Note that since i, j, k are face indices, we need to 
          // translate them to cell indices.
          int ic, jc, kc; // <-- cell indices
          str_grid_patch_translate_indices(face_patch, i, j, k,
                                           cell_patch, &ic, &jc, &kc);
          real_t U_low = (v_low.x >= 0.0) ? U[ic-1][jc][kc][0] : U[ic][jc][kc][0];
          real_t U_high = (v_high.x >= 0.0) ? U[ic][jc][kc][0] : U[ic+1][jc][kc][0];

          // Compute the low and high fluxes.
          F[i][j][k][0]   = U_low * v_low.x;
          F[i+1][j][k][0] = U_high * v_high.x;
        }
      }
    }
  }

  // y fluxes.
  pos = 0;
  while (str_grid_face_data_next_y_patch(adv->F_work, &pos, &ip, &jp, &kp, &face_patch, &bbox))
  {
    str_grid_patch_t* cell_patch = str_grid_cell_data_patch(adv->U_work, ip, jp, kp);
    DECLARE_STR_GRID_PATCH_ARRAY(F, face_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(U, cell_patch);
    point_t x_low, x_high;
    for (int i = face_patch->i1; i < face_patch->i2; ++i)
    {
      x_low.x = x_high.x = bbox.x1 + (i+0.5) * dx;
      for (int j = face_patch->j1; j < face_patch->j2; ++j)
      {
        x_low.y = bbox.y1 + j * dy;
        x_high.y = x_low.y + dy;
        for (int k = face_patch->k1; k < face_patch->k2; ++k)
        {
          x_low.z = x_high.z = bbox.z1 + (k+0.5) * dz;

          // Get the velocity at the low and high faces.
          vector_t v_low, v_high;
          st_func_eval(adv->V, &x_low, t, (real_t*)&v_low);
          st_func_eval(adv->V, &x_high, t, (real_t*)&v_high);

          // Get the upwinded values of U at the low and high faces.
          // Note that since i, j, k are face indices, we need to 
          // translate them to cell indices.
          int ic, jc, kc; // <-- cell indices
          str_grid_patch_translate_indices(face_patch, i, j, k,
                                           cell_patch, &ic, &jc, &kc);
          real_t U_low = (v_low.y >= 0.0) ? U[ic][jc-1][kc][0] : U[ic][jc][kc][0];
          real_t U_high = (v_high.y >= 0.0) ? U[ic][jc][kc][0] : U[ic][jc+1][kc][0];

          // Compute the low and high fluxes.
          F[i][j][k][0]   = U_low * v_low.y;
          F[i][j+1][k][0] = U_high * v_high.y;
        }
      }
    }
  }

  // z fluxes.
  pos = 0;
  while (str_grid_face_data_next_z_patch(adv->F_work, &pos, &ip, &jp, &kp, &face_patch, &bbox))
  {
    str_grid_patch_t* cell_patch = str_grid_cell_data_patch(adv->U_work, ip, jp, kp);
    DECLARE_STR_GRID_PATCH_ARRAY(F, face_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(U, cell_patch);
    point_t x_low, x_high;
    for (int i = face_patch->i1; i < face_patch->i2; ++i)
    {
      x_low.x = x_high.x = bbox.x1 + (i+0.5)*dx;
      for (int j = face_patch->j1; j < face_patch->j2; ++j)
      {
        x_low.y = x_high.y = bbox.y1 + (j+0.5)*dy;
        for (int k = face_patch->k1; k < face_patch->k2; ++k)
        {
          x_low.z = bbox.z1 + k*dz;
          x_high.z = x_low.z + dz;

          // Get the velocity at the low and high faces.
          vector_t v_low, v_high;
          st_func_eval(adv->V, &x_low, t, (real_t*)&v_low);
          st_func_eval(adv->V, &x_high, t, (real_t*)&v_high);

          // Get the upwinded values of U at the low and high faces.
          // Note that since i, j, k are face indices, we need to 
          // translate them to cell indices.
          int ic, jc, kc; // <-- cell indices
          str_grid_patch_translate_indices(face_patch, i, j, k,
                                           cell_patch, &ic, &jc, &kc);
          real_t U_low = (v_low.z >= 0.0) ? U[ic][jc][kc-1][0] : U[ic][jc][kc][0];
          real_t U_high = (v_high.z >= 0.0) ? U[ic][jc][kc][0] : U[ic][jc][kc+1][0];

          // Compute the low and high fluxes.
          F[i][j][k][0]   = U_low * v_low.z;
          F[i][j][k+1][0] = U_high * v_high.z;
        }
      }
    }
  }

  // Now compute dU/dt.
  str_grid_cell_data_set_buffer(adv->dUdt_work, dXdt, false);

  real_t V = dx * dy * dz;
  real_t Ax = dy * dz, Ay = dz * dx, Az = dx * dy;
  pos = 0;
  str_grid_patch_t* rhs_patch;
  while (str_grid_cell_data_next_patch(adv->dUdt_work, &pos, &ip, &jp, &kp, &rhs_patch, NULL))
  {
    str_grid_patch_t* x_face_patch = str_grid_face_data_x_patch(adv->F_work, ip, jp, kp);
    str_grid_patch_t* y_face_patch = str_grid_face_data_x_patch(adv->F_work, ip, jp, kp);
    str_grid_patch_t* z_face_patch = str_grid_face_data_x_patch(adv->F_work, ip, jp, kp);
    DECLARE_STR_GRID_PATCH_ARRAY(dUdt, rhs_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Fx, x_face_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Fy, y_face_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Fz, z_face_patch);

    for (int i = rhs_patch->i1; i < rhs_patch->i2; ++i)
    {
      for (int j = rhs_patch->j1; j < rhs_patch->j2; ++j)
      {
        for (int k = rhs_patch->k1; k < rhs_patch->k2; ++k)
        {
          // dUdt and the fluxes happen to share similar indexing, since 
          // they have no ghost values. The only difference is that the 
          // fluxes have an extra datum in their preferred direction.
          // But no index translation is needed here.
          real_t div_F = Ax * (Fx[i+1][j][k][0] - Fx[i][j][k][0]) + 
                         Ay * (Fy[i][j+1][k][0] - Fy[i][j][k][0]) + 
                         Az * (Fz[i][j][k+1][0] - Fz[i][j][k][0]);
          dUdt[i][j][k][0] = -div_F / V;
        }
      }
    }
  }

  return 0;
}

static void advect_setup(advect_t* adv)
{
  ASSERT(st_func_num_comp(adv->U0) == 1);
  ASSERT(st_func_num_comp(adv->V) == 3);

  advect_clear(adv);

  // Here we set up the machinery to fill ghost cells in the grid.
  str_grid_patch_filler_t* fill_from_east = copy_str_grid_patch_filler_new(STR_GRID_PATCH_X1_BOUNDARY, STR_GRID_PATCH_X2_BOUNDARY);
  str_grid_patch_filler_t* fill_from_west = copy_str_grid_patch_filler_new(STR_GRID_PATCH_X2_BOUNDARY, STR_GRID_PATCH_X1_BOUNDARY);
  str_grid_patch_filler_t* fill_from_south = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Y2_BOUNDARY, STR_GRID_PATCH_Y1_BOUNDARY);
  str_grid_patch_filler_t* fill_from_north = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Y1_BOUNDARY, STR_GRID_PATCH_Y2_BOUNDARY);
  str_grid_patch_filler_t* fill_from_above = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Z1_BOUNDARY, STR_GRID_PATCH_Z2_BOUNDARY);
  str_grid_patch_filler_t* fill_from_below = copy_str_grid_patch_filler_new(STR_GRID_PATCH_Z2_BOUNDARY, STR_GRID_PATCH_Z1_BOUNDARY);

  // We impose a zero flux on all boundaries, on the assumption that the 
  // solution is zero at the boundary.
  str_grid_patch_filler_t* zero_flux_x1 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_X1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_x2 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_X2_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_y1 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Y1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_y2 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Y2_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_z1 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Z1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_z2 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Z2_BOUNDARY);

  int npx, npy, npz;
  str_grid_get_extents(adv->grid, &npx, &npy, &npz);

  int pos = 0, ip, jp, kp;
  while (str_grid_next_patch(adv->grid, &pos, &ip, &jp, &kp))
  {
    if (ip > 0)
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, fill_from_west);
    else
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, zero_flux_x1);
    if (ip < npx-1)
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, fill_from_east);
    else
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, zero_flux_x2);
    if (jp > 0)
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, fill_from_south);
    else
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, zero_flux_y1);
    if (jp < npy-1)
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, fill_from_north);
    else
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, zero_flux_y2);
    if (kp > 0)
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, fill_from_below);
    else
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, zero_flux_z1);
    if (kp < npz-1)
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, fill_from_above);
    else
      str_grid_append_patch_filler(adv->grid, ip, jp, kp, zero_flux_z2);
  }

  // Set the grid spacings.
  int nx, ny, nz;
  str_grid_get_patch_size(adv->grid, &nx, &ny, &nz);
  int Nx = npx * nx;
  int Ny = npy * ny;
  int Nz = npz * nz;
  adv->dx = 1.0 / Nx;
  adv->dy = 1.0 / Ny;
  adv->dz = 1.0 / Nz;

  // Create the state vector.
  int num_comps = 1;
  int num_ghost_layers = 1;
  adv->U = str_grid_cell_data_new(adv->grid, num_comps, num_ghost_layers);

  // Create work "vectors."
  adv->U_work = str_grid_cell_data_with_buffer(adv->grid, 1, 1, NULL);
  adv->dUdt_work = str_grid_cell_data_with_buffer(adv->grid, 1, 0, NULL);
  adv->F_work = str_grid_face_data_new(adv->grid, 1);

  // Create the ARK integrator.
  int num_local_cells = str_grid_cell_data_num_cells(adv->U, false);
  int total_num_cells = str_grid_cell_data_num_cells(adv->U, true);
  int num_ghost_cells = total_num_cells - num_local_cells;
  int N_local = num_comps * num_local_cells;
  int N_remote = num_comps * num_ghost_cells;
  ode_integrator_t* ark2 = explicit_ark_ode_integrator_new(2, 
                                                           MPI_COMM_WORLD,
                                                           N_local,
                                                           N_remote,
                                                           adv,
                                                           advect_rhs,
                                                           ark_max_dt_wrapper,
                                                           NULL);
  adv->integ = str_ode_integrator_new(ark2);

  // Reset the max velocity.
  adv->V_max = -FLT_MAX;
}

static void advect_init(void* context, real_t t)
{
  advect_t* adv = context;
  advect_setup(adv);

  // Initialize the state.
  int pos = 0, ip, jp, kp;
  str_grid_patch_t* patch;
  bbox_t bbox;
  while (str_grid_cell_data_next_patch(adv->U, &pos, &ip, &jp, &kp, &patch, &bbox))
  {
    DECLARE_STR_GRID_PATCH_ARRAY(U, patch);
    point_t x;
    for (int i = patch->i1; i < patch->i2; ++i)
    {
      x.x = bbox.x1 + (i + 0.5) * adv->dx;
      for (int j = patch->j1; j < patch->j2; ++j)
      {
        x.y = bbox.y1 + (j + 0.5) * adv->dy;
        for (int k = patch->k1; k < patch->k2; ++k)
        {
          x.z = bbox.z1 + (k + 0.5) * adv->dz;
          st_func_eval(adv->U0, &x, t, &U[i][j][k][0]);
        }
      }
    }
  }
}

static real_t advect_advance(void* context, real_t max_dt, real_t t)
{
  advect_t* adv = context;

  // Integrate!
  real_t t2 = t;
  if (!str_ode_integrator_step(adv->integ, max_dt, &t2, adv->U))
    polymec_error("advect_advance: Integration failed at t = %g.", t);

  return t2 - t;
}

static void advect_load(void* context, const char* file_prefix, const char* directory, real_t* t, int step)
{
#if 0
  advect_t* adv = context;

  silo_file_t* silo = silo_file_open(MPI_COMM_WORLD, prefix, directory, 0, step, t);
  if (advect->grid != NULL)
    mesh_free(mod->mesh);
  mod->mesh = silo_file_read_mesh(silo, "mesh");
  real_t* X = silo_file_read_scalar_cell_field(silo, "X", "mesh", NULL);
  silo_file_close(silo);

  // Set everything up and copy the saved solution.
  advect_setup(adv);
  memcpy(mod->X, X, sizeof(real_t) * mod->mesh->num_cells);
  polymec_free(X);
#endif
}

static void advect_save(void* context, const char* file_prefix, const char* directory, real_t t, int step)
{
#if 0
  advect_t* adv = context;

  int rank;
  MPI_Comm_rank(adv->grid->comm, &rank);

  silo_file_t* silo = silo_file_new(adv->comm, file_prefix, directory, 1, 0, step, t);
  silo_file_write_mesh(silo, "grid", adv->grid);
  silo_file_write_scalar_cell_field(silo, "life", "grid", adv->state, NULL);
  real_t* rank_field = polymec_malloc(sizeof(real_t) * adv->grid->num_cells);
  for (int i = 0; i < adv->grid->num_cells; ++i)
    rank_field[i] = 1.0 * rank;
  silo_file_write_scalar_cell_field(silo, "rank", "grid", rank_field, NULL);
  polymec_free(rank_field);
  silo_file_close(silo);
#endif
}

static void advect_finalize(void* context, int step, real_t t)
{
//  advect_t* adv = context;
}

static void advect_dtor(void* context)
{
  advect_t* adv = context;

  advect_clear(adv);
  str_grid_free(adv->grid);
  polymec_free(adv);
}

static model_t* advect_ctor()
{
  advect_t* adv = polymec_malloc(sizeof(advect_t));
  adv->grid = NULL;
  adv->domain = NULL;
  adv->U0 = NULL;
  adv->V = NULL;
  adv->U = NULL;
  adv->F_work = NULL;
  adv->U_work = NULL;
  adv->dUdt_work = NULL;
  adv->integ = NULL;
  model_vtable vtable = {.read_input = advect_read_input,
                         .init = advect_init,
                         .max_dt = advect_max_dt,
                         .advance = advect_advance,
                         .load = advect_load,
                         .save = advect_save,
                         .finalize = advect_finalize,
                         .dtor = advect_dtor};
  docstring_t* advect_doc = docstring_from_string(advect_desc);
  model_t* model = model_new("str_advect", adv, vtable, advect_doc, MODEL_MPI);

  // Register polyamri- and advection-specific functions.
  interpreter_t* interp = model_interpreter(model);
  interpreter_register_polyamri_functions(interp);
  interpreter_register_advect_functions(interp);

  return model;
}

// Main program.
int main(int argc, char* argv[])
{
  return model_main("str_advect", advect_ctor, argc, argv);
}
