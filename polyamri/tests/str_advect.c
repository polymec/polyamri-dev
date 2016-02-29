// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/polyamri.h"
#include "core/string_utils.h"
#include "polyamri/str_grid.h"
#include "polyamri/str_grid_cell_filler_factory.h"
#include "polyamri/str_grid_cell_data.h"
#include "polyamri/str_grid_face_data.h"
#include "polyamri/ark_str_ode_integrator.h"
#include "polyamri/interpreter_register_polyamri_functions.h"
#include "polyamri/grid_to_bbox_coord_mapping.h"
#include "polyamri/silo_file.h"
#include "model/model.h"

// Fortran-style SIGN function
#define SIGN_VAL(A, B) ((B >= 0.0) ? ABS(A) : -ABS(A))

// Uncomment this to use first-order upwinding.
//#define ADVECT_FIRST_ORDER

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

// Flux limiters.
static real_t no_limiter(real_t ratio, real_t omega, real_t delta)
{
  return delta;
}

static real_t superbee_limiter(real_t ratio, real_t omega, real_t delta)
{
  real_t phi = 0.0;
  if (ratio >= 0.0)
    phi = 2.0 * ratio;
  if (ratio >= 0.5)
    phi = 1.0;
  if (ratio >= 1.0)
  {
    real_t denom = 1.0 - omega + (1.0 + omega) * ratio;
    real_t phi_r = 2.0 / denom;
    phi = MIN(2.0, MIN(phi_r, ratio));
  }
  return phi * delta;
}

static real_t van_leer_limiter(real_t ratio, real_t omega, real_t delta)
{
  real_t phi = 0.0;
  if (ratio >= 0.0)
  {
    real_t denom = 1.0 - omega + (1.0 + omega) * ratio;
    real_t phi_r = 2.0 / denom;
    phi = MIN(2.0 * ratio / (1.0 + ratio), phi_r);
  }
  return phi * delta;
}

static real_t van_albada_limiter(real_t ratio, real_t omega, real_t delta)
{
  real_t phi = 0.0;
  if (ratio >= 0.0)
  {
    real_t denom = 1.0 - omega + (1.0 + omega) * ratio;
    real_t phi_r = 2.0 / denom;
    phi = MIN(ratio * (1.0 + ratio) / (1.0 + ratio*ratio), phi_r);
  }
  return phi * delta;
}

static real_t minmod_limiter(real_t ratio, real_t omega, real_t delta)
{
  real_t phi = 0.0;
  if (ratio >= 0.0)
    phi = ratio;
  if (ratio >= 1.0)
  {
    real_t denom = 2.0 * (1.0 - omega + (1.0 + omega) * ratio);
    real_t phi_r = 4.0 / denom;
    phi = MIN(1.0, phi_r);
  }
  return phi * delta;
}

typedef struct
{
  MPI_Comm comm;

  // Grid and computational domain.
  str_grid_t* grid;
  coord_mapping_t* mapping;
  coord_mapping_t* inv_mapping;
  real_t dx, dy, dz;

  // Ghost filler.
  str_grid_cell_filler_t* ghost_filler;

  // State information.
  str_grid_cell_data_t* U;

  // Initial state.
  st_func_t* U0;

  // Velocity field.
  st_func_t* V_func;
  real_t V_max;

  // Flux limiter function.
  real_t (*limit_delta)(real_t ratio, real_t omega, real_t delta);

  // Time integrator and workspace.
  str_ode_integrator_t* integ;
  str_grid_cell_data_t* UL;
  str_grid_cell_data_t* UH;
  str_grid_face_data_t* V;
  str_grid_face_data_t* F;
  char max_dt_reason[POLYMEC_MODEL_MAXDT_REASON_SIZE];
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
  adv->V_func = interpreter_get_vector_function(interp, "velocity");
  adv->grid = interpreter_get_str_grid(interp, "grid");

  // Computational domain can either be a bounding box or a coordinate 
  // mapping.
  bbox_t* bbox = interpreter_get_bbox(interp, "domain");
  if (bbox != NULL)
    adv->mapping = grid_to_bbox_coord_mapping_new(bbox);
  else
  {
    coord_mapping_t* mapping = interpreter_get_coord_mapping(interp, "domain");
    if (mapping != NULL)
      adv->mapping = mapping;
    else
      polymec_error("domain must be a bounding box or a coordinate mapping.");
  }
  
  // Make sure the mapping is invertible.
  adv->inv_mapping = coord_mapping_inverse(adv->mapping);
  if (adv->inv_mapping == NULL)
    polymec_error("domain mapping is not invertible! Must be an invertible mapping.");
}

static real_t ark_max_dt(void* context, real_t t, str_grid_cell_data_t* U)
{
  advect_t* adv = context;
  real_t dx = MIN(adv->dx, MIN(adv->dy, adv->dz));

  if ((adv->V_max == -FLT_MAX) || !st_func_is_constant(adv->V_func))
  {
    // We compute the maximum allowable timestep using the CFL condition, 
    // measuring the velocity at each of the face centers at the given time and 
    // taking the max.
    int pos = 0, ip, jp, kp;
    str_grid_patch_t* patch;
    bbox_t bbox;
    #pragma omp parallel firstprivate(pos) private(ip, jp, kp, patch, bbox)
    while (str_grid_cell_data_next_patch(U, &pos, &ip, &jp, &kp, &patch, &bbox))
    {
      point_t x, x1;
      vector_t v, v1;

      // x face velocities.
      for (int i = patch->i1; i <= patch->i2; ++i)
      {
        x1.x = bbox.x1 + i * adv->dx;
        for (int j = patch->j1; j < patch->j2; ++j)
        {
          x1.y = bbox.y1 + (j + 0.5) * adv->dy;
          for (int k = patch->k1; k < patch->k2; ++k)
          {
            x1.z = bbox.z1 + (k + 0.5) * adv->dz;
            coord_mapping_map_point(adv->mapping, &x1, &x);
            st_func_eval(adv->V_func, &x, t, (real_t*)&v);
            coord_mapping_map_vector(adv->inv_mapping, &x, &v, &v1);
            adv->V_max = MAX(adv->V_max, vector_mag(&v1));
          }
        }
      }

      // y face velocities.
      for (int i = patch->i1; i < patch->i2; ++i)
      {
        x1.x = bbox.x1 + (i + 0.5) * adv->dx;
        for (int j = patch->j1; j <= patch->j2; ++j)
        {
          x1.y = bbox.y1 + j * adv->dy;
          for (int k = patch->k1; k < patch->k2; ++k)
          {
            x1.z = bbox.z1 + (k + 0.5) * adv->dz;
            coord_mapping_map_point(adv->mapping, &x1, &x);
            st_func_eval(adv->V_func, &x, t, (real_t*)&v);
            coord_mapping_map_vector(adv->inv_mapping, &x, &v, &v1);
            adv->V_max = MAX(adv->V_max, vector_mag(&v1));
          }
        }
      }

      // z face velocities.
      for (int i = patch->i1; i < patch->i2; ++i)
      {
        x1.x = bbox.x1 + (i + 0.5) * adv->dx;
        for (int j = patch->j1; j < patch->j2; ++j)
        {
          x1.y = bbox.y1 + (j + 0.5) * adv->dy;
          for (int k = patch->k1; k <= patch->k2; ++k)
          {
            x1.z = bbox.z1 + k * adv->dz;
            coord_mapping_map_point(adv->mapping, &x1, &x);
            st_func_eval(adv->V_func, &x, t, (real_t*)&v);
            coord_mapping_map_vector(adv->inv_mapping, &x, &v, &v1);
            adv->V_max = MAX(adv->V_max, vector_mag(&v1));
          }
        }
      }
    }
  }
  snprintf(adv->max_dt_reason, POLYMEC_MODEL_MAXDT_REASON_SIZE-1, 
           "CFL constraint (V_max = %g, dx = %g)", adv->V_max, dx);

  if (adv->V_max > 0.0)
    return dx / adv->V_max;
  else 
    return FLT_MAX;
}

static real_t advect_max_dt(void* context, real_t t, char* reason)
{
  advect_t* adv = context;
  real_t dt = ark_max_dt(context, t, adv->U);
  strcpy(reason, adv->max_dt_reason);
  return dt;
}

static void advect_clear(advect_t* adv)
{
  if (adv->integ != NULL)
    str_ode_integrator_free(adv->integ);
  if (adv->U != NULL)
    str_grid_cell_data_free(adv->U);
  if (adv->UL != NULL)
    str_grid_cell_data_free(adv->UL);
  if (adv->UH != NULL)
    str_grid_cell_data_free(adv->UH);
  if (adv->F != NULL)
    str_grid_face_data_free(adv->F);
  if (adv->V != NULL)
    str_grid_face_data_free(adv->V);
  if (adv->ghost_filler != NULL)
    str_grid_cell_filler_free(adv->ghost_filler);
}

static void extrapolate_U_to_faces(advect_t* adv,
                                   real_t t,
                                   str_grid_cell_data_t* U,
                                   str_grid_cell_data_t* UL,
                                   str_grid_cell_data_t* UH,
                                   str_grid_face_data_t* V)
{
  real_t dx = adv->dx, dy = adv->dy, dz = adv->dz;

  int pos, ip, jp, kp;
  bbox_t bbox;
  str_grid_patch_t* V_patch;

  // Traverse the patches and extrapolate U to x-faces.
  pos = 0;
  #pragma omp parallel firstprivate(pos) private(ip, jp, kp, V_patch, bbox)
  while (str_grid_face_data_next_x_patch(V, &pos, &ip, &jp, &kp, &V_patch, &bbox))
  {
    str_grid_patch_t* U_patch = str_grid_cell_data_patch(U, ip, jp, kp);
    str_grid_patch_t* UL_patch = str_grid_cell_data_patch(UL, ip, jp, kp);
    str_grid_patch_t* UH_patch = str_grid_cell_data_patch(UH, ip, jp, kp);

    // Compute "low" and "high" values of U on the x-faces of the cells.
    DECLARE_STR_GRID_PATCH_ARRAY(U, U_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(UL, UL_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(UH, UH_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Vx, V_patch);

    point_t x_low, x1_low, x_high, x1_high;
    for (int i = V_patch->i1; i < V_patch->i2-1; ++i) 
    {
      x1_low.x = bbox.x1 + i * dx;
      x1_high.x = x1_low.x + dx;
      for (int j = V_patch->j1; j < V_patch->j2; ++j)
      {
        x1_low.y = x1_high.y = bbox.y1 + (j+0.5) * dy;
        for (int k = V_patch->k1; k < V_patch->k2; ++k)
        {
          x1_low.z = x1_high.z = bbox.z1 + (k+0.5) * dz;

          // Get the velocity at the low and high faces, recording 
          // the velocities.
          vector_t v_low, v1_low, v_high, v1_high;
          coord_mapping_map_point(adv->mapping, &x1_low, &x_low);
          st_func_eval(adv->V_func, &x_low, t, (real_t*)&v_low);
          coord_mapping_map_vector(adv->inv_mapping, &x_low, &v_low, &v1_low);
          coord_mapping_map_point(adv->mapping, &x1_high, &x_high);
          st_func_eval(adv->V_func, &x_high, t, (real_t*)&v_high);
          coord_mapping_map_vector(adv->inv_mapping, &x_high, &v_high, &v1_high);
          Vx[i][j][k][0] = v1_low.x;
          Vx[i+1][j][k][0] = v1_high.x;

          // Note that since i, j, k are face indices, we need to 
          // translate them to cell indices.
          int ic, jc, kc; // <-- cell indices
          str_grid_patch_translate_indices(V_patch, i, j, k,
                                           U_patch, &ic, &jc, &kc);
#ifdef ADVECT_FIRST_ORDER 
          // Get the values of U at the low and high faces.
          UL[ic][jc][kc][0] = (v_low.x >= 0.0) ? U[ic-1][jc][kc][0] : U[ic][jc][kc][0];
          UH[ic][jc][kc][0] = (v_high.x >= 0.0) ? U[ic][jc][kc][0] : U[ic+1][jc][kc][0];
#else
          // Use the MUSCL-Hancock method to get a 2nd-order accurate flux.
          real_t U_minus = U[ic-1][jc][kc][0];
          real_t U_0     = U[ic][jc][kc][0];
          real_t U_plus  = U[ic+1][jc][kc][0];
          real_t delta_upwind = (v_low.x >= 0.0) ? U_0 - U_minus : U_plus - U_0;
          real_t delta_loc    = (v_low.x >= 0.0) ? U_plus - U_0 : U_0 - U_minus;

          static const real_t omega = 0.0;
          real_t delta = 0.5 * ((1.0 + omega) * delta_upwind + 
                                (1.0 - omega) * delta_loc);

          // Reset small deltas, preserving sign.
          static const real_t tol = 1e-6;
          if (ABS(delta_upwind) < tol)
            delta_upwind = tol * SIGN_VAL(1.0, delta_upwind);
          if (ABS(delta_loc) < tol)
            delta_loc = tol * SIGN_VAL(1.0, delta_loc);

          // Compute the ratio of slopes and apply the appropriate limiter.
          real_t ratio = delta_upwind/delta_loc;
          delta = adv->limit_delta(ratio, omega, delta);

          // Boundary extrapolated values.
          UL[ic][jc][kc][0] = U_0 - 0.5 * delta;
          UH[ic][jc][kc][0] = U_0 + 0.5 * delta;
#endif
        }
      }
    }
  }

  // Traverse the patches and extrapolate U to y-faces.
  pos = 0;
  #pragma omp parallel firstprivate(pos) private(ip, jp, kp, V_patch, bbox)
  while (str_grid_face_data_next_y_patch(V, &pos, &ip, &jp, &kp, &V_patch, &bbox))
  {
    str_grid_patch_t* U_patch = str_grid_cell_data_patch(U, ip, jp, kp);
    str_grid_patch_t* UL_patch = str_grid_cell_data_patch(UL, ip, jp, kp);
    str_grid_patch_t* UH_patch = str_grid_cell_data_patch(UH, ip, jp, kp);

    // Compute "low" and "high" values of U on the y-faces of the cells.
    DECLARE_STR_GRID_PATCH_ARRAY(U, U_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(UL, UL_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(UH, UH_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Vy, V_patch);

    point_t x_low, x1_low, x_high, x1_high;
    for (int i = V_patch->i1; i < V_patch->i2; ++i) 
    {
      x1_low.x = x1_high.x = bbox.x1 + (i+0.5) * dx;
      for (int j = V_patch->j1; j < V_patch->j2-1; ++j)
      {
        x1_low.y = bbox.y1 + j * dy;
        x1_high.y = x1_low.y + dy;
        for (int k = V_patch->k1; k < V_patch->k2; ++k)
        {
          x1_low.z = x1_high.z = bbox.z1 + (k+0.5) * dz;

          // Get the velocity at the low and high faces, recording 
          // the y-velocities.
          vector_t v_low, v1_low, v_high, v1_high;
          coord_mapping_map_point(adv->mapping, &x1_low, &x_low);
          st_func_eval(adv->V_func, &x_low, t, (real_t*)&v_low);
          coord_mapping_map_vector(adv->inv_mapping, &x_low, &v_low, &v1_low);
          coord_mapping_map_point(adv->mapping, &x1_high, &x_high);
          st_func_eval(adv->V_func, &x_high, t, (real_t*)&v_high);
          coord_mapping_map_vector(adv->inv_mapping, &x_high, &v_high, &v1_high);
          Vy[i][j][k][0] = v1_low.y;
          Vy[i][j+1][k][0] = v1_high.y;

          // Note that since i, j, k are face indices, we need to 
          // translate them to cell indices.
          int ic, jc, kc; // <-- cell indices
          str_grid_patch_translate_indices(V_patch, i, j, k,
                                           U_patch, &ic, &jc, &kc);
#ifdef ADVECT_FIRST_ORDER 
          // Get the values of U at the low and high faces.
          UL[ic][jc][kc][1] = (v_low.y >= 0.0) ? U[ic][jc-1][kc][0] : U[ic][jc][kc][0];
          UH[ic][jc][kc][1] = (v_high.y >= 0.0) ? U[ic][jc][kc][0] : U[ic][jc+1][kc][0];
#else
          // Use the MUSCL-Hancock method to get a 2nd-order accurate flux.
          real_t U_minus = U[ic][jc-1][kc][0];
          real_t U_0     = U[ic][jc][kc][0];
          real_t U_plus  = U[ic][jc+1][kc][0];
          real_t delta_upwind = (v_low.y >= 0.0) ? U_0 - U_minus : U_plus - U_0;
          real_t delta_loc    = (v_low.y >= 0.0) ? U_plus - U_0 : U_0 - U_minus;

          static const real_t omega = 0.0;
          real_t delta = 0.5 * ((1.0 + omega) * delta_upwind + 
                                (1.0 - omega) * delta_loc);

          // Reset small deltas, preserving sign.
          static const real_t tol = 1e-6;
          if (ABS(delta_upwind) < tol)
            delta_upwind = tol * SIGN_VAL(1.0, delta_upwind);
          if (ABS(delta_loc) < tol)
            delta_loc = tol * SIGN_VAL(1.0, delta_loc);

          // Compute the ratio of slopes and apply the appropriate limiter.
          real_t ratio = delta_upwind/delta_loc;
          delta = adv->limit_delta(ratio, omega, delta);

          // Boundary extrapolated values.
          UL[ic][jc][kc][1] = U_0 - 0.5 * delta;
          UH[ic][jc][kc][1] = U_0 + 0.5 * delta;
#endif
        }
      }
    }
  }

  // Traverse the patches and extrapolate U to z-faces.
  pos = 0;
  #pragma omp parallel firstprivate(pos) private(ip, jp, kp, V_patch, bbox)
  while (str_grid_face_data_next_z_patch(V, &pos, &ip, &jp, &kp, &V_patch, &bbox))
  {
    str_grid_patch_t* U_patch = str_grid_cell_data_patch(U, ip, jp, kp);
    str_grid_patch_t* UL_patch = str_grid_cell_data_patch(UL, ip, jp, kp);
    str_grid_patch_t* UH_patch = str_grid_cell_data_patch(UH, ip, jp, kp);

    // Compute "low" and "high" values of U on the x-faces of the cells.
    DECLARE_STR_GRID_PATCH_ARRAY(U, U_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(UL, UL_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(UH, UH_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Vz, V_patch);

    point_t x_low, x1_low, x_high, x1_high;
    for (int i = V_patch->i1; i < V_patch->i2; ++i) 
    {
      x1_low.x = x1_high.x = bbox.x1 + (i+0.5) * dx;
      for (int j = V_patch->j1; j < V_patch->j2; ++j)
      {
        x1_low.y = x1_high.y = bbox.y1 + (j+0.5) * dy;
        for (int k = V_patch->k1; k < V_patch->k2-1; ++k)
        {
          x1_low.z = bbox.z1 + k * dz;
          x1_high.z = x1_low.z + dz;

          // Get the velocity at the low and high faces, recording 
          // the x-velocities.
          vector_t v_low, v1_low, v_high, v1_high;
          coord_mapping_map_point(adv->mapping, &x1_low, &x_low);
          st_func_eval(adv->V_func, &x_low, t, (real_t*)&v_low);
          coord_mapping_map_vector(adv->inv_mapping, &x_low, &v_low, &v1_low);
          coord_mapping_map_point(adv->mapping, &x1_high, &x_high);
          st_func_eval(adv->V_func, &x_high, t, (real_t*)&v_high);
          coord_mapping_map_vector(adv->inv_mapping, &x_high, &v_high, &v1_high);
          Vz[i][j][k][0] = v1_low.z;
          Vz[i][j][k+1][0] = v1_high.z;

          // Note that since i, j, k are face indices, we need to 
          // translate them to cell indices.
          int ic, jc, kc; // <-- cell indices
          str_grid_patch_translate_indices(V_patch, i, j, k,
                                           U_patch, &ic, &jc, &kc);
#ifdef ADVECT_FIRST_ORDER 
          // Get the values of U at the low and high faces.
          UL[ic][jc][kc][2] = (v_low.z >= 0.0) ? U[ic][jc][kc-1][0] : U[ic][jc][kc][0];
          UH[ic][jc][kc][2] = (v_high.z >= 0.0) ? U[ic][jc][kc][0] : U[ic][jc][kc+1][0];
#else
          // Use the MUSCL-Hancock method to get a 2nd-order accurate flux.
          real_t U_minus = U[ic][jc][kc-1][0];
          real_t U_0     = U[ic][jc][kc][0];
          real_t U_plus  = U[ic][jc][kc+1][0];
          real_t delta_upwind = (v_low.z >= 0.0) ? U_0 - U_minus : U_plus - U_0;
          real_t delta_loc    = (v_low.z >= 0.0) ? U_plus - U_0 : U_0 - U_minus;

          static const real_t omega = 0.0;
          real_t delta = 0.5 * ((1.0 + omega) * delta_upwind + 
                                (1.0 - omega) * delta_loc);

          // Reset small deltas, preserving sign.
          static const real_t tol = 1e-6;
          if (ABS(delta_upwind) < tol)
            delta_upwind = tol * SIGN_VAL(1.0, delta_upwind);
          if (ABS(delta_loc) < tol)
            delta_loc = tol * SIGN_VAL(1.0, delta_loc);

          // Compute the ratio of slopes and apply the appropriate limiter.
          real_t ratio = delta_upwind/delta_loc;
          delta = adv->limit_delta(ratio, omega, delta);

          // Boundary extrapolated values.
          UL[ic][jc][kc][2] = U_0 - 0.5 * delta;
          UH[ic][jc][kc][2] = U_0 + 0.5 * delta;
#endif
        }
      }
    }
  }
}
                      
static void compute_fluxes(advect_t* adv,
                           real_t t,
                           str_grid_cell_data_t* UL,
                           str_grid_cell_data_t* UH,
                           str_grid_face_data_t* V,
                           str_grid_face_data_t* F)
{
  real_t dx = adv->dx, dy = adv->dy, dz = adv->dz;
  real_t Ax = dy * dz, Ay = dz * dx, Az = dx * dy;

  int pos, ip, jp, kp;
  str_grid_patch_t* F_patch;
  bbox_t bbox;

  // Traverse the patches and compute the fluxes through x-faces.
  pos = 0;
  #pragma omp parallel firstprivate(pos) private(ip, jp, kp, F_patch, bbox)
  while (str_grid_face_data_next_x_patch(F, &pos, &ip, &jp, &kp, &F_patch, &bbox))
  {
    str_grid_patch_t* UL_patch = str_grid_cell_data_patch(UL, ip, jp, kp);
    str_grid_patch_t* UH_patch = str_grid_cell_data_patch(UH, ip, jp, kp);
    str_grid_patch_t* Vx_patch = str_grid_face_data_x_patch(V, ip, jp, kp);

    DECLARE_STR_GRID_PATCH_ARRAY(UL, UL_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(UH, UH_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Vx, Vx_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Fx, F_patch);

    for (int i = F_patch->i1; i < F_patch->i2; ++i)
    {
      for (int j = F_patch->j1; j < F_patch->j2; ++j)
      {
        for (int k = F_patch->k1; k < F_patch->k2; ++k)
        {
          int ic, jc, kc; // <-- cell indices
          str_grid_patch_translate_indices(F_patch, i, j, k,
                                           UL_patch, &ic, &jc, &kc);
          real_t V = Vx[i][j][k][0];
          real_t U_flux = (V >= 0.0) ? UH[ic-1][jc][kc][0] : UL[ic][jc][kc][0];
          Fx[i][j][k][0] = Ax * V * U_flux;
        }
      }
    }
  }

  // Now compute the fluxes through y-faces.
  pos = 0;
  #pragma omp parallel firstprivate(pos) private(ip, jp, kp, F_patch, bbox)
  while (str_grid_face_data_next_y_patch(F, &pos, &ip, &jp, &kp, &F_patch, &bbox))
  {
    str_grid_patch_t* UL_patch = str_grid_cell_data_patch(UL, ip, jp, kp);
    str_grid_patch_t* UH_patch = str_grid_cell_data_patch(UH, ip, jp, kp);
    str_grid_patch_t* Vy_patch = str_grid_face_data_y_patch(V, ip, jp, kp);

    DECLARE_STR_GRID_PATCH_ARRAY(UL, UL_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(UH, UH_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Vy, Vy_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Fy, F_patch);

    for (int i = F_patch->i1; i < F_patch->i2; ++i)
    {
      for (int j = F_patch->j1; j < F_patch->j2; ++j)
      {
        for (int k = F_patch->k1; k < F_patch->k2; ++k)
        {
          int ic, jc, kc; // <-- cell indices
          str_grid_patch_translate_indices(F_patch, i, j, k,
                                           UL_patch, &ic, &jc, &kc);
          real_t V = Vy[i][j][k][0];
          real_t U_flux = (V >= 0.0) ? UH[ic][jc-1][kc][1] : UL[ic][jc][kc][1];
          Fy[i][j][k][0] = Ay * V * U_flux;
        }
      }
    }
  }

  // Now through z-faces.
  pos = 0;
  #pragma omp parallel firstprivate(pos) private(ip, jp, kp, F_patch, bbox)
  while (str_grid_face_data_next_z_patch(F, &pos, &ip, &jp, &kp, &F_patch, &bbox))
  {
    str_grid_patch_t* UL_patch = str_grid_cell_data_patch(UL, ip, jp, kp);
    str_grid_patch_t* UH_patch = str_grid_cell_data_patch(UH, ip, jp, kp);
    str_grid_patch_t* Vz_patch = str_grid_face_data_z_patch(V, ip, jp, kp);

    DECLARE_STR_GRID_PATCH_ARRAY(UL, UL_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(UH, UH_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Vz, Vz_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Fz, UL_patch);
    for (int i = F_patch->i1; i < F_patch->i2; ++i)
    {
      for (int j = F_patch->j1; j < F_patch->j2; ++j)
      {
        for (int k = F_patch->k1; k < F_patch->k2; ++k)
        {
          int ic, jc, kc; // <-- cell indices
          str_grid_patch_translate_indices(F_patch, i, j, k,
                                           UL_patch, &ic, &jc, &kc);
          real_t V = Vz[i][j][k][0];
          real_t U_flux = (V >= 0.0) ? UH[ic][jc][kc-1][2] : UL[ic][jc][kc][2];
          Fz[i][j][k][0] = Az * V * U_flux;
        }
      }
    }
  }
}

// Here's the right-hand side for the ARK integrator.
static int advect_rhs(void* context, 
                      real_t t, 
                      str_grid_cell_data_t* U, 
                      str_grid_cell_data_t* dUdt)
{
  advect_t* adv = context;
  real_t dx = adv->dx, dy = adv->dy, dz = adv->dz;

  // Make sure ghost cells are filled.
  str_grid_cell_filler_fill(adv->ghost_filler, U);

  // Extrapolate U to faces and compute face velocities.
  extrapolate_U_to_faces(adv, t, U, adv->UL, adv->UH, adv->V);

  // Fill ghost cells for UL and UH, the extrapolated values of U.
  str_grid_cell_filler_fill(adv->ghost_filler, adv->UL);
  str_grid_cell_filler_fill(adv->ghost_filler, adv->UH);

  // Compute fluxes.
  compute_fluxes(adv, t, adv->UL, adv->UH, adv->V, adv->F);

  // Finally, compute the flux divergence.
  real_t V = dx * dy * dz;
  int pos = 0, ip, jp, kp;
  str_grid_patch_t* rhs_patch;
  #pragma omp parallel firstprivate(pos) private(ip, jp, kp, rhs_patch)
  while (str_grid_cell_data_next_patch(dUdt, &pos, &ip, &jp, &kp, &rhs_patch, NULL))
  {
    str_grid_patch_t* Fx_patch = str_grid_face_data_x_patch(adv->F, ip, jp, kp);
    str_grid_patch_t* Fy_patch = str_grid_face_data_y_patch(adv->F, ip, jp, kp);
    str_grid_patch_t* Fz_patch = str_grid_face_data_z_patch(adv->F, ip, jp, kp);
    DECLARE_STR_GRID_PATCH_ARRAY(dUdt, rhs_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Fx, Fx_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Fy, Fy_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(Fz, Fz_patch);

    for (int i = rhs_patch->i1; i < rhs_patch->i2; ++i)
    {
      for (int j = rhs_patch->j1; j < rhs_patch->j2; ++j)
      {
        for (int k = rhs_patch->k1; k < rhs_patch->k2; ++k)
        {
          // Translate i, j, k (cell indices) to the various x-, y-, and z-face indices.
          int ifx, jfx, kfx;
          str_grid_patch_translate_indices(rhs_patch, i, j, k,
                                           Fx_patch, &ifx, &jfx, &kfx);
          int ify, jfy, kfy;
          str_grid_patch_translate_indices(rhs_patch, i, j, k,
                                           Fy_patch, &ify, &jfy, &kfy);
          int ifz, jfz, kfz;
          str_grid_patch_translate_indices(rhs_patch, i, j, k,
                                           Fz_patch, &ifz, &jfz, &kfz);

          // Compute the right hand side.
          real_t div_F = Fx[ifx+1][jfx][kfx][0] - Fx[ifx][jfx][kfx][0] + 
                         Fy[ify][jfy+1][kfy][0] - Fy[ify][jfy][kfy][0] + 
                         Fz[ifz][jfz][kfz+1][0] - Fz[ifz][jfz][kfz][0];
          dUdt[i][j][k][0] = -div_F / V;
        }
      }
    }
  }

  return 0;
}

static void override_algorithm_options(options_t* options,
                                       real_t (**limiter_func)(real_t, real_t, real_t))
{
  *limiter_func = van_leer_limiter;

  // Override flux limiter?
  char* limiter_str = options_value(options, "limiter");
  if (limiter_str != NULL) 
  {
    const char* list[] = {"none", "superbee", "van_leer", "van_albada", 
                          "minmod", NULL};
    int l = string_find_in_list(limiter_str, list, false);
    if (l == -1)
      polymec_error("Invalid limiter: %s (must be one of none, superbee, van_leer, van_albada, minmod).");
    switch(l) 
    {
      case 0: *limiter_func = no_limiter; break;
      case 1: *limiter_func = superbee_limiter; break;
      case 2: *limiter_func = van_leer_limiter; break;
      case 3: *limiter_func = van_albada_limiter; break;
      case 4: *limiter_func = minmod_limiter;
    }
  }
}

static void advect_setup(advect_t* adv)
{
  ASSERT(st_func_num_comp(adv->U0) == 1);
  ASSERT(st_func_num_comp(adv->V_func) == 3);

  advect_clear(adv);

  // Adjust our algorithm options if needed.
  options_t* options = options_argv();
  override_algorithm_options(options, &adv->limit_delta);

  // We impose a zero flux on all boundaries, on the assumption that the 
  // solution is zero at the boundary.
  str_grid_patch_filler_t* zero_flux_x1 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_X1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_x2 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_X2_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_y1 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Y1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_y2 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Y2_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_z1 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Z1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_z2 = zero_flux_str_grid_patch_filler_new(STR_GRID_PATCH_Z2_BOUNDARY);

  str_grid_cell_filler_factory_t* factory = str_grid_cell_filler_factory_new(MPI_COMM_WORLD);
  adv->ghost_filler = str_grid_cell_filler_factory_ghost_filler(factory, adv->grid, 
                                                                zero_flux_x1, zero_flux_x2,
                                                                zero_flux_y1, zero_flux_y2,
                                                                zero_flux_z1, zero_flux_z2);
  factory = NULL;

  // Set the grid spacings.
  int npx, npy, npz;
  str_grid_get_extents(adv->grid, &npx, &npy, &npz);
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
  adv->F = str_grid_face_data_new(adv->grid, 1);
  adv->UL = str_grid_cell_data_new(adv->grid, 3*num_comps, num_ghost_layers);
  adv->UH = str_grid_cell_data_new(adv->grid, 3*num_comps, num_ghost_layers);
  adv->V = str_grid_face_data_new(adv->grid, 1);

  // Create the ARK integrator.
  adv->integ = explicit_ark_str_ode_integrator_new(2,
                                                   MPI_COMM_WORLD,
                                                   adv->grid,
                                                   num_comps,
                                                   num_ghost_layers,
                                                   adv,
                                                   advect_rhs,
                                                   ark_max_dt,
                                                   NULL);

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
  #pragma omp parallel firstprivate(pos) private(ip, jp, kp, patch, bbox)
  while (str_grid_cell_data_next_patch(adv->U, &pos, &ip, &jp, &kp, &patch, &bbox))
  {
    DECLARE_STR_GRID_PATCH_ARRAY(U, patch);
    point_t x, x1;
    for (int i = patch->i1; i < patch->i2; ++i)
    {
      x1.x = bbox.x1 + (i + 0.5) * adv->dx;
      for (int j = patch->j1; j < patch->j2; ++j)
      {
        x1.y = bbox.y1 + (j + 0.5) * adv->dy;
        for (int k = patch->k1; k < patch->k2; ++k)
        {
          x1.z = bbox.z1 + (k + 0.5) * adv->dz;
          coord_mapping_map_point(adv->mapping, &x1, &x);
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
  polymec_suspend_fpe();
  if (!str_ode_integrator_step(adv->integ, max_dt, &t2, adv->U))
    polymec_error("advect_advance: Integration failed at t = %g.", t);
  polymec_restore_fpe();

  return t2 - t;
}

static void advect_plot(void* context, const char* prefix, const char* directory, real_t t, int step)
{
  advect_t* adv = context;

  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, prefix, directory, 1, 0, step, t);
  silo_file_write_str_grid(silo, "grid", adv->grid, adv->mapping);
  const char* U_name[] = {"U"};
  silo_file_write_str_grid_cell_data(silo, U_name, "grid", adv->U, NULL, adv->mapping);

  // Compute and plot the velocity field.
  {
    const char* V_names[] = {"vx", "vy", "vz"};
    str_grid_cell_data_t* velocity = str_grid_cell_data_new(adv->grid, 3, 0);
    int pos = 0, ip, jp, kp;
    str_grid_patch_t* patch;
    bbox_t bbox;
    #pragma omp parallel firstprivate(pos) private(ip, jp, kp, patch, bbox)
    while (str_grid_cell_data_next_patch(velocity, &pos, &ip, &jp, &kp, &patch, &bbox))
    {
      DECLARE_STR_GRID_PATCH_ARRAY(V, patch);
      point_t x1, x;
      for (int i = patch->i1; i < patch->i2; ++i)
      {
        x1.x = bbox.x1 + (i + 0.5) * adv->dx;
        for (int j = patch->j1; j < patch->j2; ++j)
        {
          x1.y = bbox.y1 + (j + 0.5) * adv->dy;
          for (int k = patch->k1; k < patch->k2; ++k)
          {
            x1.z = bbox.z1 + (k + 0.5) * adv->dz;
            coord_mapping_map_point(adv->mapping, &x1, &x);
            st_func_eval(adv->V_func, &x, t, &V[i][j][k][0]);
          }
        }
      }
    }
    silo_file_write_str_grid_cell_data(silo, V_names, "grid", velocity, NULL, adv->mapping);
    silo_file_write_vector_expression(silo, "velocity", "{<vx>, <vy>, <vz>}");
    str_grid_cell_data_free(velocity);
  }

  // Compute and plot the time derivative of the solution.
  {
    const char* dUdt_name[] = {"dUdt"};
    str_grid_cell_data_t* dUdt = str_grid_cell_data_new(adv->grid, 1, 1);
    advect_rhs(context, t, adv->U, dUdt);
    silo_file_write_str_grid_cell_data(silo, dUdt_name, "grid", dUdt, NULL, adv->mapping);
    str_grid_cell_data_free(dUdt);
  }

  // Wrap it up.
  silo_file_close(silo);
}

static void advect_load(void* context, const char* prefix, const char* directory, real_t* t, int step)
{
  advect_t* adv = context;
  POLYMEC_NOT_IMPLEMENTED;
}

static void advect_save(void* context, const char* prefix, const char* directory, real_t t, int step)
{
  advect_t* adv = context;

  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, prefix, directory, 1, 0, step, t);
  silo_file_write_str_grid(silo, "grid", adv->grid, adv->mapping);
  const char* comp_names[] = {"U"};
  silo_file_write_str_grid_cell_data(silo, comp_names, "grid", adv->U, NULL, adv->mapping);
  silo_file_close(silo);
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
  adv->mapping = NULL;
  adv->inv_mapping = NULL;
  adv->U0 = NULL;
  adv->V_func = NULL;
  adv->V = NULL;
  adv->U = NULL;
  adv->UL = NULL;
  adv->UH = NULL;
  adv->F = NULL;
  adv->ghost_filler = NULL;
  adv->integ = NULL;
  model_vtable vtable = {.read_input = advect_read_input,
                         .init = advect_init,
                         .max_dt = advect_max_dt,
                         .advance = advect_advance,
                         .plot = advect_plot,
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
