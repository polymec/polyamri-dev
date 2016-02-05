// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/polyamri.h"
#include "polyamri/str_grid.h"
#include "polyamri/str_grid_cell_data.h"
#include "polyamri/str_grid_face_data.h"
#include "polyamri/str_ode_integrator.h"
#include "model/model.h"

typedef struct
{
  MPI_Comm comm;

  // Problem metadata.
  int nx, ny, nz;     // cells
  int npx, npy, npz;  // patches
  bool x_periodic, y_periodic, z_periodic;
  bbox_t bbox;
  real_t dx, dy, dz;

  // State information.
  str_grid_t* grid;
  str_grid_cell_data_t* U;

  // Initial state.
  st_func_t* U0;

  // Velocity field.
  st_func_t* V;

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
                              interpreter_t* interpreter, 
                              options_t* options)
{
}

static void advect_init(void* context, real_t t)
{
  advect_t* adv = context;
  ASSERT(adv->nx > 0);
  ASSERT(adv->ny > 0);
  ASSERT(adv->nz > 0);
  ASSERT(adv->npx > 0);
  ASSERT(adv->npy > 0);
  ASSERT(adv->npz > 0);
  ASSERT(!bbox_is_empty_set(&adv->bbox));
  ASSERT(!bbox_is_point(&adv->bbox));
  ASSERT(!bbox_is_line(&adv->bbox));
  ASSERT(!bbox_is_plane(&adv->bbox));
  ASSERT(st_func_num_comp(adv->U0) == 1);
  ASSERT(st_func_num_comp(adv->V) == 3);

  // Set the grid spacings.
  int Nx = adv->npx * adv->nx;
  int Ny = adv->npy * adv->ny;
  int Nz = adv->npz * adv->nz;
  adv->dx = (adv->bbox.x2 - adv->bbox.x1) / Nx;
  adv->dy = (adv->bbox.y2 - adv->bbox.y1) / Ny;
  adv->dz = (adv->bbox.z2 - adv->bbox.z1) / Nz;

  // Create the grid and the state vector. 
  log_detail("str_advect: Creating %d x %d x %d grid...", Nx, Ny, Nz);
  adv->grid = str_grid_new(adv->npx, adv->npy, adv->npz, 
                           adv->nx, adv->ny, adv->nz, 
                           adv->x_periodic, adv->y_periodic, adv->z_periodic);
  adv->U = str_grid_cell_data_new(adv->grid, 1, 1);

  // Create work "vectors."
  adv->U_work = str_grid_cell_data_with_buffer(adv->grid, 1, 1, NULL);
  adv->dUdt_work = str_grid_cell_data_with_buffer(adv->grid, 1, 1, NULL);
  adv->F_work = str_grid_face_data_new(adv->grid, 1);

  // Initialize the state.
  int pos = 0, ip, jp, kp;
  str_grid_patch_t* patch;
  while (str_grid_cell_data_next_patch(adv->U, &pos, &ip, &jp, &kp, &patch))
  {
    DECLARE_STR_GRID_PATCH_ARRAY(U, patch);
    point_t x;
    for (int i = patch->i1; i < patch->i2; ++i)
    {
      int I = i - patch->i1;
      x.x = (ip * adv->nx + (I + 0.5)) * adv->dx;
      for (int j = patch->j1; j < patch->j2; ++j)
      {
        int J = j - patch->j1;
        x.y = (jp * adv->ny + (J + 0.5)) * adv->dy;
        for (int k = patch->k1; k < patch->k2; ++k)
        {
          int K = j - patch->k1;
          x.z = (kp * adv->nz + (K + 0.5)) * adv->dz;
          st_func_eval(adv->U0, &x, t, &U[i][j][k][0]);
        }
      }
    }
  }
}

static real_t advect_max_dt(void* context, real_t t, char* reason)
{
  advect_t* adv = context;
  real_t dx = MIN(adv->dx, MIN(adv->dy, adv->dz));

  // We compute the maximum allowable timestep using the CFL condition, 
  // measuring the velocity at each of the face centers at the given time and 
  // taking the max.
  real_t V_max = 0.0;
  int pos = 0, ip, jp, kp;
  str_grid_patch_t* patch;
  while (str_grid_cell_data_next_patch(adv->U, &pos, &ip, &jp, &kp, &patch))
  {
    point_t x;
    vector_t v;

    // x face velocities.
    for (int i = 0; i <= adv->nx; ++i)
    {
      x.x = (ip * adv->nx + i) * adv->dx;
      for (int j = 0; j < adv->ny; ++j)
      {
        x.y = (jp * adv->ny + (j + 0.5)) * adv->dy;
        for (int k = 0; k < adv->nz; ++k)
        {
          x.z = (kp * adv->nz + (k + 0.5)) * adv->dz;
          st_func_eval(adv->V, &x, t, (real_t*)&v);
          V_max = MAX(V_max, vector_mag(&v));
        }
      }
    }

    // y face velocities.
    for (int i = 0; i < adv->nx; ++i)
    {
      x.x = (ip * adv->nx + (i + 0.5)) * adv->dx;
      for (int j = 0; j <= adv->ny; ++j)
      {
        x.y = (jp * adv->ny + j) * adv->dy;
        for (int k = 0; k < adv->nz; ++k)
        {
          x.z = (kp * adv->nz + (k + 0.5)) * adv->dz;
          st_func_eval(adv->V, &x, t, (real_t*)&v);
          V_max = MAX(V_max, vector_mag(&v));
        }
      }
    }

    // z face velocities.
    for (int i = 0; i < adv->nx; ++i)
    {
      x.x = (ip * adv->nx + (i + 0.5)) * adv->dx;
      for (int j = 0; j < adv->ny; ++j)
      {
        x.y = (jp * adv->ny + (j + 0.5)) * adv->dy;
        for (int k = 0; k <= adv->nz; ++k)
        {
          x.z = (kp * adv->nz + k) * adv->dz;
          st_func_eval(adv->V, &x, t, (real_t*)&v);
          V_max = MAX(V_max, vector_mag(&v));
        }
      }
    }
  }

  if (V_max > 0.0)
    return dx / V_max;
  else 
    return FLT_MAX;
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
  
  // x fluxes.
  pos = 0;
  real_t dx = adv->dx, dy = adv->dy, dz = adv->dz;
  while (str_grid_face_data_next_x_patch(adv->F_work, &pos, &ip, &jp, &kp, &face_patch))
  {
    str_grid_patch_t* cell_patch = str_grid_cell_data_patch(adv->U_work, ip, jp, kp);
    DECLARE_STR_GRID_PATCH_ARRAY(F, face_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(U, cell_patch);
    point_t x_low, x_high;
    for (int i = face_patch->i1; i < face_patch->i2; ++i)
    {
      int I = i - face_patch->i1;
      x_low.x = (ip * adv->nx + I) * dx;
      x_high.x = x_low.x + dx;
      for (int j = face_patch->j1; j < face_patch->j2; ++j)
      {
        int J = j - face_patch->j1;
        x_low.y = x_high.y = (jp * adv->ny + J) * dy;
        for (int k = face_patch->k1; k < face_patch->k2; ++k)
        {
          int K = k - face_patch->k1;
          x_low.z = x_high.z = (kp * adv->nz + K) * dz;

          // Get the velocity at the low and high faces.
          vector_t v_low, v_high;
          st_func_eval(adv->V, &x_low, t, (real_t*)&v_low);
          st_func_eval(adv->V, &x_high, t, (real_t*)&v_high);

          // Get the upwinded values of U at the low and high faces.
          real_t U_low = (v_low.x >= 0.0) ? U[i-1][j][k][0] : U[i][j][k][0];
          real_t U_high = (v_high.x >= 0.0) ? U[i][j][k][0] : U[i+1][j][k][0];

          // Compute the low and high fluxes.
          F[i][j][k][0]   = U_low * v_low.x;
          F[i+1][j][k][0] = U_high * v_high.x;
        }
      }
    }
  }

  // y fluxes.
  pos = 0;
  while (str_grid_face_data_next_y_patch(adv->F_work, &pos, &ip, &jp, &kp, &face_patch))
  {
    str_grid_patch_t* cell_patch = str_grid_cell_data_patch(adv->U_work, ip, jp, kp);
    DECLARE_STR_GRID_PATCH_ARRAY(F, face_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(U, cell_patch);
    point_t x_low, x_high;
    for (int i = face_patch->i1; i < face_patch->i2; ++i)
    {
      int I = i - face_patch->i1;
      x_low.x = x_high.x = (ip * adv->nx + I) * dx;
      for (int j = face_patch->j1; j < face_patch->j2; ++j)
      {
        int J = j - face_patch->j1;
        x_low.y = (jp * adv->ny + (J+0.5)) * dy;
        x_high.y = x_low.y + dy;
        for (int k = face_patch->k1; k < face_patch->k2; ++k)
        {
          int K = k - face_patch->k1;
          x_low.z = x_high.z = (kp * adv->nz + (K+0.5)) * dz;

          // Get the velocity at the low and high faces.
          vector_t v_low, v_high;
          st_func_eval(adv->V, &x_low, t, (real_t*)&v_low);
          st_func_eval(adv->V, &x_high, t, (real_t*)&v_high);

          // Get the upwinded values of U at the low and high faces.
          real_t U_low = (v_low.y >= 0.0) ? U[i][j-1][k][0] : U[i][j][k][0];
          real_t U_high = (v_high.y >= 0.0) ? U[i][j][k][0] : U[i][j+1][k][0];

          // Compute the low and high fluxes.
          F[i][j][k][0]   = U_low * v_low.y;
          F[i+1][j][k][0] = U_high * v_high.y;
        }
      }
    }
  }

  // z fluxes.
  pos = 0;
  while (str_grid_face_data_next_z_patch(adv->F_work, &pos, &ip, &jp, &kp, &face_patch))
  {
    str_grid_patch_t* cell_patch = str_grid_cell_data_patch(adv->U_work, ip, jp, kp);
    DECLARE_STR_GRID_PATCH_ARRAY(F, face_patch);
    DECLARE_STR_GRID_PATCH_ARRAY(U, cell_patch);
    point_t x_low, x_high;
    for (int i = face_patch->i1; i < face_patch->i2; ++i)
    {
      int I = i - face_patch->i1;
      x_low.x = x_high.x = (ip * adv->nx + (I+0.5)) * dx;
      for (int j = face_patch->j1; j < face_patch->j2; ++j)
      {
        int J = j - face_patch->j1;
        x_low.y = x_high.y = (jp * adv->ny + (J+0.5)) * dy;
        for (int k = face_patch->k1; k < face_patch->k2; ++k)
        {
          int K = k - face_patch->k1;
          x_low.z = (kp * adv->nz + K) * dz;
          x_high.z = x_low.z + dz;

          // Get the velocity at the low and high faces.
          vector_t v_low, v_high;
          st_func_eval(adv->V, &x_low, t, (real_t*)&v_low);
          st_func_eval(adv->V, &x_high, t, (real_t*)&v_high);

          // Get the upwinded values of U at the low and high faces.
          real_t U_low = (v_low.z >= 0.0) ? U[i][j][k-1][0] : U[i][j][k][0];
          real_t U_high = (v_high.z >= 0.0) ? U[i][j][k][0] : U[i][j][k+1][0];

          // Compute the low and high fluxes.
          F[i][j][k][0]   = U_low * v_low.z;
          F[i+1][j][k][0] = U_high * v_high.z;
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
  while (str_grid_cell_data_next_patch(adv->dUdt_work, &pos, &ip, &jp, &kp, &rhs_patch))
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
  silo_file_t* silo = silo_file_open(adv->grid->comm, file_prefix, directory, 0, step, t);
  *t = 1.0 * step;

  ASSERT(adv->grid == NULL);
  ASSERT(adv->state == NULL);
  adv->grid = silo_file_read_mesh(silo, "grid");
  adv->stencil = cell_halo_stencil_new(adv->grid);
  adv->state = silo_file_read_scalar_cell_field(silo, "state", "grid", NULL);
  silo_file_close(silo);
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
  str_grid_free(adv->grid);
  str_grid_cell_data_free(adv->U);
  str_grid_face_data_free(adv->F_work);
  str_grid_cell_data_free(adv->U_work);
  str_grid_cell_data_free(adv->dUdt_work);
  polymec_free(adv);
}

static model_t* advect_ctor()
{
  advect_t* adv = polymec_malloc(sizeof(advect_t));
  adv->grid = NULL;
  adv->U = NULL;
  model_vtable vtable = {.read_input = advect_read_input,
                         .init = advect_init,
                         .max_dt = advect_max_dt,
                         .advance = advect_advance,
                         .load = advect_load,
                         .save = advect_save,
                         .finalize = advect_finalize,
                         .dtor = advect_dtor};
  docstring_t* advect_doc = docstring_from_string(advect_desc);
  return model_new("str_advect", adv, vtable, advect_doc, MODEL_MPI);
}

// Main program.
int main(int argc, char* argv[])
{
  return model_main("str_advect", advect_ctor, argc, argv);
}
