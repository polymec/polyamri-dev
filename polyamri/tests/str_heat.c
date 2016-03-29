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
#include "polyamri/str_grid_solver.h"
#include "polyamri/ark_str_ode_integrator.h"
#include "polyamri/interpreter_register_polyamri_functions.h"
#include "polyamri/grid_to_bbox_coord_mapping.h"
#include "polyamri/silo_file.h"
#include "model/model.h"

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
  real_t current_time;

  // Initial state.
  st_func_t* U0;

  // Thermal diffusivity function and heat source.
  st_func_t* D_func;
  st_func_t* source_func;

  // Linear solver.
  str_grid_solver_t* solver;

  // Time integrator and workspace.
  str_ode_integrator_t* integ;
  char max_dt_reason[POLYMEC_MODEL_MAXDT_REASON_SIZE];
} heat_t;

static const char heat_desc[] = "PolyAMRI heat equation (str_heat)\n"
  "This model demonstrates a parallel implementation of a finite volume\n"
  "discretization of the heat equation on a uniform grid,\n"
  "integrated implicitly in time using the method of lines.\n";

static void heat_read_input(void* context, 
                            interpreter_t* interp, 
                            options_t* options)
{
  heat_t* heat = context;

  heat->U0 = interpreter_get_scalar_function(interp, "initial_state");
  heat->D_func = interpreter_get_vector_function(interp, "diffusivity");
  heat->source_func = interpreter_get_vector_function(interp, "source");
  heat->grid = interpreter_get_str_grid(interp, "grid");

  // Peel the computational domain off the grid.
  heat->mapping = str_grid_property(heat->grid, "mapping");
  if (heat->mapping == NULL)
  {
    // We use an "unmapped" domain.
    bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
    heat->mapping = grid_to_bbox_coord_mapping_new(&bbox);
  }
  
  // Make sure the mapping is invertible.
  heat->inv_mapping = coord_mapping_inverse(heat->mapping);
  if (heat->inv_mapping == NULL)
    polymec_error("domain mapping is not invertible! Must be an invertible mapping.");
}

static void heat_clear(heat_t* heat)
{
  if (heat->solver != NULL)
    str_grid_solver_free(heat->solver);
  if (heat->integ != NULL)
    str_ode_integrator_free(heat->integ);
  if (heat->U != NULL)
    str_grid_cell_data_free(heat->U);
  if (heat->ghost_filler != NULL)
    str_grid_cell_filler_free(heat->ghost_filler);
}

// Here's the right-hand side for the ARK integrator.
static int heat_rhs(void* context, 
                    real_t t, 
                    str_grid_cell_data_t* U, 
                    str_grid_cell_data_t* dUdt)
{
  heat_t* heat = context;
  real_t dx = heat->dx, dy = heat->dy, dz = heat->dz;

  // Make sure ghost cells are filled.
  str_grid_cell_filler_fill(heat->ghost_filler, U);

  // Solve for the next solution, storing it in dU/dt for now.
  hypre_smg_helmholtz_str_grid_cell_solver_set_operator(heat->solver, t, 1.0, dt, heat->D_func);
  hypre_smg_helmholtz_str_grid_cell_solver_set_rhs(heat->solver, t, 0.0, 0.0, NULL, 1.0, heat->source_func);
  bool solved = str_grid_solver_solve(heat->solver, dUdt);
  if (!solved) return -1;

  // Compute the real dU/dt using a simple finite difference.
  real_t dt = t - heat->current_time;
  int pos = 0, ip, jp, kp;
  str_grid_patch_t* patch;
  bbox_t bbox;
  #pragma omp parallel firstprivate(pos) private(ip, jp, kp, patch, bbox)
  while (str_grid_cell_data_next_patch(heat->U, &pos, &ip, &jp, &kp, &patch, &bbox))
  {
    DECLARE_STR_GRID_PATCH_ARRAY(U, patch);
    DECLARE_STR_GRID_PATCH_ARRAY(dUdt, patch);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int j = patch->j1; j < patch->j2; ++j)
        for (int k = patch->k1; k < patch->k2; ++k)
          dUdt[i][j][k][0] = (dUdt[i][j][k][0] - U[i][j][k][0]) / dt;
  }

  return 0;
}

static void override_algorithm_options(options_t* options)
{
}

static void heat_setup(heat_t* heat, real_t t)
{
  ASSERT(st_func_num_comp(heat->U0) == 1);
  ASSERT(st_func_num_comp(heat->D_func) == 1);
  ASSERT(st_func_num_comp(heat->source_func) == 1);

  heat_clear(heat);
  heat->current_time = t;

  // Adjust our algorithm options if needed.
  options_t* options = options_argv();
  override_algorithm_options(options);

  // We impose a zero flux on all boundaries, on the assumption that the 
  // solution is zero at the boundary.
  str_grid_patch_filler_t* zero_flux_x1 = neumann_bc_str_grid_patch_filler_new(1.0, 0.0, 1.0, 0, STR_GRID_PATCH_X1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_x2 = neumann_bc_str_grid_patch_filler_new(1.0, 0.0, 1.0, 0, STR_GRID_PATCH_X2_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_y1 = neumann_bc_str_grid_patch_filler_new(1.0, 0.0, 1.0, 0, STR_GRID_PATCH_Y1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_y2 = neumann_bc_str_grid_patch_filler_new(1.0, 0.0, 1.0, 0, STR_GRID_PATCH_Y2_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_z1 = neumann_bc_str_grid_patch_filler_new(1.0, 0.0, 1.0, 0, STR_GRID_PATCH_Z1_BOUNDARY);
  str_grid_patch_filler_t* zero_flux_z2 = neumann_bc_str_grid_patch_filler_new(1.0, 0.0, 1.0, 0, STR_GRID_PATCH_Z2_BOUNDARY);

  str_grid_cell_filler_factory_t* factory = str_grid_cell_filler_factory_new(MPI_COMM_WORLD);
  heat->ghost_filler = str_grid_cell_filler_factory_ghost_filler(factory, heat->grid, 
                                                                zero_flux_x1, zero_flux_x2,
                                                                zero_flux_y1, zero_flux_y2,
                                                                zero_flux_z1, zero_flux_z2);
  factory = NULL;

  // Set the grid spacings.
  int npx, npy, npz;
  str_grid_get_extents(heat->grid, &npx, &npy, &npz);
  int nx, ny, nz;
  str_grid_get_patch_size(heat->grid, &nx, &ny, &nz);
  int Nx = npx * nx;
  int Ny = npy * ny;
  int Nz = npz * nz;
  heat->dx = 1.0 / Nx;
  heat->dy = 1.0 / Ny;
  heat->dz = 1.0 / Nz;

  // Create the state vector.
  int num_comps = 1;
  int num_ghost_layers = 1;
  heat->U = str_grid_cell_data_new(heat->grid, num_comps, num_ghost_layers);

  // Create the linear solver.
  heat->solver = hypre_smg_helmholtz_str_grid_cell_solver_new(heat->grid, num_comps);

  // Create the implicit functional ARK integrator.
  heat->integ = functional_ark_str_ode_integrator_new(2,
                                                      MPI_COMM_WORLD,
                                                      heat->grid,
                                                      num_comps,
                                                      num_ghost_layers,
                                                      heat,
                                                      NULL,
                                                      heat_rhs,
                                                      NULL,
                                                      NULL,
                                                      2);
}

static void heat_init(void* context, real_t t)
{
  heat_t* heat = context;
  heat_setup(heat, t);

  // Initialize the state.
  int pos = 0, ip, jp, kp;
  str_grid_patch_t* patch;
  bbox_t bbox;
  #pragma omp parallel firstprivate(pos) private(ip, jp, kp, patch, bbox)
  while (str_grid_cell_data_next_patch(heat->U, &pos, &ip, &jp, &kp, &patch, &bbox))
  {
    DECLARE_STR_GRID_PATCH_ARRAY(U, patch);
    point_t x, x1;
    for (int i = patch->i1; i < patch->i2; ++i)
    {
      x1.x = bbox.x1 + (i + 0.5) * heat->dx;
      for (int j = patch->j1; j < patch->j2; ++j)
      {
        x1.y = bbox.y1 + (j + 0.5) * heat->dy;
        for (int k = patch->k1; k < patch->k2; ++k)
        {
          x1.z = bbox.z1 + (k + 0.5) * heat->dz;
          coord_mapping_map_point(heat->mapping, &x1, &x);
          st_func_eval(heat->U0, &x, t, &U[i][j][k][0]);
        }
      }
    }
  }
}

static real_t heat_advance(void* context, real_t max_dt, real_t t)
{
  heat_t* heat = context;

  // Integrate!
  real_t t2 = t;
  polymec_suspend_fpe();
  if (!str_ode_integrator_step(heat->integ, max_dt, &t2, heat->U))
    polymec_error("heat_advance: Integration failed at t = %g.", t);
  polymec_restore_fpe();

  real_t dt = t2 - t;
  heat->current_time += dt;
  return dt;
}

static void heat_plot(void* context, const char* prefix, const char* directory, real_t t, int step)
{
  heat_t* heat = context;

  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, prefix, directory, 1, 0, step, t);
  silo_file_write_str_grid(silo, "grid", heat->grid, heat->mapping);
  const char* U_name[] = {"U"};
  silo_file_write_str_grid_cell_data(silo, U_name, "grid", heat->U, NULL, heat->mapping);

  // Compute and plot the diffusivity.
  {
    const char* D_name[] = {"diffusivity"};
    str_grid_cell_data_t* D = str_grid_cell_data_new(heat->grid, 1, 0);
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
        x1.x = bbox.x1 + (i + 0.5) * heat->dx;
        for (int j = patch->j1; j < patch->j2; ++j)
        {
          x1.y = bbox.y1 + (j + 0.5) * heat->dy;
          for (int k = patch->k1; k < patch->k2; ++k)
          {
            x1.z = bbox.z1 + (k + 0.5) * heat->dz;
            coord_mapping_map_point(heat->mapping, &x1, &x);
            st_func_eval(heat->D_func, &x, t, &D[i][j][k][0]);
          }
        }
      }
    }
    silo_file_write_str_grid_cell_data(silo, U_name, "grid", D, NULL, heat->mapping);
    str_grid_cell_data_free(D);
  }

  // Compute and plot the time derivative of the solution.
  {
    const char* dUdt_name[] = {"dUdt"};
    str_grid_cell_data_t* dUdt = str_grid_cell_data_new(heat->grid, 1, 1);
    heat_rhs(context, t, heat->U, dUdt);
    silo_file_write_str_grid_cell_data(silo, dUdt_name, "grid", dUdt, NULL, heat->mapping);
    str_grid_cell_data_free(dUdt);
  }

  // Wrap it up.
  silo_file_close(silo);
}

static void heat_load(void* context, const char* prefix, const char* directory, real_t* t, int step)
{
  heat_t* heat = context;
  POLYMEC_NOT_IMPLEMENTED;
}

static void heat_save(void* context, const char* prefix, const char* directory, real_t t, int step)
{
  heat_t* heat = context;

  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, prefix, directory, 1, 0, step, t);
  silo_file_write_str_grid(silo, "grid", heat->grid, heat->mapping);
  const char* comp_names[] = {"U"};
  silo_file_write_str_grid_cell_data(silo, comp_names, "grid", heat->U, NULL, heat->mapping);
  silo_file_close(silo);
}

static void heat_finalize(void* context, int step, real_t t)
{
//  heat_t* heat = context;
}

static void heat_dtor(void* context)
{
  heat_t* heat = context;

  heat_clear(heat);
  str_grid_free(heat->grid);
  polymec_free(heat);
}

static model_t* heat_ctor()
{
  heat_t* heat = polymec_malloc(sizeof(heat_t));
  heat->grid = NULL;
  heat->mapping = NULL;
  heat->inv_mapping = NULL;
  heat->U0 = NULL;
  heat->D_func = NULL;
  heat->source_func = NULL;
  heat->U = NULL;
  heat->ghost_filler = NULL;
  heat->integ = NULL;
  model_vtable vtable = {.read_input = heat_read_input,
                         .init = heat_init,
                         .max_dt = heat_max_dt,
                         .advance = heat_advance,
                         .plot = heat_plot,
                         .load = heat_load,
                         .save = heat_save,
                         .finalize = heat_finalize,
                         .dtor = heat_dtor};
  docstring_t* heat_doc = docstring_from_string(heat_desc);
  model_t* model = model_new("str_heat", heat, vtable, heat_doc, MODEL_MPI);

  // Register polyamri-specific functions.
  interpreter_t* interp = model_interpreter(model);
  interpreter_register_polyamri_functions(interp);

  return model;
}

// Main program.
int main(int argc, char* argv[])
{
  return model_main("str_heat", heat_ctor, argc, argv);
}
