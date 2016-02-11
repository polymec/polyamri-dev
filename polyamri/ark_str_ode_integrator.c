// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/ark_str_ode_integrator.h"

typedef struct
{
  ode_integrator_t* ark; // Underlying ARK integrator.

  // Methods.
  void* context;
  int (*fe_func)(void* context, real_t t, str_grid_cell_data_t* x, str_grid_cell_data_t* fe);
  int (*fi_func)(void* context, real_t t, str_grid_cell_data_t* x, str_grid_cell_data_t* fi);
  real_t (*stable_dt_func)(void* context, real_t t, str_grid_cell_data_t* x);
  void (*dtor)(void* context);

  // Bookkeeping.
  int nc, ng;
  str_grid_cell_data_t* x;
  str_grid_cell_data_t* fi;
  str_grid_cell_data_t* fe;
} str_ark_t;

static int ark_fe_func(void* context, 
                       real_t t, 
                       real_t* x,
                       real_t* fe)
{
  str_ark_t* ark = context;
  str_grid_cell_data_set_buffer(ark->x, x, false);
  str_grid_cell_data_set_buffer(ark->fe, fe, false);
  return ark->fe_func(ark->context, t, ark->x, ark->fe);
}

static int ark_fi_func(void* context, 
                       real_t t, 
                       real_t* x,
                       real_t* fi)
{
  str_ark_t* ark = context;
  str_grid_cell_data_set_buffer(ark->x, x, false);
  str_grid_cell_data_set_buffer(ark->fi, fi, false);
  return ark->fi_func(ark->context, t, ark->x, ark->fi);
}

static real_t ark_stable_dt_func(void* context, 
                                 real_t t, 
                                 real_t* x)
{
  str_ark_t* ark = context;
  str_grid_cell_data_set_buffer(ark->x, x, false);
  return ark->stable_dt_func(ark->context, t, ark->x);
}

static bool ark_step(void* context, real_t max_dt, real_t* t, str_grid_cell_data_t* x)
{
  str_ark_t* ark = context;
  real_t* X = str_grid_cell_data_buffer(x);
  return ode_integrator_step(ark->ark, max_dt, t, X);
}

static bool ark_advance(void* context, real_t t1, real_t t2, str_grid_cell_data_t* x)
{
  str_ark_t* ark = context;
  real_t* X = str_grid_cell_data_buffer(x);
  return ode_integrator_advance(ark->ark, t1, t2, X);
}

static void ark_reset(void* context, real_t t, str_grid_cell_data_t* x)
{
  str_ark_t* ark = context;
  real_t* X = str_grid_cell_data_buffer(x);
  ode_integrator_reset(ark->ark, t, X);
}

static void ark_dtor(void* context)
{
  str_ark_t* ark = context;
  str_grid_cell_data_free(ark->fi);
  str_grid_cell_data_free(ark->fe);
  str_grid_cell_data_free(ark->x);
  ode_integrator_free(ark->ark);
  if ((ark->context != NULL) && (ark->dtor != NULL))
    ark->dtor(ark->context);
  polymec_free(ark);
}

str_ode_integrator_t* explicit_ark_str_ode_integrator_new(int order, 
                                                          MPI_Comm comm,
                                                          str_grid_t* grid,
                                                          int num_components,
                                                          int num_ghost_layers,
                                                          void* context, 
                                                          int (*fe_func)(void* context, real_t t, str_grid_cell_data_t* x, str_grid_cell_data_t* fe),
                                                          real_t (*stable_dt_func)(void* context, real_t t, str_grid_cell_data_t* x),
                                                          void (*dtor)(void* context))
{
  return functional_ark_str_ode_integrator_new(order, comm, grid, num_components, num_ghost_layers, 
                                               context, fe_func, NULL, stable_dt_func, dtor, 0);
}

str_ode_integrator_t* functional_ark_str_ode_integrator_new(int order, 
                                                            MPI_Comm comm,
                                                            str_grid_t* grid,
                                                            int num_components,
                                                            int num_ghost_layers,
                                                            void* context, 
                                                            int (*fe_func)(void* context, real_t t, str_grid_cell_data_t* x, str_grid_cell_data_t* fe),
                                                            int (*fi_func)(void* context, real_t t, str_grid_cell_data_t* x, str_grid_cell_data_t* fi),
                                                            real_t (*stable_dt_func)(void* context, real_t t, str_grid_cell_data_t* x),
                                                            void (*dtor)(void* context),
                                                            int max_anderson_accel_dim)
{
  str_ark_t* ark = polymec_malloc(sizeof(str_ark_t));
  ark->context = context;
  ark->fe_func = fe_func;
  ark->fi_func = fi_func;
  ark->stable_dt_func = stable_dt_func;
  ark->dtor = dtor;
  ark->x = str_grid_cell_data_with_buffer(grid, num_components, num_ghost_layers, NULL);
  ark->fe = str_grid_cell_data_with_buffer(grid, num_components, num_ghost_layers, NULL);
  ark->fi = str_grid_cell_data_with_buffer(grid, num_components, num_ghost_layers, NULL);
  ark->nc = num_components;
  ark->ng = num_ghost_layers;

  // Create the underlying ARK integrator.
  int num_values = num_components * str_grid_cell_data_num_cells(ark->x, true); // includes ghost cells!
  ark->ark = functional_ark_ode_integrator_new(order, 
                                               comm, 
                                               num_values, 
                                               0, 
                                               ark, 
                                               (fe_func != NULL) ? ark_fe_func : NULL, 
                                               (fi_func != NULL) ? ark_fi_func : NULL,
                                               (stable_dt_func) ? ark_stable_dt_func : NULL, 
                                               NULL, 
                                               max_anderson_accel_dim);

  // Pull everything together.
  str_ode_integrator_vtable vtable = {.step = ark_step, .advance = ark_advance, .reset = ark_reset, .dtor = ark_dtor};
  char name[1024];
  if (ark->fi != NULL)
    snprintf(name, 1024, "Additive Runge-Kutta (fixed-point, order %d)", order);
  else
    snprintf(name, 1024, "Explicit Runge-Kutta (order %d)", order);
  str_ode_integrator_t* integ = str_ode_integrator_new(name, ark, vtable, order,
                                                       grid, num_components, num_ghost_layers);
  return integ;
}

str_ode_integrator_t* jfnk_ark_str_ode_integrator_new(int order, 
                                                      MPI_Comm comm,
                                                      str_grid_t* grid,
                                                      int num_components,
                                                      int num_ghost_layers,
                                                      void* context, 
                                                      int (*fe_func)(void* context, real_t t, str_grid_cell_data_t* x, str_grid_cell_data_t* fe),
                                                      int (*fi_func)(void* context, real_t t, str_grid_cell_data_t* x, str_grid_cell_data_t* fi),
                                                      bool fi_is_linear,
                                                      bool fi_is_time_dependent,
                                                      real_t (*stable_dt_func)(void* context, real_t, str_grid_cell_data_t* x),
                                                      void (*dtor)(void* context),
                                                      newton_pc_t* precond,
                                                      jfnk_ark_krylov_t solver_type,
                                                      int max_krylov_dim)
{
  str_ark_t* ark = polymec_malloc(sizeof(str_ark_t));
  ark->context = context;
  ark->fe_func = fe_func;
  ark->fi_func = fi_func;
  ark->stable_dt_func = stable_dt_func;
  ark->dtor = dtor;
  ark->x = str_grid_cell_data_with_buffer(grid, num_components, num_ghost_layers, NULL);
  ark->fe = str_grid_cell_data_with_buffer(grid, num_components, num_ghost_layers, NULL);
  ark->fi = str_grid_cell_data_with_buffer(grid, num_components, num_ghost_layers, NULL);
  ark->nc = num_components;
  ark->ng = num_ghost_layers;

  // Create the underlying ARK integrator.
  int num_values = num_components * str_grid_cell_data_num_cells(ark->x, true); // includes ghost cells!
  ark->ark = jfnk_ark_ode_integrator_new(order, 
                                         comm, 
                                         num_values, 
                                         0, 
                                         ark, 
                                         (fe_func != NULL) ? ark_fe_func : NULL, 
                                         (fi_func != NULL) ? ark_fi_func : NULL,
                                         fi_is_linear, 
                                         fi_is_time_dependent, 
                                         (stable_dt_func) ? ark_stable_dt_func : NULL, 
                                         NULL, 
                                         NULL, 
                                         precond, 
                                         solver_type, 
                                         max_krylov_dim);

  // Pull everything together.
  str_ode_integrator_vtable vtable = {.step = ark_step, .advance = ark_advance, .reset = ark_reset, .dtor = ark_dtor};
  char name[1024];
  if (ark->fe != NULL)
    snprintf(name, 1024, "JFNK IMEX Additive Runge-Kutta (fixed-point, order %d)", order);
  else
    snprintf(name, 1024, "JFNK implicit Runge-Kutta (order %d)", order);
  str_ode_integrator_t* integ = str_ode_integrator_new(name, ark, vtable, order,
                                                       grid, num_components, num_ghost_layers);
  return integ;
}

void ark_str_ode_integrator_set_step_controls(str_ode_integrator_t* integrator,
                                              real_t max_growth,
                                              real_t max_initial_growth,
                                              real_t max_convergence_cut_factor,
                                              real_t max_accuracy_cut_factor,
                                              real_t safety_factor,
                                              real_t cfl_fraction)
{
  str_ark_t* ark = str_ode_integrator_context(integrator);
  ark_ode_integrator_set_step_controls(ark->ark, max_growth, max_initial_growth,
                                       max_convergence_cut_factor, max_accuracy_cut_factor,
                                       safety_factor, cfl_fraction);
}

void ark_str_ode_integrator_set_predictor(str_ode_integrator_t* integrator, 
                                          ark_predictor_t predictor)
{
  str_ark_t* ark = str_ode_integrator_context(integrator);
  ark_ode_integrator_set_predictor(ark->ark, predictor);
}

void ark_str_ode_integrator_set_tolerances(str_ode_integrator_t* integrator,
                                           real_t relative_tol, real_t absolute_tol)
{
  str_ark_t* ark = str_ode_integrator_context(integrator);
  ark_ode_integrator_set_tolerances(ark->ark, relative_tol, absolute_tol);
}

void ark_str_ode_integrator_set_max_err_test_failures(str_ode_integrator_t* integrator,
                                                      int max_failures)
{
  str_ark_t* ark = str_ode_integrator_context(integrator);
  ark_ode_integrator_set_max_err_test_failures(ark->ark, max_failures);
}

void ark_str_ode_integrator_set_max_nonlinear_iterations(str_ode_integrator_t* integrator,
                                                         int max_iterations)
{
  str_ark_t* ark = str_ode_integrator_context(integrator);
  ark_ode_integrator_set_max_nonlinear_iterations(ark->ark, max_iterations);
}

void ark_str_ode_integrator_set_nonlinear_convergence_coeff(str_ode_integrator_t* integrator,
                                                            real_t coefficient)
{
  str_ark_t* ark = str_ode_integrator_context(integrator);
  ark_ode_integrator_set_nonlinear_convergence_coeff(ark->ark, coefficient);
}

void ark_str_ode_integrator_eval_fe(str_ode_integrator_t* integrator, real_t t, str_grid_cell_data_t* X, str_grid_cell_data_t* fe)
{
  str_ark_t* ark = str_ode_integrator_context(integrator);
  ASSERT(str_grid_cell_data_num_components(X) == ark->nc);
  ASSERT(str_grid_cell_data_num_ghost_layers(X) == ark->ng);
  ASSERT(str_grid_cell_data_num_components(fe) == ark->nc);
  ASSERT(str_grid_cell_data_num_ghost_layers(fe) == ark->ng);

  real_t* x = str_grid_cell_data_buffer(X);
  real_t* f = str_grid_cell_data_buffer(fe);
  ark_ode_integrator_eval_fi(ark->ark, t, x, f);
}

void ark_str_ode_integrator_eval_fi(str_ode_integrator_t* integrator, real_t t, str_grid_cell_data_t* X, str_grid_cell_data_t* fi)
{
  str_ark_t* ark = str_ode_integrator_context(integrator);
  ASSERT(str_grid_cell_data_num_components(X) == ark->nc);
  ASSERT(str_grid_cell_data_num_ghost_layers(X) == ark->ng);
  ASSERT(str_grid_cell_data_num_components(fi) == ark->nc);
  ASSERT(str_grid_cell_data_num_ghost_layers(fi) == ark->ng);

  real_t* x = str_grid_cell_data_buffer(X);
  real_t* f = str_grid_cell_data_buffer(fi);
  ark_ode_integrator_eval_fi(ark->ark, t, x, f);
}

newton_pc_t* ark_str_ode_integrator_preconditioner(str_ode_integrator_t* integrator)
{
  str_ark_t* ark = str_ode_integrator_context(integrator);
  return ark_ode_integrator_preconditioner(ark->ark);
}

void ark_str_ode_integrator_get_diagnostics(str_ode_integrator_t* integrator, 
                                            ark_ode_integrator_diagnostics_t* diagnostics)
{
  str_ark_t* ark = str_ode_integrator_context(integrator);
  ark_ode_integrator_get_diagnostics(ark->ark, diagnostics);
}

