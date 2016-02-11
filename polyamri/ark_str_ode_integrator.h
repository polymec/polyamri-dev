// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_ARK_STR_ODE_INTEGRATOR_H
#define POLYAMRI_ARK_STR_ODE_INTEGRATOR_H

#include "integrators/ark_ode_integrator.h"
#include "polyamri/str_ode_integrator.h"

// This type of structured ODE integrator integrates a multi-scale set of 
// ordinary differential equations using semi-implicit Additive Runge Kutte 
// (ARK) methods. The equations integrated take the form
//
// dx/dt = fe(t, x) + fi(t, x)
//
// where fe is a slowly-varying function that can be integrated explicitly, 
// and fi is a quickly-varying function that should be integrated implicitly.
// This is based on polymec's ARK integrators and has the same features.

// This constructs an integrator for the slowly-varying system
// dx/dt = fe(t, x) with the desired order of time accuracy.
// The optional stable_dt_func argument supplies a function that computes the 
// next timestep subject to stability constraints.
str_ode_integrator_t* explicit_ark_str_ode_integrator_new(int order, 
                                                          MPI_Comm comm,
                                                          str_grid_t* grid,
                                                          int num_components,
                                                          int num_ghost_layers,
                                                          void* context, 
                                                          int (*fe_func)(void* context, real_t t, str_grid_cell_data_t* x, str_grid_cell_data_t* fe),
                                                          real_t (*stable_dt_func)(void* context, real_t, str_grid_cell_data_t* x),
                                                          void (*dtor)(void* context));

// This constructs an integrator for the system with given fe, fi,
// with the desired order of time accuracy, using a fixed-point (functional) 
// iteration method that does not require a Newton-Krylov solver. At least 
// one of fe and fi must be non-NULL. If fi is NULL, the system is assumed to 
// vary "slowly" and will be explicitly integrated; if fe is NULL, the system 
// is assumed to be stiff, and will be implicitly integrated; if both are
// non-NULL, they will be integrated using an adaptive IMEX method.
// max_anderson_accel_dim is the maximum dimension of the underlying Anderson 
// acceleration subspace.
str_ode_integrator_t* functional_ark_str_ode_integrator_new(int order, 
                                                            MPI_Comm comm,
                                                            str_grid_t* grid,
                                                            int num_components,
                                                            int num_ghost_layers,
                                                            void* context, 
                                                            int (*fe_func)(void* context, real_t t, str_grid_cell_data_t* x, str_grid_cell_data_t* fe),
                                                            int (*fi_func)(void* context, real_t t, str_grid_cell_data_t* x, str_grid_cell_data_t* fi),
                                                            real_t (*stable_dt_func)(void* context, real_t, str_grid_cell_data_t* x),
                                                            void (*dtor)(void* context),
                                                            int max_anderson_accel_dim);

// This constructs an integrator for the system with given fe, fi,
// with the desired order of time accuracy, using the Jacobian-Free 
// Newton-Krylov solver of the given type. The Jacobian-vector 
// product will be approximated using difference quotients. fi must 
// be non-NULL, and fe can be either non-NULL or NULL. If fe is NULL, 
// the system is assumed to be stiff, and will be implicitly integrated. 
// The fi_is_linear and fi_is_time_dependent flags provide information 
// about fi in order to help the integrator optimize its performance 
// (and fi_is_time_dependent is ignored if fi is not linear in x).
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
                                                      int max_krylov_dim);

// Sets control parameters for adaptive time stepping. Specifically:
// max_growth (> 1) - The maximum factor by which the step is allowed to grow 
//                    between consecutive steps.
// max_initial_growth (> 1) - The maximum factor by which the step is allowed to grow 
//                            after the first step.
// max_convergence_cut_factor (< 1) - The maximum factor by which the time step is 
//                                    modified after successive convergence failures.
// max_accuracy_cut_factor (< 1) - The maximum factor by which the time step is 
//                                 modified after successive accuracy failures.
// safety_factor - The factor by which the accuracy-based time step is multiplied
//                 by the integrator in taking a step.
// cfl_fraction (<= 1) - The factor by which the stability-based time step is multiplied
//                       by the integrator in taking a step (ignored if there's no fe).
void ark_str_ode_integrator_set_step_controls(str_ode_integrator_t* integrator,
                                              real_t max_growth,
                                              real_t max_initial_growth,
                                              real_t max_convergence_cut_factor,
                                              real_t max_accuracy_cut_factor,
                                              real_t safety_factor,
                                              real_t cfl_fraction);

// Sets the method to use for predicting implicit terms in the RK integration.
void ark_str_ode_integrator_set_predictor(str_ode_integrator_t* integrator, 
                                          ark_predictor_t predictor);

// Sets the relative and absolute tolerances for integrated quantities.
void ark_str_ode_integrator_set_tolerances(str_ode_integrator_t* integrator,
                                           real_t relative_tol, real_t absolute_tol);

// Sets the maximum number of error test failures permitted in attempting 
// a single time step. By default, this value is 7.
void ark_str_ode_integrator_set_max_err_test_failures(str_ode_integrator_t* integrator,
                                                      int max_failures);

// Sets the maximum number of nonlinear solver iterations per time step.
// By default, this value is 3.
void ark_str_ode_integrator_set_max_nonlinear_iterations(str_ode_integrator_t* integrator,
                                                         int max_iterations);

// Sets the safety factor (coefficient) used in the nonlinear convergence test.
// By default, this value is 0.1.
void ark_str_ode_integrator_set_nonlinear_convergence_coeff(str_ode_integrator_t* integrator,
                                                            real_t coefficient);

// Evaluates the explicit part of the system at the given time and with the 
// given solution X, placing the results in fe.
void ark_str_ode_integrator_eval_fe(str_ode_integrator_t* integrator, real_t t, str_grid_cell_data_t* X, str_grid_cell_data_t* fe);

// Evaluates the implicit part of the system at the given time and with the 
// given solution X, placing the results in fi.
void ark_str_ode_integrator_eval_fi(str_ode_integrator_t* integrator, real_t t, str_grid_cell_data_t* X, str_grid_cell_data_t* fi);

// Returns an internal pointer to the preconditioner passed to this 
// integrator during construction time.
newton_pc_t* ark_str_ode_integrator_preconditioner(str_ode_integrator_t* integrator);

// Retrieve diagnostics for the time integrator. Note that this uses the same
// diagnostics class as ark_ode_integrator.
void ark_str_ode_integrator_get_diagnostics(str_ode_integrator_t* integrator, 
                                            ark_ode_integrator_diagnostics_t* diagnostics);

#endif

