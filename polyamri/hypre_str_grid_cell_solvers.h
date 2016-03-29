// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_HYPRE_STR_GRID_CELL_SOLVERS_H
#define POLYAMRI_HYPRE_STR_GRID_CELL_SOLVERS_H

#include "core/st_func.h"
#include "polyamri/str_grid_cell_solver.h"

//------------------------------------------------------------------------
//                   HYPRE linear solvers for structured grids
//------------------------------------------------------------------------
// All of these functions create HYPRE-enabled linear solvers. If the HYPRE 
// library exists on the system within the directory given in hypre_path, and 
// is properly installed, each function returns a pointer to a newly allocated 
// HYPRE-enabled str_grid_cell_solver. Otherwise, the function returns NULL.
//------------------------------------------------------------------------

// Creates a solver that can solve a Helmholtz-like equation on the cells of 
// a structured grid using HYPRE's SMG (simple multigrid) solver. The form of 
// the linear system is 
//
// (alpha * I + beta * A) X = (gamma * I + delta * A) B + epsilon * C
// 
// where A is a linear operator, I is the identity operator, B and C are 
// vectors, alpha, beta, gamma, delta, epsilon are scalar coefficients, and 
// X is the solution vector.
str_grid_cell_solver_t* hypre_smg_helmholtz_str_grid_cell_solver_new(const char* hypre_dir,
                                                                     str_grid_t* grid,
                                                                     int num_comps);

// Sets the coefficients and operator information for the left hand side of 
// the Helmholtz equation.
void hypre_smg_helmholtz_str_grid_cell_solver_set_operator(str_grid_cell_solver_t* solver,
                                                           real_t time, 
                                                           real_t alpha, 
                                                           real_t beta, 
                                                           st_func_t* A);

// Sets the coefficients and operator information for the right hand side of 
// the Helmholtz equation.
void hypre_smg_helmholtz_str_grid_cell_solver_set_rhs(str_grid_cell_solver_t* solver, 
                                                      real_t time,
                                                      real_t gamma, 
                                                      real_t delta, 
                                                      st_func_t* A, 
                                                      st_func_t* B, 
                                                      real_t epsilon, 
                                                      st_func_t* C);

#endif

