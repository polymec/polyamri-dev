// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/hypre_str_grid_cell_solvers.h"

typedef struct
{
  int hi_there;
} smg_t;

static bool smg_solve(void* context, str_grid_cell_data_t* X)
{
  return false;
}

str_grid_cell_solver_t* hypre_smg_helmholtz_str_grid_cell_solver_new(str_grid_t* grid,
                                                                     int num_comps)
{
  smg_t* smg = polymec_malloc(sizeof(smg_t));
  str_grid_cell_solver_vtable vtable = {.solve = smg_solve, .dtor = polymec_free};
  return str_grid_cell_solver_new("HYPRE SMG Helmholtz", 
                                  smg, vtable, grid, num_comps);
}

void hypre_smg_helmholtz_str_grid_cell_solver_set_operator(str_grid_cell_solver_t* solver,
                                                           real_t time, 
                                                           real_t alpha, 
                                                           real_t beta, 
                                                           st_func_t* A)
{
}

void hypre_smg_helmholtz_str_grid_cell_solver_set_rhs(str_grid_cell_solver_t* solver, 
                                                      real_t time,
                                                      real_t gamma, 
                                                      real_t delta, 
                                                      st_func_t* A, 
                                                      st_func_t* B, 
                                                      real_t epsilon, 
                                                      st_func_t* C)
{
}

