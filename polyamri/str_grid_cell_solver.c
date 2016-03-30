// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/str_grid_cell_solver.h"

struct str_grid_cell_solver_t 
{
  char* name;
  MPI_Comm comm;
  void* context;
  str_grid_cell_solver_vtable vtable;
  str_grid_t* grid;
  int num_comps;
};

str_grid_cell_solver_t* str_grid_cell_solver_new(const char* solver_name,
                                                 MPI_Comm comm,
                                                 void* context,
                                                 str_grid_cell_solver_vtable vtable,
                                                 str_grid_t* grid,
                                                 int num_comps)
{
  ASSERT(vtable.solve != NULL);
  ASSERT(num_comps > 0);
  str_grid_cell_solver_t* solver = polymec_malloc(sizeof(str_grid_cell_solver_t));
  solver->name = string_dup(solver_name);
  solver->comm = comm;
  solver->context = context;
  solver->vtable = vtable;
  solver->grid = grid;
  solver->num_comps = num_comps;
  return solver;
}

void str_grid_cell_solver_free(str_grid_cell_solver_t* solver)
{
  if ((solver->vtable.dtor != NULL) && (solver->context != NULL))
    solver->vtable.dtor(solver->context);
  string_free(solver->name);
  polymec_free(solver);
}

str_grid_t* str_grid_cell_solver_grid(str_grid_cell_solver_t* solver)
{
  return solver->grid;
}

int str_grid_cell_solver_num_comps(str_grid_cell_solver_t* solver)
{
  return solver->num_comps;
}

bool str_grid_cell_solver_solve(str_grid_cell_solver_t* solver, 
                                str_grid_cell_data_t* X)
{
  return solver->vtable.solve(solver->context, X);
}

