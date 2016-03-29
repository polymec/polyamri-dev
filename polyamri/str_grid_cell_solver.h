// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_GRID_CELL_SOLVER_H
#define POLYAMRI_STR_GRID_CELL_SOLVER_H

#include "polyamri/str_grid_cell_data.h"

// A structured grid cell solver computes solutions to linear equations on 
// the cells of an associated structured grid.
typedef struct str_grid_cell_solver_t str_grid_cell_solver_t;

// Every str_grid_cell_solver must override the methods in this virtual table.
typedef struct
{
  bool (*solve)(void* context, str_grid_cell_data_t* X);
  void (*dtor)(void* context);
} str_grid_cell_solver_vtable;

// Creates a new grid solver using the given context and vtable, on the given 
// grid, with a solution vector having the given number of components per cell.
str_grid_cell_solver_t* str_grid_cell_solver_new(const char* solver_name,
                                                 void* context,
                                                 str_grid_cell_solver_vtable vtable,
                                                 str_grid_t* grid,
                                                 int num_comps);

// Destroys the given grid solver.
void str_grid_cell_solver_free(str_grid_cell_solver_t* solver);

// Returns an internal pointer to the underlying grid.
str_grid_t* str_grid_cell_solver_grid(str_grid_cell_solver_t* solver);

// Returns the number of components for the solution vector on each cell.
int str_grid_cell_solver_num_comps(str_grid_cell_solver_t* solver);

// Solves the linear system underlying the given grid solver, returning true 
// on success and false on failure, and placing the solution in X.
bool str_grid_cell_solver_solve(str_grid_cell_solver_t* solver, 
                                str_grid_cell_data_t* X);

#endif

