// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/str_ode_integrator.h"

struct str_ode_integrator_t 
{
  char* name;
  str_grid_t* grid;
  int nc, ng;
  int npx, npy, npz;

  ode_integrator_t* integ; // underlying ODE integrator

  int X_size;
};

static void str_copy_in(void* context, real_t* solution_data, real_t* x)
{
  str_ode_integrator_t* integ = context;
printf("fake copy in\n");
  memcpy(x, solution_data, sizeof(real_t) * integ->X_size);
}

static void str_copy_out(void* context, real_t* x, real_t* solution_data)
{
  str_ode_integrator_t* integ = context;
printf("fake copy out\n");
  memcpy(solution_data, x, sizeof(real_t) * integ->X_size);
}

str_ode_integrator_t* str_ode_integrator_new(ode_integrator_t* integrator,
                                             str_grid_t* grid,
                                             int num_comps,
                                             int num_ghost_layers)
{
  ASSERT(num_comps > 0);
  ASSERT(num_ghost_layers >= 0);

  str_ode_integrator_t* integ = polymec_malloc(sizeof(str_ode_integrator_t));
  char name[1025];
  snprintf(name, 1024, "Structured ODE integrator[%s]", ode_integrator_name(integrator));
  integ->name = string_dup(name);
  integ->grid = grid;
  integ->nc = num_comps;
  integ->ng = num_ghost_layers;

  // Count up patch cells.
  int num_patches = str_grid_num_patches(integ->grid);
  str_grid_get_patch_size(integ->grid, &integ->npx, &integ->npy, &integ->npz);
  int num_patch_cells = (integ->npx + 2*integ->ng) * (integ->npy + 2*integ->ng) * (integ->npz + 2*integ->ng);
  integ->X_size = integ->nc * num_patch_cells * num_patches;
  ASSERT(integ->X_size == ode_integrator_solution_vector_size(integrator));

  // Set up a "decorated integrator" that adds special copy_in / copy_out 
  // methods to the base one.
  ode_integrator_vtable vtable = {.copy_in = str_copy_in,
                                  .copy_out = str_copy_out};
  integ->integ = decorated_ode_integrator_new(integrator, integ, vtable);

  return integ;
}

void str_ode_integrator_free(str_ode_integrator_t* integ)
{
  ode_integrator_free(integ->integ);
  string_free(integ->name);
  polymec_free(integ);
}

char* str_ode_integrator_name(str_ode_integrator_t* integ)
{
  return integ->name;
}

ode_integrator_t* str_ode_integrator_underlying(str_ode_integrator_t* integ)
{
  return integ->integ;
}

void str_ode_integrator_set_max_dt(str_ode_integrator_t* integ, real_t max_dt)
{
  ode_integrator_set_max_dt(integ->integ, max_dt);
}

void str_ode_integrator_set_stop_time(str_ode_integrator_t* integ, real_t stop_time)
{
  ode_integrator_set_stop_time(integ->integ, stop_time);
}

bool str_ode_integrator_step(str_ode_integrator_t* integ, 
                             real_t max_dt, 
                             real_t* t, 
                             str_grid_cell_data_t* X)
{
  ASSERT(str_grid_cell_data_grid(X) == integ->grid);
  ASSERT(str_grid_cell_data_num_components(X) == integ->nc);
  ASSERT(str_grid_cell_data_num_ghost_layers(X) == integ->ng);

  real_t* x = str_grid_cell_data_buffer(X);
  bool result = ode_integrator_step(integ->integ, max_dt, t, x);
  return result;
}

bool str_ode_integrator_advance(str_ode_integrator_t* integ, 
                                real_t t1, real_t t2, 
                                str_grid_cell_data_t* X)
{
  ASSERT(str_grid_cell_data_grid(X) == integ->grid);
  ASSERT(str_grid_cell_data_num_components(X) == integ->nc);
  ASSERT(str_grid_cell_data_num_ghost_layers(X) == integ->ng);

  real_t* x = str_grid_cell_data_buffer(X);
  bool result = ode_integrator_advance(integ->integ, t1, t2, x);
  return result;
}

void str_ode_integrator_reset(str_ode_integrator_t* integ, 
                              real_t t, 
                              str_grid_cell_data_t* X)
{
  ASSERT(str_grid_cell_data_grid(X) == integ->grid);
  ASSERT(str_grid_cell_data_num_components(X) == integ->nc);
  ASSERT(str_grid_cell_data_num_ghost_layers(X) == integ->ng);

  real_t* x = str_grid_cell_data_buffer(X);
  ode_integrator_reset(integ->integ, t, x);
}

real_t str_ode_integrator_current_time(str_ode_integrator_t* integ)
{
  return ode_integrator_current_time(integ->integ);
}

