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
  void* context;
  str_ode_integrator_vtable vtable;

  int order;
  str_grid_t* grid;
  int nc, ng;

  bool initialized;
  real_t current_time;
  real_t max_dt, stop_time;
};

str_ode_integrator_t* str_ode_integrator_new(const char* name, 
                                             void* context,
                                             str_ode_integrator_vtable vtable,
                                             int order,
                                             str_grid_t* grid,
                                             int num_components,
                                             int num_ghost_layers)
{
  ASSERT(order > 0);
  ASSERT(num_components > 0);
  ASSERT(num_ghost_layers >= 0);

  str_ode_integrator_t* integ = polymec_malloc(sizeof(str_ode_integrator_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->vtable = vtable;
  integ->order = order;
  integ->grid = grid;
  integ->nc = num_components;
  integ->ng = num_ghost_layers;

  integ->current_time = 0.0;
  integ->initialized = false;
  integ->max_dt = FLT_MAX;
  integ->stop_time = FLT_MAX;

  return integ;
}

void str_ode_integrator_free(str_ode_integrator_t* integ)
{
  if ((integ->context != NULL) && (integ->vtable.dtor != NULL))
    integ->vtable.dtor(integ->context);
  string_free(integ->name);
  polymec_free(integ);
}

char* str_ode_integrator_name(str_ode_integrator_t* integ)
{
  return integ->name;
}

void* str_ode_integrator_context(str_ode_integrator_t* integ)
{
  return integ->context;
}

int str_ode_integrator_order(str_ode_integrator_t* integ)
{
  return integ->order;
}

void str_ode_integrator_set_max_dt(str_ode_integrator_t* integ, real_t max_dt)
{
  ASSERT(max_dt > 0.0);
  integ->max_dt = max_dt;
}

void str_ode_integrator_set_stop_time(str_ode_integrator_t* integ, real_t stop_time)
{
  integ->stop_time = stop_time;
}

bool str_ode_integrator_step(str_ode_integrator_t* integ, 
                             real_t max_dt, 
                             real_t* t, 
                             str_grid_cell_data_t* x)
{
  ASSERT(str_grid_cell_data_grid(x) == integ->grid);
  ASSERT(str_grid_cell_data_num_components(x) == integ->nc);
  ASSERT(str_grid_cell_data_num_ghost_layers(x) == integ->ng);

  // Figure out the actual maximum time.
  real_t dt = MIN(max_dt, MIN(integ->max_dt, integ->stop_time - *t));

  // Integrate.
  return integ->vtable.step(integ->context, dt, t, x);
}

bool str_ode_integrator_advance(str_ode_integrator_t* integ, 
                                real_t t1, real_t t2, 
                                str_grid_cell_data_t* x)
{
  ASSERT(str_grid_cell_data_grid(x) == integ->grid);
  ASSERT(str_grid_cell_data_num_components(x) == integ->nc);
  ASSERT(str_grid_cell_data_num_ghost_layers(x) == integ->ng);

  // Figure out the actual end time.
  t2 = MIN(t2, integ->stop_time);

  // Advance.
  bool result = integ->vtable.advance(integ->context, t1, t2, x);

  // After a full integration, we must be reset.
  integ->initialized = false;

  return result;
}

void str_ode_integrator_reset(str_ode_integrator_t* integ, 
                              real_t t, 
                              str_grid_cell_data_t* x)
{
  ASSERT(str_grid_cell_data_grid(x) == integ->grid);
  ASSERT(str_grid_cell_data_num_components(x) == integ->nc);
  ASSERT(str_grid_cell_data_num_ghost_layers(x) == integ->ng);

  integ->current_time = t;
  if (integ->vtable.reset != NULL)
    integ->vtable.reset(integ->context, t, x);

  integ->initialized = true;
}

real_t str_ode_integrator_current_time(str_ode_integrator_t* integ)
{
  return integ->current_time;
}

