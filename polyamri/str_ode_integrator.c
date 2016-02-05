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
  ode_integrator_t* integ;
  real_t* X;
  int X_size;
};

str_ode_integrator_t* str_ode_integrator_new(ode_integrator_t* integrator)
{
  str_ode_integrator_t* integ = polymec_malloc(sizeof(str_ode_integrator_t));
  char name[1025];
  snprintf(name, 1024, "Structured ODE integrator[%s]", ode_integrator_name(integrator));
  integ->name = string_dup(name);
  integ->integ = integrator;
  integ->X = NULL;
  integ->X_size = 0;
  return integ;
}

void str_ode_integrator_free(str_ode_integrator_t* integ)
{
  ode_integrator_free(integ->integ);
  string_free(integ->name);
  if (integ->X != NULL)
    polymec_free(integ->X);
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

static void pack_data(str_ode_integrator_t* integ, 
                      str_grid_cell_data_t* X)
{
  int num_comps = str_grid_cell_data_num_components(X);
  int data_size = str_grid_cell_data_num_cells(X, true) * num_comps;

  // Resize our buffer if needed.
  if (integ->X_size == 0)
  {
    integ->X_size = data_size;
    integ->X = polymec_malloc(sizeof(real_t) * integ->X_size);
  }
  else if (integ->X_size < data_size)
  {
    while (integ->X_size < data_size)
      integ->X_size *= 2;
    integ->X = polymec_realloc(integ->X, sizeof(real_t) * integ->X_size);
  }

  // Copy the contents in.
  str_grid_cell_data_pack(X, integ->X);
}

static void unpack_data(str_ode_integrator_t* integ, 
                        str_grid_cell_data_t* X)
{
  // Copy the data back into the patches.
  str_grid_cell_data_unpack(X, integ->X);
}

bool str_ode_integrator_step(str_ode_integrator_t* integ, 
                             real_t max_dt, 
                             real_t* t, 
                             str_grid_cell_data_t* X)
{
  pack_data(integ, X);
  bool result = ode_integrator_step(integ->integ, max_dt, t, integ->X);
  if (result)
    unpack_data(integ, X);
  return result;
}

bool str_ode_integrator_advance(str_ode_integrator_t* integ, 
                                real_t t1, real_t t2, 
                                str_grid_cell_data_t* X)
{
  pack_data(integ, X);
  bool result = ode_integrator_advance(integ->integ, t1, t2, integ->X);
  if (result)
    unpack_data(integ, X);
  return result;
}

void str_ode_integrator_reset(str_ode_integrator_t* integ, 
                              real_t t, 
                              str_grid_cell_data_t* X)
{
  pack_data(integ, X);
  ode_integrator_reset(integ->integ, t, integ->X);
  unpack_data(integ, X);
}

real_t str_ode_integrator_current_time(str_ode_integrator_t* integ)
{
  return ode_integrator_current_time(integ->integ);
}

