// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/str_ode_integrator.h"

struct str_ode_integrator_t 
{
  ode_integrator_t* integ;
};

str_ode_integrator_t* str_ode_integrator_new(ode_integrator_t* integrator)
{
}

void str_ode_integrator_free(str_ode_integrator_t* integ)
{
}

char* str_ode_integrator_name(str_ode_integrator_t* integ)
{
}

ode_integrator_t* str_ode_integrator_underlying(str_ode_integrator_t* integ)
{
}

void str_ode_integrator_set_max_dt(str_ode_integrator_t* integ, real_t max_dt)
{
}

void str_ode_integrator_set_stop_time(str_ode_integrator_t* integ, real_t stop_time)
{
}

bool str_ode_integrator_step(str_ode_integrator_t* integ, 
                             real_t max_dt, 
                             real_t* t, 
                             str_grid_cell_data_t* X)
{
}

bool str_ode_integrator_advance(str_ode_integrator_t* integ, 
                                real_t t1, real_t t2, 
                                str_grid_cell_data_t* X)
{
}

void str_ode_integrator_reset(str_ode_integrator_t* integ, 
                              real_t t, 
                              str_grid_cell_data_t* X)
{
}

real_t str_ode_integrator_current_time(str_ode_integrator_t* integ)
{
}

