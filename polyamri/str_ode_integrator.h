// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_STR_ODE_INTEGRATOR_H
#define POLYAMRI_STR_ODE_INTEGRATOR_H

#include "integrators/ode_integrator.h"
#include "polyamri/str_grid_cell_data.h"

// This class is an adaptor for polymec's ode_integrator class that allows
// cell-centered data on structured grids to be integrated.
typedef struct str_ode_integrator_t str_ode_integrator_t;

// Creates a str_ode_integrator that uses the given ode_integrator instance
// to integrate its cell-centered data. The underlying integrator is consumed
// by this new structured integrator.
str_ode_integrator_t* str_ode_integrator_new(ode_integrator_t* integrator);

// Frees the structured ODE integrator.
void str_ode_integrator_free(str_ode_integrator_t* integ);

// Returns the name of the structured integrator (internally stored).
char* str_ode_integrator_name(str_ode_integrator_t* integ);

// Returns the underlying ODE integrator.
ode_integrator_t* str_ode_integrator_underlying(str_ode_integrator_t* integ);

// Sets the maximum time step size for the next integration step.
void str_ode_integrator_set_max_dt(str_ode_integrator_t* integ, real_t max_dt);

// Sets the time past which the integrator will not step.
void str_ode_integrator_set_stop_time(str_ode_integrator_t* integ, real_t stop_time);

// Integrates the given cell-centered grid data X in place, taking a single step with 
// a maximum size of max_dt, starting at time *t and storing the new time there as well. 
// Returns true if the step succeeded, false if it failed for some reason. 
// If a step fails, both *t and X remain unchanged.
bool str_ode_integrator_step(str_ode_integrator_t* integ, 
                             real_t max_dt, 
                             real_t* t, 
                             str_grid_cell_data_t* X);

// Integrates the given cell-centered grid data X in place from time t1 to t2, taking 
// as many steps as are needed. The amount of work done by this method depends on the 
// number of steps taken, so do not call it in contexts that require a 
// predictable amount of work. Returns true if the integration succeeded, 
// false if it failed for some reason. If a step fails, X remains unchanged.
bool str_ode_integrator_advance(str_ode_integrator_t* integ, 
                                real_t t1, real_t t2, 
                                str_grid_cell_data_t* X);

// Resets the internal state of the integrator to use the given cell-centered grid
// data X at the given time t. It is necessary to call this function when the 
// solution data has been altered since the last step.
void str_ode_integrator_reset(str_ode_integrator_t* integ, 
                              real_t t, 
                              str_grid_cell_data_t* X);

// Returns the current time at which the integrator sits.
real_t str_ode_integrator_current_time(str_ode_integrator_t* integ);

#endif

