
// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_INTERPRETER_REGISTER_POLYAMRI_FUNCTIONS_H
#define POLYAMRI_INTERPRETER_REGISTER_POLYAMRI_FUNCTIONS_H

#include "model/interpreter.h"
#include "polyamri/str_grid.h"
#include "polyamri/str_grid_patch_filler.h"

// This file contains additional functions for polymec's interpreter 
// that augment its capabilities for use with polyamri.

// Register the polyamri-specific constructors with the given interpreter.
void interpreter_register_polyamri_functions(interpreter_t* interp);

// Fetches the given structured grid from the interpreter, returning NULL if 
// it is not found or if the given variable is not a structured grid.
str_grid_t* interpreter_get_str_grid(interpreter_t* interp, const char* name);

// Fetches the given patch filler from the interpreter, returning NULL if it 
// is not found or if the given variable is not a patch filler. Patch fillers 
// can be used to implement boundary conditions, for example.
str_grid_patch_filler_t* interpreter_get_str_grid_patch_filler(interpreter_t* interp, const char* name);

#endif
