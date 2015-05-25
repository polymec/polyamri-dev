// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_WRITE_SILO_DATA_HIERARCHY_H
#define POLYAMRI_WRITE_SILO_DATA_HIERARCHY_H

#include "silo.h"
#include "polyamri/data_hierarchy.h"

// Writes the given data hierarchy to a the given Silo file.
void write_silo_data_hierarchy(DBfile* file, data_hierarchy_t* hierarchy);

#endif

