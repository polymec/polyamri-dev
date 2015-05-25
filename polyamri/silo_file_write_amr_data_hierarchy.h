// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_AMR_SILO_FILE_H
#define POLYAMRI_AMR_SILO_FILE_H

#include "core/silo_file.h"
#include "polyamri/amr_data_hierarchy.h"

// Writes the given AMR data hierarchy to a the given Silo file.
void silo_file_write_amr_data_hierarchy(silo_file_t* file, 
                                        const char* hierarchy_name,
                                        amr_data_hierarchy_t* hierarchy,
                                        silo_field_metadata_t** field_metadata);

#endif

