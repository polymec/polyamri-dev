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

// Writes the given AMR patch to the given Silo file, optionally using 
// the given bounding box to provide coordinates. If bbox is NULL, 
// the patch will sit in the space [0,1] x [0,1] x [0,1].
void silo_file_write_amr_patch(silo_file_t* file, 
                               const char* patch_name,
                               amr_patch_t* patch,
                               bbox_t* bbox);

// Writes the given AMR patch to the given Silo file, using the given
// spatial function to provide a mapping from [0,1] x [0,1] x [0,1] to 
// a deformed region of space.
void silo_file_write_mapped_amr_patch(silo_file_t* file, 
                                      const char* patch_name,
                                      amr_patch_t* patch,
                                      sp_func_t* mapping);

// Writes the given AMR grid to the given Silo file.
void silo_file_write_amr_grid(silo_file_t* file, 
                              const char* grid_name,
                              amr_grid_t* grid);

// Writes the given AMR grid data to the given Silo file.
void silo_file_write_amr_grid_data(silo_file_t* file, 
                                   const char* data_name,
                                   amr_grid_data_t* grid_data,
                                   silo_field_metadata_t** field_metadata);

// Writes the given AMR grid hierarchy to the given Silo file.
void silo_file_write_amr_grid_hierarchy(silo_file_t* file, 
                                        const char* hierarchy_name,
                                        amr_grid_hierarchy_t* hierarchy);

// Writes the given AMR data hierarchy to the given Silo file.
void silo_file_write_amr_data_hierarchy(silo_file_t* file, 
                                        const char* hierarchy_name,
                                        amr_data_hierarchy_t* hierarchy,
                                        silo_field_metadata_t** field_metadata);

#endif

