// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_SILO_FILE_H
#define POLYAMRI_SILO_FILE_H

#include "core/silo_file.h"
#include "geometry/coord_mapping.h"
#include "polyamri/str_grid_cell_data.h"
#include "polyamri/str_grid_face_data.h"

// Writes the given structured grid patch to the given Silo file. If a non-NULL 
// mapping is provided, the coordinates of the patch will be mapped from 
// [0,1] x [0,1] x [0,1] to a deformed region of space. This is mostly useful 
// for debugging purposes, since AMR grids write patches out using different 
// machinery. If data_mapping is non-NULL, data will also be mapped.
void silo_file_write_str_grid_patch(silo_file_t* file, 
                                    const char** field_component_names,
                                    const char* patch_name,
                                    str_grid_patch_t* patch,
                                    coord_mapping_t* mapping,
                                    coord_mapping_t* data_mapping);

// Returns true if the Silo file contains a structured grid patch with the 
// given name, false if not. 
bool silo_file_contains_str_grid_patch(silo_file_t* file, 
                                       const char* patch_name);

// Writes the given structured grid to the given Silo file. The grid's 
// coordinates will be mapped from [0,1] x [0,1] x [0,1] if a non-NULL mapping 
// is given.
void silo_file_write_str_grid(silo_file_t* file, 
                              const char* grid_name,
                              str_grid_t* grid,
                              coord_mapping_t* mapping);

// Returns true if the Silo file contains a structured grid with the given name, 
// false if not. 
bool silo_file_contains_str_grid(silo_file_t* file, 
                                 const char* grid_name);

// Writes the given cell-centered structured grid data to the given Silo file, 
// associating it with the entry for the grid with the given name. If a 
// non-NULL mapping is given, the data will be mapped accordingly.
void silo_file_write_str_grid_cell_data(silo_file_t* file, 
                                        const char** field_component_names,
                                        const char* grid_name,
                                        str_grid_cell_data_t* cell_data,
                                        silo_field_metadata_t** field_metadata,
                                        coord_mapping_t* mapping);

// Returns true if the Silo file contains cell-centered structured grid data 
// with the given name, associated with the structured grid with the given 
// name, and false if not. 
bool silo_file_contains_str_grid_cell_data(silo_file_t* file, 
                                           const char* cell_data_name,
                                           const char* grid_name);

#endif

