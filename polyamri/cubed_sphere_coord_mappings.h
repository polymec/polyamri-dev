// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_GNOMONIC_CUBED_SPHERE_COORD_MAPPING_H
#define POLYAMRI_GNOMONIC_CUBED_SPHERE_COORD_MAPPING_H

#include "geometry/coord_mapping.h"

// Returns a coordinate mapping from [0,1]x[0,1]x[0,1] to gnomonic coordinates 
// (r, X, Y, index) within a spherical shell with inner radius r1 and outer 
// radius r2, and the particular mapping determined by the given index. 
// Possible index values are:
// 0-3: one of 4 "equatorial" maps
// 4: a map whose center is the North pole
// 5: a map whose center is the South pole
coord_mapping_t* gnomonic_cubed_sphere_coord_mapping_new(real_t r1,
                                                         real_t r2,
                                                         int index);

// Returns a coordinate mapping from [0,1]x[0,1]x[0,1] to equiangular 
// coordinates (r, alpha, beta, index) within a spherical shell with inner
// radius r1 and outer radius r2, and the particular mapping determined by 
// the given index.
coord_mapping_t* equiangular_cubed_sphere_coord_mapping_new(real_t r1,
                                                            real_t r2,
                                                            int index);

// Returns a coordinate mapping from [0,1]x[0,1]x[0,1] to spherical 
// coordinates (r, lambda, theta, index) within a spherical shell with inner
// radius r1 and outer radius r2, and the particular mapping determined by 
// the given index. Here, lambda is longitude and theta is latitude.
// The sphere is assumed to be centered at (0, 0, 0) in Cartesian coordinates.
coord_mapping_t* spherical_cubed_sphere_coord_mapping_new(real_t r1,
                                                          real_t r2,
                                                          int index);

#endif

