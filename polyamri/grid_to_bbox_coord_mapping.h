// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_GRID_TO_BBOX_COORD_MAPPING_H
#define POLYAMRI_GRID_TO_BBOX_COORD_MAPPING_H

#include "geometry/coord_mapping.h"

// Returns a coordinate mapping from [0,1]x[0,1]x[0,1] to a given bounding box.
coord_mapping_t* grid_to_bbox_coord_mapping_new(bbox_t* bbox);

#endif

