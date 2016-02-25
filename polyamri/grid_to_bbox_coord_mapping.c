// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/grid_to_bbox_coord_mapping.h"

// Forward mapping.
static void g2b_map_point(void* context, point_t* x, point_t* y)
{
  bbox_t* bbox = context;
  y->x = bbox->x1 + x->x * (bbox->x2 - bbox->x1);
  y->y = bbox->y1 + x->y * (bbox->y2 - bbox->y1);
  y->z = bbox->z1 + x->z * (bbox->z2 - bbox->z1);
}

static void g2b_map_vector(void* context, point_t* x, vector_t* v, vector_t* v1)
{
  bbox_t* bbox = context;
  v1->x = v->x * (bbox->x2 - bbox->x1);
  v1->y = v->y * (bbox->y2 - bbox->y1);
  v1->z = v->z * (bbox->z2 - bbox->z1);
}

static void g2b_jacobian(void* context, point_t* x, real_t* J)
{
  bbox_t* bbox = context;
  J[0] = bbox->x2 - bbox->x1;
  J[4] = bbox->y2 - bbox->y1;
  J[8] = bbox->z2 - bbox->z1;
  J[1] = J[2] = J[3] = J[5] = J[6] = J[7] = 0.0; // off-diagonals.
}

// Inverse mapping.
static void b2g_map_point(void* context, point_t* y, point_t* x)
{
  bbox_t* bbox = context;
  x->x = (y->x - bbox->x1) / (bbox->x2 - bbox->x1);
  x->y = (y->y - bbox->y1) / (bbox->y2 - bbox->y1);
  x->z = (y->z - bbox->z1) / (bbox->z2 - bbox->z1);
}

static void b2g_map_vector(void* context, point_t* x, vector_t* v, vector_t* v1)
{
  bbox_t* bbox = context;
  v1->x = v->x / (bbox->x2 - bbox->x1);
  v1->y = v->y / (bbox->y2 - bbox->y1);
  v1->z = v->z / (bbox->z2 - bbox->z1);
}

static void b2g_jacobian(void* context, point_t* x, real_t* J)
{
  bbox_t* bbox = context;
  J[0] = 1.0 / (bbox->x2 - bbox->x1);
  J[4] = 1.0 / (bbox->y2 - bbox->y1);
  J[8] = 1.0 / (bbox->z2 - bbox->z1);
  J[1] = J[2] = J[3] = J[5] = J[6] = J[7] = 0.0; // off-diagonals.
}

// Creation of inverse mappings.
static coord_mapping_t* g2b_inverse(void* context)
{
  bbox_t* bbox = context;
  char name[1024];
  snprintf(name, 1023, "[%g,%g]x[%g,%g]x[%g,%g] -> [0,1]x[0,1]x[0,1]",
           bbox->x1, bbox->x2, bbox->y1, bbox->y2, bbox->z1, bbox->z2);
  coord_mapping_vtable vtable = {.map_point = b2g_map_point,
                                 .map_vector = b2g_map_vector,
                                 .jacobian = b2g_jacobian,
                                 .inverse = (coord_mapping_inverse_func)grid_to_bbox_coord_mapping_new};
  return coord_mapping_new(name, bbox_clone(bbox), vtable);
}

coord_mapping_t* grid_to_bbox_coord_mapping_new(bbox_t* bbox)
{
  char name[1024];
  snprintf(name, 1023, "[0,1]x[0,1]x[0,1] -> [%g,%g]x[%g,%g]x[%g,%g]",
           bbox->x1, bbox->x2, bbox->y1, bbox->y2, bbox->z1, bbox->z2);
  coord_mapping_vtable vtable = {.map_point = g2b_map_point,
                                 .map_vector = g2b_map_vector,
                                 .jacobian = g2b_jacobian,
                                 .inverse = g2b_inverse};
  return coord_mapping_new(name, bbox_clone(bbox), vtable);
}

