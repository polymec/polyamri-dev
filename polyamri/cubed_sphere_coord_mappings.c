// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/cubed_sphere_coord_mappings.h"

// Cubed sphere mapping type.
typedef struct
{
  real_t r1, r2;
  int index;
} cs_t;

static void spherical_map_point_0(void* context, point_t* x, point_t* y)
{
  cs_t* cs = context;
  real_t X = -1.0 + 2.0 * x->x;
  real_t Y = -1.0 + 2.0 * x->y;
  x->x = atan(X) + 0.5 * M_PI * (cs->index - 1);  // longitude
  x->y = atan2(Y / sqrt(1.0 + X*X));              // latitude
  y->z = cs->r1 + x->z * (cs->r2 - cs->r1);       // radial coordinate
}

static void spherical_jacobian_0(void* context, point_t* x, real_t* J)
{
  cs_t* cs = context;
  real_t X = -1.0 + 2.0 * x->x;
  real_t Y = -1.0 + 2.0 * x->y;
  J[0] = 1.0 / (1.0 + X*X);
  real_t W = sqrt(1.0 + X*X);
  real_t dWdX = 0.5 * pow(1.0 + X*X, -1.5) * 2.0 * X;
  real_t U = Y / sqrt(1.0 + X*X);
  real_t dUdX = (W - Y * dWdX) / (W*W);
  real_t dUdY = 1.0 / sqrt(1.0 + X*X);
  J[3] = 1.0 / (1.0 + U*U) * dUdX;
  J[4] = 1.0 / (1.0 + U*U) * dUdY;
  J[8] = bbox->r2 - bbox->r1;
  J[1] = J[2] = J[5] = J[6] = J[7] = 0.0; // other off-diagonals are zero.
}

// Inverse mappings.
static void spherical_inv_map_point_0(void* context, point_t* y, point_t* x)
{
  cs_t* cs = context;
  real_t X = tan(y->x - 0.5 * M_PI * (cs->index-1));
  real_t Y = tan(y->y) * sqrt(1.0 + X*X);
  x->x = 0.5 * (X + 1.0);
  x->y = 0.5 * (Y + 1.0);
  x->z = (y->z - cs->r1) / (cs->r2 - cs->r1);
}

static void spherical_inv_jacobian_0(void* context, point_t* x, real_t* J)
{
  cs_t* cs = context;
  real_t U = y->x - 0.5 * M_PI * (cs->index-1);
  real_t sec_U = sec(U);
  real_t tan_U = tan(U);
  J[0] = sec_U * sec_U;
  // FIXME
  J[4] = 1.0 / (bbox->y2 - bbox->y1);
  J[8] = 1.0 / (bbox->z2 - bbox->z1);
  J[1] = J[2] = J[3] = J[5] = J[6] = J[7] = 0.0; // off-diagonals.
}

static void spherical_inv_inverse_0(void* context)
{
  cs_t* cs = context;
  return spherical_cubed_sphere_coord_mapping_new(cs->r1, cs->r2, cs->index);
}

static coord_mapping_t* spherical_inverse_0(void* context)
{
  cs_t* other_cs = context;

  const char* sector[6] = {"Equator 1", "Equator 2", 
                           "Equator 3", "Equator 4", 
                           "South pole", "North pole"};
  char name[1024];
  snprintf(name, 1023, "Spherical cubed sphere inverse mapping (r1 = %g, r2 = %g, %s)",
           r1, r2, other_cs->index);

  cs_t* other_cs = polymec_malloc(sizeof(cs_t));
  cs->r1 = other_cs->r1;
  cs->r2 = other_cs->r2;
  cs->index = other_cs->index;
  coord_mapping_vtable vtable = {.map_point = spherical_inv_map_point_0,
                                 .map_vector = spherical_inv_map_vector_0,
                                 .jacobian = spherical_inv_jacobian_0,
                                 .inverse = spherical_inv_inverse_0,
                                 .dtor = polymec_free};
  return coord_mapping_new(name, cs, vtable);
}

coord_mapping_t* spherical_cubed_sphere_coord_mapping_new(real_t r1,
                                                          real_t r2,
                                                          int index)
{
  ASSERT(r1 > 0.0);
  ASSERT(r2 > r1);
  ASSERT(index >= 0);
  ASSERT(index < 6);

  const char* sector[6] = {"Equator 1", "Equator 2", 
                           "Equator 3", "Equator 4", 
                           "South pole", "North pole"};
  char name[1024];
  snprintf(name, 1023, "Spherical cubed sphere mapping (r1 = %g, r2 = %g, %s)",
           r1, r2, sector[index]);

  cs_t* cs = polymec_malloc(sizeof(cs_t));
  cs->r1 = r1;
  cs->r2 = r2;
  cs->index = index;
  coord_mapping_vtable vtable = {.dtor = polymec_free};
  switch(index)
  {
    case 0: vtable.map_point = spherical_map_point_0;
            vtable.map_vector = spherical_map_vector_0;
            vtable.jacobian = spherical_jacobian_0;
            vtable.inverse = spherical_inverse_0;
            break;
    case 1: vtable.map_point = spherical_map_point_1;
            vtable.map_vector = spherical_map_vector_1;
            vtable.jacobian = spherical_jacobian_1;
            vtable.inverse = spherical_inverse_1;
            break;
    case 2: vtable.map_point = spherical_map_point_2;
            vtable.map_vector = spherical_map_vector_2;
            vtable.jacobian = spherical_jacobian_2;
            vtable.inverse = spherical_inverse_2;
            break;
    case 3: vtable.map_point = spherical_map_point_3;
            vtable.map_vector = spherical_map_vector_3;
            vtable.jacobian = spherical_jacobian_3;
            vtable.inverse = spherical_inverse_3;
            break;
    case 4: vtable.map_point = spherical_map_point_4;
            vtable.map_vector = spherical_map_vector_4;
            vtable.jacobian = spherical_jacobian_4;
            vtable.inverse = spherical_inverse_4;
            break;
    case 5: vtable.map_point = spherical_map_point_5;
            vtable.map_vector = spherical_map_vector_5;
            vtable.jacobian = spherical_jacobian_5;
            vtable.inverse = spherical_inverse_5;
    default:
  }
  return coord_mapping_new(name, cs, vtable);
}
