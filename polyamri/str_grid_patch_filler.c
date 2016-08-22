// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "core/unordered_map.h"
#include "core/tuple.h"
#include "polyamri/str_grid_patch_filler.h"

struct str_grid_patch_filler_t 
{
  char* name;
  void* context;
  str_grid_patch_filler_vtable vtable;
};

static void str_grid_patch_filler_free(void* ctx, void* dummy)
{
  str_grid_patch_filler_t* filler = ctx;
  string_free(filler->name);
  if ((filler->vtable.dtor != NULL) && (filler->context != NULL))
    filler->vtable.dtor(filler->context);
}

str_grid_patch_filler_t* str_grid_patch_filler_new(const char* name,
                                                   void* context,
                                                   str_grid_patch_filler_vtable vtable)
{
  ASSERT(vtable.start_filling_cells != NULL);
  str_grid_patch_filler_t* filler = GC_MALLOC(sizeof(str_grid_patch_filler_t));
  filler->name = string_dup(name);
  filler->context = context;
  filler->vtable = vtable;
  GC_register_finalizer(filler, str_grid_patch_filler_free, filler, NULL, NULL);
  return filler;
}

char* str_grid_patch_filler_name(str_grid_patch_filler_t* filler)
{
  return filler->name;
}

int str_grid_patch_filler_start(str_grid_patch_filler_t* filler, 
                                int i, int j, int k,
                                str_grid_cell_data_t* cell_data)
{
  return filler->vtable.start_filling_cells(filler->context, i, j, k, cell_data);
}

void str_grid_patch_filler_finish(str_grid_patch_filler_t* filler,
                                  int token)
{
  ASSERT(token >= 0);
  if (filler->vtable.finish_filling_cells != NULL)
    filler->vtable.finish_filling_cells(filler->context, token);
}

//------------------------------------------------------------------------
//                  Local copying patch ghost filler
//------------------------------------------------------------------------

static int copy_x1_to_x2(void* context, 
                         int i, int j, int k, 
                         str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* src_patch = str_grid_cell_data_patch(cell_data, i+1, j, k);
  str_grid_patch_t* dest_patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = dest_patch->nc;
  int ng = dest_patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(src_data, src_patch);
  DECLARE_STR_GRID_PATCH_ARRAY(dest_data, dest_patch);
  for (int jj = dest_patch->j1; jj < dest_patch->j2; ++jj)
    for (int kk = dest_patch->k1; kk < dest_patch->k2; ++kk)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[dest_patch->i2+g][jj][kk][c] = src_data[dest_patch->i1+g][jj][kk][c];

  return -1;
}

static int copy_x2_to_x1(void* context, 
                         int i, int j, int k, 
                         str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* src_patch = str_grid_cell_data_patch(cell_data, i-1, j, k);
  str_grid_patch_t* dest_patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = dest_patch->nc;
  int ng = dest_patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(src_data, src_patch);
  DECLARE_STR_GRID_PATCH_ARRAY(dest_data, dest_patch);
  for (int jj = dest_patch->j1; jj < dest_patch->j2; ++jj)
    for (int kk = dest_patch->k1; kk < dest_patch->k2; ++kk)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[dest_patch->i1-ng+g][jj][kk][c] = src_data[dest_patch->i2-ng+g][jj][kk][c];

  return -1;
}

static int copy_y1_to_y2(void* context, 
                         int i, int j, int k, 
                         str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* src_patch = str_grid_cell_data_patch(cell_data, i, j+1, k);
  str_grid_patch_t* dest_patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = dest_patch->nc;
  int ng = dest_patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(src_data, src_patch);
  DECLARE_STR_GRID_PATCH_ARRAY(dest_data, dest_patch);
  for (int ii = dest_patch->i1; ii < dest_patch->i2; ++ii)
    for (int kk = dest_patch->k1; kk < dest_patch->k2; ++kk)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[ii][dest_patch->j2+g][kk][c] = src_data[ii][dest_patch->j1+g][kk][c];

  return -1;
}

static int copy_y2_to_y1(void* context, 
                         int i, int j, int k, 
                         str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* src_patch = str_grid_cell_data_patch(cell_data, i, j-1, k);
  str_grid_patch_t* dest_patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = dest_patch->nc;
  int ng = dest_patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(src_data, src_patch);
  DECLARE_STR_GRID_PATCH_ARRAY(dest_data, dest_patch);
  for (int ii = dest_patch->i1; ii < dest_patch->i2; ++ii)
    for (int kk = dest_patch->k1; kk < dest_patch->k2; ++kk)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[ii][dest_patch->j1-ng+g][kk][c] = src_data[ii][dest_patch->j2-ng+g][kk][c];

  return -1;
}

static int copy_z1_to_z2(void* context, 
                         int i, int j, int k, 
                         str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* src_patch = str_grid_cell_data_patch(cell_data, i, j, k+1);
  str_grid_patch_t* dest_patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = dest_patch->nc;
  int ng = dest_patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(src_data, src_patch);
  DECLARE_STR_GRID_PATCH_ARRAY(dest_data, dest_patch);
  for (int ii = dest_patch->i1; ii < dest_patch->i2; ++ii)
    for (int jj = dest_patch->j1; jj < dest_patch->j2; ++jj)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[ii][jj][dest_patch->k2+g][c] = src_data[ii][jj][dest_patch->k1+g][c];

  return -1;
}

static int copy_z2_to_z1(void* context, 
                         int i, int j, int k, 
                         str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* src_patch = str_grid_cell_data_patch(cell_data, i, j, k-1);
  str_grid_patch_t* dest_patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = dest_patch->nc;
  int ng = dest_patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(src_data, src_patch);
  DECLARE_STR_GRID_PATCH_ARRAY(dest_data, dest_patch);
  for (int ii = dest_patch->i1; ii < dest_patch->i2; ++ii)
    for (int jj = dest_patch->j1; jj < dest_patch->j2; ++jj)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[ii][jj][dest_patch->k1-ng+g][c] = src_data[ii][jj][dest_patch->k2-ng+g][c];

  return -1;
}

str_grid_patch_filler_t* copy_str_grid_patch_filler_new(str_grid_patch_boundary_t src_boundary,
                                                        str_grid_patch_boundary_t dest_boundary)
{
  ASSERT(src_boundary != dest_boundary);

  char name[1025];
  str_grid_patch_filler_vtable vtable;
  if (src_boundary == STR_GRID_PATCH_X1_BOUNDARY)
  {
    ASSERT(dest_boundary == STR_GRID_PATCH_X2_BOUNDARY);
    strcpy(name, "copy(x1 -> x2)");
    vtable.start_filling_cells = copy_x1_to_x2;
  }
  else if (src_boundary == STR_GRID_PATCH_X2_BOUNDARY)
  {
    ASSERT(dest_boundary == STR_GRID_PATCH_X1_BOUNDARY);
    strcpy(name, "copy(x2 -> x1)");
    vtable.start_filling_cells = copy_x2_to_x1;
  }
  else if (src_boundary == STR_GRID_PATCH_Y1_BOUNDARY)
  {
    ASSERT(dest_boundary == STR_GRID_PATCH_Y2_BOUNDARY);
    strcpy(name, "copy(y1 -> y2)");
    vtable.start_filling_cells = copy_y1_to_y2;
  }
  else if (src_boundary == STR_GRID_PATCH_Y2_BOUNDARY)
  {
    ASSERT(dest_boundary == STR_GRID_PATCH_Y1_BOUNDARY);
    strcpy(name, "copy(y2 -> y1)");
    vtable.start_filling_cells = copy_y2_to_y1;
  }
  else if (src_boundary == STR_GRID_PATCH_Z1_BOUNDARY)
  {
    ASSERT(dest_boundary == STR_GRID_PATCH_Z2_BOUNDARY);
    strcpy(name, "copy(z1 -> z2)");
    vtable.start_filling_cells = copy_z1_to_z2;
  }
  else if (src_boundary == STR_GRID_PATCH_Z2_BOUNDARY)
  {
    ASSERT(dest_boundary == STR_GRID_PATCH_Z1_BOUNDARY);
    strcpy(name, "copy(z2 -> z1)");
    vtable.start_filling_cells = copy_z2_to_z1;
  }

  return str_grid_patch_filler_new(name, NULL, vtable);
}

//------------------------------------------------------------------------
//                      Robin BC patch filler
//------------------------------------------------------------------------

typedef struct
{
  real_t A;
  real_t B;
  real_t C;
  real_t h;
  int component;
} robin_bc_t;

static int robin_x1(void* context, 
                    int i, int j, int k,
                    str_grid_cell_data_t* cell_data)
{
  robin_bc_t* robin = context;
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int ng = patch->ng, c = robin->component;
  ASSERT(c < patch->nc);
  real_t num = 0.5 * robin->A - robin->B/robin->h;
  real_t denom = 0.5 * robin->A + robin->B/robin->h;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int jj = patch->j1; jj < patch->j2; ++jj)
    for (int kk = patch->k1; kk < patch->k2; ++kk)
      for (int g = 0; g < ng; ++g)
        data[patch->i1-g][jj][kk][c] = (robin->C - num * data[patch->i1+ng-g][jj][kk][c]) / denom;

  return -1;
}

static int robin_x2(void* context, 
                    int i, int j, int k,
                    str_grid_cell_data_t* cell_data)
{
  robin_bc_t* robin = context;
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int ng = patch->ng, c = robin->component;
  ASSERT(c < patch->nc);
  real_t num = 0.5 * robin->A - robin->B/robin->h;
  real_t denom = 0.5 * robin->A + robin->B/robin->h;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int jj = patch->j1; jj < patch->j2; ++jj)
    for (int kk = patch->k1; kk < patch->k2; ++kk)
      for (int g = 0; g < ng; ++g)
        data[patch->i2+g][jj][kk][c] = (robin->C - num * data[patch->i2-ng+g][jj][kk][c]) / denom;

  return -1;
}

static int robin_y1(void* context, 
                    int i, int j, int k,
                    str_grid_cell_data_t* cell_data)
{
  robin_bc_t* robin = context;
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int ng = patch->ng, c = robin->component;
  ASSERT(c < patch->nc);
  real_t num = 0.5 * robin->A - robin->B/robin->h;
  real_t denom = 0.5 * robin->A + robin->B/robin->h;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int ii = patch->i1; ii < patch->i2; ++ii)
    for (int kk = patch->k1; kk < patch->k2; ++kk)
      for (int g = 0; g < ng; ++g)
        data[ii][patch->j1-g][kk][c] = (robin->C - num * data[ii][patch->j1+ng-g][kk][c]) / denom;

  return -1;
}

static int robin_y2(void* context, 
                    int i, int j, int k,
                    str_grid_cell_data_t* cell_data)
{
  robin_bc_t* robin = context;
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int ng = patch->ng, c = robin->component;
  ASSERT(c < patch->nc);
  real_t num = 0.5 * robin->A - robin->B/robin->h;
  real_t denom = 0.5 * robin->A + robin->B/robin->h;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int ii = patch->i1; ii < patch->i2; ++ii)
    for (int kk = patch->k1; kk < patch->k2; ++kk)
      for (int g = 0; g < ng; ++g)
        data[ii][patch->j2+g][kk][c] = (robin->C - num * data[ii][patch->j2-ng+g][kk][c]) / denom;

  return -1;
}

static int robin_z1(void* context, 
                    int i, int j, int k,
                    str_grid_cell_data_t* cell_data)
{
  robin_bc_t* robin = context;
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int ng = patch->ng, c = robin->component;
  ASSERT(c < patch->nc);
  real_t num = 0.5 * robin->A - robin->B/robin->h;
  real_t denom = 0.5 * robin->A + robin->B/robin->h;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int ii = patch->i1; ii < patch->i2; ++ii)
    for (int jj = patch->j1; jj < patch->j2; ++jj)
      for (int g = 0; g < ng; ++g)
        data[ii][jj][patch->k2-g][c] = (robin->C - num * data[ii][jj][patch->k2+ng-g][c]) / denom;

  return -1;
}

static int robin_z2(void* context, 
                    int i, int j, int k,
                    str_grid_cell_data_t* cell_data)
{
  robin_bc_t* robin = context;
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int ng = patch->ng, c = robin->component;
  ASSERT(c < patch->nc);
  real_t num = 0.5 * robin->A - robin->B/robin->h;
  real_t denom = 0.5 * robin->A + robin->B/robin->h;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int ii = patch->i1; ii < patch->i2; ++ii)
    for (int jj = patch->j1; jj < patch->j2; ++jj)
      for (int g = 0; g < ng; ++g)
        data[ii][jj][patch->k2+g][c] = (robin->C - num * data[ii][jj][patch->k2-ng+g][c]) / denom;

  return -1;
}

static robin_bc_t* robin_bc_new(real_t A, real_t B, real_t C, real_t h, int component)
{
  robin_bc_t* robin = polymec_malloc(sizeof(robin_bc_t));
  robin->A = A;
  robin->B = B;
  robin->C = C;
  robin->h = h;
  robin->component = component;
  return robin;
}

str_grid_patch_filler_t* robin_bc_str_grid_patch_filler_new(real_t A, 
                                                            real_t B, 
                                                            real_t C,
                                                            real_t h,
                                                            int component,
                                                            str_grid_patch_boundary_t patch_boundary)
{
  ASSERT(reals_equal(A, 0.0) || reals_equal(B, 0.0) || reals_equal(C, 0.0));
  ASSERT(h > 0.0);
  ASSERT(component >= 0);

  robin_bc_t* robin = robin_bc_new(A, B, C, h, component);
  str_grid_patch_filler_vtable vtable = {.dtor = polymec_free};
  switch(patch_boundary)
  {
    case STR_GRID_PATCH_X1_BOUNDARY: vtable.start_filling_cells = robin_x1; break;
    case STR_GRID_PATCH_X2_BOUNDARY: vtable.start_filling_cells = robin_x2; break;
    case STR_GRID_PATCH_Y1_BOUNDARY: vtable.start_filling_cells = robin_y1; break;
    case STR_GRID_PATCH_Y2_BOUNDARY: vtable.start_filling_cells = robin_y2; break;
    case STR_GRID_PATCH_Z1_BOUNDARY: vtable.start_filling_cells = robin_z1; break;
    case STR_GRID_PATCH_Z2_BOUNDARY: vtable.start_filling_cells = robin_z2; break;
  }
  char name[1025];
  const char* boundaries[] = {"x1", "x2", "y1", "y2", "z1", "z2"};
  if (reals_equal(A, 0.0))
  {
    snprintf(name, 1024, "Neumann BC (%g * dU/dn = %g, h = %g, component %d) on %s boundary", 
             robin->B, robin->C, robin->h, robin->component, boundaries[patch_boundary]);
  }
  else if (reals_equal(robin->B, 0.0))
  {
    snprintf(name, 1024, "Dirichlet BC (%g * U = %g, component %d) on %s boundary", 
             robin->A, robin->C, robin->component, boundaries[patch_boundary]);
  }
  else
  {
    snprintf(name, 1024, "Robin BC (%g * U + %g * dU/dn = %g, h = %g, component %d) on %s boundary", 
             robin->A, robin->B, robin->C, robin->h, robin->component, boundaries[patch_boundary]);
  }
  return str_grid_patch_filler_new(name, robin, vtable);
}

str_grid_patch_filler_t* dirichlet_bc_str_grid_patch_filler_new(real_t F, 
                                                                int component,
                                                                str_grid_patch_boundary_t patch_boundary)
{
  return robin_bc_str_grid_patch_filler_new(1.0, 0.0, F, 1.0, component, patch_boundary);
}

str_grid_patch_filler_t* neumann_bc_str_grid_patch_filler_new(real_t A, 
                                                              real_t B,
                                                              real_t h,
                                                              int component,
                                                              str_grid_patch_boundary_t patch_boundary)
{
  return robin_bc_str_grid_patch_filler_new(0.0, A, B, h, component, patch_boundary);
}

