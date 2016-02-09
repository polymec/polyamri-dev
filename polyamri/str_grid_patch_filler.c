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
  for (int j = dest_patch->j1; j < dest_patch->j2; ++j)
    for (int k = dest_patch->k1; k < dest_patch->k2; ++k)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[dest_patch->i2+g][j][k][c] = src_data[dest_patch->i1+g][j][k][c];

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
  for (int j = dest_patch->j1; j < dest_patch->j2; ++j)
    for (int k = dest_patch->k1; k < dest_patch->k2; ++k)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[dest_patch->i1-ng+g][j][k][c] = src_data[dest_patch->i2-ng+g][j][k][c];

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
  for (int i = dest_patch->i1; i < dest_patch->i2; ++i)
    for (int k = dest_patch->k1; k < dest_patch->k2; ++k)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[i][dest_patch->j2+g][k][c] = src_data[i][dest_patch->j1+g][k][c];

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
  for (int i = dest_patch->i1; i < dest_patch->i2; ++i)
    for (int k = dest_patch->k1; k < dest_patch->k2; ++k)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[i][dest_patch->j1-ng+g][k][c] = src_data[i][dest_patch->j2-ng+g][k][c];

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
  for (int i = dest_patch->i1; i < dest_patch->i2; ++i)
    for (int j = dest_patch->j1; j < dest_patch->j2; ++j)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[i][j][dest_patch->k2+g][c] = src_data[i][j][dest_patch->k1+g][c];

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
  for (int i = dest_patch->i1; i < dest_patch->i2; ++i)
    for (int j = dest_patch->j1; j < dest_patch->j2; ++j)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          dest_data[i][j][dest_patch->k1-ng+g][c] = src_data[i][j][dest_patch->k2-ng+g][c];

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
//                      Zero-flux patch filler
//------------------------------------------------------------------------

static int zf_x1(void* context, 
                 int i, int j, int k,
                 str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = patch->nc, ng = patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int j = patch->j1; j < patch->j2; ++j)
    for (int k = patch->k1; k < patch->k2; ++k)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          data[patch->i1-g][j][k][c] = data[patch->i1+ng-g][j][k][c];

  return -1;
}

static int zf_x2(void* context, 
                 int i, int j, int k,
                 str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = patch->nc, ng = patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int j = patch->j1; j < patch->j2; ++j)
    for (int k = patch->k1; k < patch->k2; ++k)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          data[patch->i2+g][j][k][c] = data[patch->i2-ng+g][j][k][c];

  return -1;
}

static int zf_y1(void* context, 
                 int i, int j, int k,
                 str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = patch->nc, ng = patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int i = patch->i1; i < patch->i2; ++i)
    for (int k = patch->k1; k < patch->k2; ++k)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          data[i][patch->j1-g][k][c] = data[i][patch->j1+ng-g][k][c];

  return -1;
}

static int zf_y2(void* context, 
                 int i, int j, int k,
                 str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = patch->nc, ng = patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int i = patch->i1; i < patch->i2; ++i)
    for (int k = patch->k1; k < patch->k2; ++k)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          data[i][patch->j2+g][k][c] = data[i][patch->j2-ng+g][k][c];

  return -1;
}

static int zf_z1(void* context, 
                 int i, int j, int k,
                 str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = patch->nc, ng = patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int i = patch->i1; i < patch->i2; ++i)
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          data[i][j][patch->k2-g][c] = data[i][j][patch->k2+ng-g][c];

  return -1;
}

static int zf_z2(void* context, 
                 int i, int j, int k,
                 str_grid_cell_data_t* cell_data)
{
  str_grid_patch_t* patch = str_grid_cell_data_patch(cell_data, i, j, k);
  int nc = patch->nc, ng = patch->ng;

  DECLARE_STR_GRID_PATCH_ARRAY(data, patch);
  for (int i = patch->i1; i < patch->i2; ++i)
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          data[i][j][patch->k2+g][c] = data[i][j][patch->k2-ng+g][c];

  return -1;
}

str_grid_patch_filler_t* zero_flux_str_grid_patch_filler_new(str_grid_patch_boundary_t patch_boundary)
{
  str_grid_patch_filler_vtable vtable;
  switch(patch_boundary)
  {
    case STR_GRID_PATCH_X1_BOUNDARY: vtable.start_filling_cells = zf_x1; break;
    case STR_GRID_PATCH_X2_BOUNDARY: vtable.start_filling_cells = zf_x2; break;
    case STR_GRID_PATCH_Y1_BOUNDARY: vtable.start_filling_cells = zf_y1; break;
    case STR_GRID_PATCH_Y2_BOUNDARY: vtable.start_filling_cells = zf_y2; break;
    case STR_GRID_PATCH_Z1_BOUNDARY: vtable.start_filling_cells = zf_z1; break;
    case STR_GRID_PATCH_Z2_BOUNDARY: vtable.start_filling_cells = zf_z2; break;
  }
  return str_grid_patch_filler_new("zero-flux patch filler", NULL, vtable);
}

