// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
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
  str_grid_patch_filler_t* filler = polymec_malloc(sizeof(str_grid_patch_filler_t));
  filler->name = string_dup(name);
  filler->context = context;
  filler->vtable = vtable;
  GC_register_finalizer(filler, str_grid_patch_filler_free, filler, NULL, NULL);
  return filler;
}

int str_grid_patch_filler_start(str_grid_patch_filler_t* filler,
                                str_grid_patch_t* patch)
{
  return filler->vtable.start_filling_cells(filler->context, patch);
}

void str_grid_patch_filler_finish(str_grid_patch_filler_t* filler,
                                  int token)
{
  if (filler->vtable.finish_filling_cells != NULL)
    filler->vtable.finish_filling_cells(filler->context, token);
}

