// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "polyamri/amr_grid_assembly.h"

struct amr_grid_assembly_t 
{
  string_ptr_unordered_map_t* blocks;
};

amr_grid_assembly_t* amr_grid_assembly_new()
{
  amr_grid_assembly_t* assembly = polymec_malloc(sizeof(amr_grid_assembly_t));
  assembly->blocks = string_ptr_unordered_map_new();
  return assembly;
}

void amr_grid_assembly_free(amr_grid_assembly_t* assembly)
{
  string_ptr_unordered_map_free(assembly->blocks);
  polymec_free(assembly);
}

int amr_grid_assembly_num_blocks(amr_grid_assembly_t* assembly)
{
  return assembly->blocks->size;
}

void amr_grid_assembly_add_block(amr_grid_assembly_t* assembly,
                                 const char* block_name,       
                                 amr_grid_hierarchy_t* block)
{
  string_ptr_unordered_map_insert_with_kv_dtors(assembly->blocks, 
                                                string_dup(block_name),
                                                block,
                                                string_free,
                                                DTOR(amr_grid_hierarchy_free));
}

amr_grid_hierarchy_t* amr_grid_assembly_block(amr_grid_assembly_t* assembly,
                                              const char* block_name)
{
  amr_grid_hierarchy_t** block = (amr_grid_hierarchy_t**)
                                 string_ptr_unordered_map_get(assembly->blocks, (char*)block_name);
  if (block == NULL)
    return NULL;
  else
    return *block;
}

bool amr_grid_assembly_next_block(amr_grid_assembly_t* assembly, 
                                  int* pos, 
                                  char** block_name,
                                  amr_grid_hierarchy_t** block)
{
  return string_ptr_unordered_map_next(assembly->blocks, pos, 
                                       block_name,
                                       (void**)block);
}

