// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "polyamri/str_grid_assembly.h"

typedef struct
{
  str_grid_t* block;
  str_grid_t* neighbors[6];
} block_entry_t;

static void block_entry_free(block_entry_t* entry)
{
  str_grid_free(entry->block);
  polymec_free(entry);
}

struct str_grid_assembly_t 
{
  // Blocks in the assembly.
  string_ptr_unordered_map_t* blocks;

  // Metadata for each block.
  int npx, npy, npz, nx, ny, nz;
};

str_grid_assembly_t* str_grid_assembly_new()
{
  str_grid_assembly_t* assembly = polymec_malloc(sizeof(str_grid_assembly_t));
  assembly->blocks = string_ptr_unordered_map_new();
  return assembly;
}

void str_grid_assembly_free(str_grid_assembly_t* assembly)
{
  string_ptr_unordered_map_free(assembly->blocks);
  polymec_free(assembly);
}

int str_grid_assembly_num_blocks(str_grid_assembly_t* assembly)
{
  return assembly->blocks->size;
}

void str_grid_assembly_add_block(str_grid_assembly_t* assembly,
                                 const char* block_name,       
                                 str_grid_t* block)
{
  // If this is the first block, copy its metadata. 
  if (string_ptr_unordered_map_empty(assembly->blocks))
  {
    str_grid_get_extents(block, &assembly->npx, &assembly->npy, &assembly->npz);
    str_grid_get_patch_size(block, &assembly->nx, &assembly->ny, &assembly->nz);
  }
  // Otherwise, validate this block against existing metadata.
  else
  {
    int npx, npy, npz, nx, ny, nz;
    str_grid_get_extents(block, &npx, &npy, &npz);
    str_grid_get_patch_size(block, &nx, &ny, &nz);
    if (npx != assembly->npx)
    {
      polymec_error("str_grid_assembly_add_block: Number of patches in x "
                    "direction (%d) of block does not match that of other "
                    "blocks in the assembly (%d).", npx, assembly->npx);
    }
    if (npy != assembly->npy)
    {
      polymec_error("str_grid_assembly_add_block: Number of patches in y "
                    "direction (%d) of block does not match that of other "
                    "blocks in the assembly (%d).", npy, assembly->npy);
    }
    if (npz != assembly->npz)
    {
      polymec_error("str_grid_assembly_add_block: Number of patches in z "
                    "direction (%d) of block does not match that of other "
                    "blocks in the assembly (%d).", npz, assembly->npz);
    }
    if (nx != assembly->nx)
    {
      polymec_error("str_grid_assembly_add_block: Patch x size of (%d) of "
                    "block does not match that of other blocks in the "
                    "assembly (%d).", nx, assembly->nx);
    }
    if (ny != assembly->ny)
    {
      polymec_error("str_grid_assembly_add_block: Patch y size of (%d) of "
                    "block does not match that of other blocks in the "
                    "assembly (%d).", ny, assembly->ny);
    }
    if (nz != assembly->nz)
    {
      polymec_error("str_grid_assembly_add_block: Patch z size of (%d) of "
                    "block does not match that of other blocks in the "
                    "assembly (%d).", nz, assembly->nz);
    }
  }

  block_entry_t* entry = polymec_malloc(sizeof(block_entry_t));
  entry->block = block;
  memset(entry->neighbors, 0, sizeof(str_grid_t*) * 6);
  string_ptr_unordered_map_insert_with_kv_dtors(assembly->blocks, 
                                                string_dup(block_name), 
                                                entry,
                                                string_free,
                                                DTOR(block_entry_free));
}

void str_grid_assembly_connect(str_grid_assembly_t* assembly,
                               const char* block1_name,
                               str_grid_boundary_t block1_boundary,
                               const char* block2_name,
                               str_grid_boundary_t block2_boundary)
{
  block_entry_t** entry1_p = (block_entry_t**)string_ptr_unordered_map_get(assembly->blocks, (char*)block1_name);
  ASSERT(entry1_p != NULL);
  block_entry_t** entry2_p = (block_entry_t**)string_ptr_unordered_map_get(assembly->blocks, (char*)block2_name);
  ASSERT(entry2_p != NULL);
  block_entry_t* entry1 = *entry1_p;
  block_entry_t* entry2 = *entry2_p;
  int i1 = block1_boundary;
  int i2 = block2_boundary;
  entry1->neighbors[i1] = entry2->block;
  entry2->neighbors[i2] = entry1->block;
}

str_grid_t* str_grid_assembly_block(str_grid_assembly_t* assembly,
                                    const char* block_name)
{
  block_entry_t** entry_p = (block_entry_t**)string_ptr_unordered_map_get(assembly->blocks, (char*)block_name);
  if (entry_p != NULL)
    return (*entry_p)->block;
  else
    return NULL;
}

bool str_grid_assembly_next_block(str_grid_assembly_t* assembly, 
                                  int* pos, 
                                  char** block_name,
                                  str_grid_t** block,
                                  str_grid_t** block_x1_neighbor,
                                  str_grid_t** block_x2_neighbor,
                                  str_grid_t** block_y1_neighbor,
                                  str_grid_t** block_y2_neighbor,
                                  str_grid_t** block_z1_neighbor,
                                  str_grid_t** block_z2_neighbor)
{
  block_entry_t* entry;
  char* bname;
  bool result = string_ptr_unordered_map_next(assembly->blocks, pos, 
                                              &bname, (void**)&entry);
  if (result)
  {
    if (block_name != NULL)
      *block_name = bname;
    *block = entry->block;
    if (block_x1_neighbor != NULL)
      *block_x1_neighbor = entry->neighbors[STR_GRID_X1_BOUNDARY];
    if (block_x2_neighbor != NULL)
      *block_x2_neighbor = entry->neighbors[STR_GRID_X2_BOUNDARY];
    if (block_y1_neighbor != NULL)
      *block_y1_neighbor = entry->neighbors[STR_GRID_Y1_BOUNDARY];
    if (block_y2_neighbor != NULL)
      *block_y2_neighbor = entry->neighbors[STR_GRID_Y2_BOUNDARY];
    if (block_z1_neighbor != NULL)
      *block_z1_neighbor = entry->neighbors[STR_GRID_Z1_BOUNDARY];
    if (block_z2_neighbor != NULL)
      *block_z2_neighbor = entry->neighbors[STR_GRID_Z2_BOUNDARY];
  }
  return result;
}

bool str_grid_assembly_next_patch(str_grid_assembly_t* assembly,
                                  int* pos, 
                                  str_grid_t** block, 
                                  int* i, int* j, int* k)
{
  int num_grid_patches = assembly->npx * assembly->npy * assembly->npz;
  int block_pos = *pos / num_grid_patches;
  int patch_pos = *pos % num_grid_patches;
  ++(*pos);

  bool result = str_grid_assembly_next_block(assembly, &block_pos, NULL, block,
                                             NULL, NULL, NULL, NULL, NULL, NULL);
  if (result)
    return str_grid_next_patch(*block, &patch_pos, i, j, k);
  else
    return false;
}
