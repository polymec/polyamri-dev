// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/silo_file_write_amr_data_hierarchy.h"
#include "silo.h"
#include "pmpio.h"

// These functions are implemented in polymec/core/silo_file.c, and used 
// here, even though they are not part of polymec's API.

// Here's the silo_file_t struct, replicated here from the source in 
// polymec/core/silo_file.c. Note that whenever that source is changed, 
// we must update it here as well. This logic is tightly coupled to that 
// in polymec, anyway, so this is a calculated cost.
struct silo_file_t 
{
  // File data.
  DBfile* dbfile;

  // Metadata.
  char prefix[FILENAME_MAX], directory[FILENAME_MAX], filename[FILENAME_MAX];
  int cycle;
  real_t time;
  int mode; // Open for reading (DB_READ) or writing (DB_CLOBBER)? 
  string_ptr_unordered_map_t* expressions;

#if POLYMEC_HAVE_MPI
  // Stuff for poor man's parallel I/O.
  PMPIO_baton_t* baton;
  MPI_Comm comm;
  int num_files, mpi_tag, nproc, rank, group_rank, rank_in_group;
  ptr_array_t* multimeshes;
  ptr_array_t* multivars;
#endif
};

void silo_file_write_amr_data_hierarchy(silo_file_t* file, 
                                        const char* hierarchy_name,
                                        amr_data_hierarchy_t* hierarchy,
                                        silo_field_metadata_t** field_metadata)
{
  // Create a new MRG tree that can hold an AMR hierarchy.
  int num_levels = amr_data_hierarchy_num_levels(hierarchy);

  // The following stuff is taken largely from the add_amr_mrgtree.c example
  // provided by the Silo folks.
  DBmrgtree* tree = DBMakeMrgtree(DB_MULTIMESH, 0, 2, NULL);
  DBAddRegion(tree, "amr_decomp", 0, 2, NULL, 0, NULL, NULL, NULL, NULL); 
  DBSetCwr(tree, "amr_decomp");
  DBAddRegion(tree, "levels", 0, num_levels, NULL, 0, NULL, NULL, NULL, NULL); 
  DBSetCwr(tree, "levels");

  // Define each AMR refinement level in the tree.
  int num_patches = 0;
  {
    const char* level_maps_name = "mesh_level_maps";
    char* level_region_names[1];
    int* seg_types = polymec_malloc(sizeof(int) * num_levels);
    int* seg_ids = polymec_malloc(sizeof(int) * num_levels);
    int* num_level_patches = polymec_malloc(sizeof(int) * num_levels);
    int pos = 0, i = 0;
    amr_grid_data_t* grid_data;
    while (amr_data_hierarchy_next_coarsest(hierarchy, &pos, &grid_data))
    {
      seg_ids[i] = i;
      seg_types[i] = DB_BLOCKCENT;
      num_level_patches[i] = amr_grid_data_num_patches(grid_data);
      num_patches += num_level_patches[i];
      ++i;
    }
    level_region_names[0] = "@level%d@n";
    DBAddRegionArray(tree, num_levels, (const char* const*)level_region_names, 
                     0, level_maps_name, 1, seg_ids, num_level_patches, seg_types, NULL);
  }
  DBSetCwr(tree, "..");

  // Now define the patches.
  DBAddRegion(tree, "patches", 0, num_patches, NULL, 0, NULL, NULL, NULL, NULL); 
  DBSetCwr(tree, "patches");

  // FIXME: Fancy things go here.

  // Write some metadata.
  char *mrgv_onames[5], temp_name[FILENAME_MAX];
  const char* mesh_name = "mesh";
	sprintf(temp_name, "%s_wmrgtree_lvlRatios", mesh_name);
  mrgv_onames[0] = string_dup(temp_name);
	sprintf(temp_name, "%s_wmrgtree_ijkExts", mesh_name);
  mrgv_onames[1] = string_dup(temp_name);
	sprintf(temp_name, "%s_wmrgtree_xyzExts", mesh_name);
  mrgv_onames[2] = string_dup(temp_name);
  mrgv_onames[3] = string_dup("rank");
  mrgv_onames[4] = NULL;
  DBoptlist* options = DBMakeOptlist(10);
  DBAddOption(options, DBOPT_MRGV_ONAMES, mrgv_onames);

  // Deposit the MRG tree into the Silo file.
  DBPutMrgtree(file->dbfile, "mrgTree", "amr_mesh", tree, options);
  DBFreeMrgtree(tree);
}

