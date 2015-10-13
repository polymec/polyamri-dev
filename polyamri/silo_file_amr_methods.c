// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "silo.h"
#include "core/silo_file.h"
#include "polyamri/silo_file_amr_methods.h"

#if POLYMEC_HAVE_DOUBLE_PRECISION
#define SILO_FLOAT_TYPE DB_DOUBLE
#else
#define SILO_FLOAT_TYPE DB_FLOAT
#endif

// These functions are implemented in polymec/core/silo_file.c, and used 
// here, even though they are not part of polymec's API.
extern DBfile* silo_file_dbfile(silo_file_t* file);
extern DBoptlist* optlist_from_metadata(silo_field_metadata_t* metadata);
extern void optlist_free(DBoptlist* optlist);
extern void silo_file_add_subdomain_mesh(silo_file_t* file, const char* mesh_name, int silo_mesh_type, DBoptlist* optlist);
extern void silo_file_add_subdomain_field(silo_file_t* file, const char* mesh_name, const char* field_name, int silo_field_type, DBoptlist* optlist);

void silo_file_write_amr_patch(silo_file_t* file, 
                               const char** field_component_names,
                               const char* patch_name,
                               amr_patch_t* patch,
                               bbox_t* bbox)
{
  // If there's no bounding box, call the mapped version.
  if (bbox == NULL)
  {
    silo_file_write_mapped_amr_patch(file, field_component_names, patch_name, patch, NULL);
    return;
  }

  DBfile* dbfile = silo_file_dbfile(file);

  // Name the coordinate axes.
  const char* const coord_names[3] = {"x", "y", "z"};

  // Provide coordinates.
  int nx = patch->i2 - patch->i1, 
      ny = patch->j2 - patch->j1, 
      nz = patch->k2 - patch->k1;
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  real_t x_node[nx+1], y_node[ny+1], z_node[nz+1];
  int dimensions[3];
  dimensions[0] = nx+1;
  dimensions[1] = ny+1;
  dimensions[2] = nz+1;
  real_t dx = (bbox->x2 - bbox->x1) / nx;
  for (int i = patch->i1; i <= patch->i2; ++i)
    x_node[i] = bbox->x1 + i*dx;
  real_t dy = (bbox->y2 - bbox->y1) / ny;
  for (int j = patch->j1; j <= patch->j2; ++j)
    y_node[j] = bbox->y1 + j*dy;
  real_t dz = (bbox->z2 - bbox->z1) / nz;
  for (int k = patch->k1; k <= patch->k2; ++k)
    z_node[k] = bbox->z1 + k*dz;
  real_t* coords[3];
  coords[0] = x_node;
  coords[1] = y_node;
  coords[2] = z_node;

  // Write the grid.
  char patch_mesh_name[FILENAME_MAX];
  snprintf(patch_mesh_name, FILENAME_MAX-1, "%s_grid", patch_name);
  DBPutQuadmesh(dbfile, patch_mesh_name, coord_names, coords, dimensions, 3, 
                SILO_FLOAT_TYPE, DB_COLLINEAR, NULL);

  // Write the data.
  real_t* data[patch->nc];
  for (int c = 0; c < patch->nc; ++c)
    data[c] = polymec_malloc(sizeof(real_t) * nx*ny*nz);
  DECLARE_AMR_PATCH_ARRAY(a, patch);
  int l = 0;
  for (int i = patch->i1; i < patch->i2; ++i)
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int k = patch->k1; k < patch->k2; ++k, ++l)
        for (int c = 0; c < patch->nc; ++c)
          data[c][l] = a[i][j][k][c];

  int cell_dims[3] = {nx, ny, nz};
  DBPutQuadvar(dbfile, patch_name, patch_mesh_name, patch->nc, 
               (const char* const *)field_component_names, data, cell_dims, 
               3, NULL, 0, SILO_FLOAT_TYPE, DB_ZONECENT, NULL);

  // Clean up.
  for (int c = 0; c < patch->nc; ++c)
    polymec_free(data[c]);
}

static void write_amr_patch_grid(silo_file_t* file,
                                 const char* patch_grid_name,
                                 int i1, int i2, int j1, int j2, int k1, int k2,
                                 sp_func_t* mapping)
{
  ASSERT(i2 > i1);
  ASSERT(j2 > j1);
  ASSERT(k2 > k1);
  ASSERT((mapping == NULL) || (sp_func_num_comp(mapping) == 3));

  DBfile* dbfile = silo_file_dbfile(file);

  // Name the coordinate axes.
  const char* const coord_names[3] = {"x1", "x2", "x3"};

  // Provide coordinates.
  int n1 = i2 - i1, n2 = j2 - j1, n3 = k2 - k1;
  ASSERT(n1 > 0);
  ASSERT(n2 > 0);
  ASSERT(n3 > 0);
  real_t x1_node[n1+1], x2_node[n2+1], x3_node[n3+1];
  int dimensions[3] = {n1+1, n2+1, n3+1};
  if (mapping != NULL)
  {
    real_t dx1 = 1.0 / n1, dx2 = 1.0 / n2, dx3 = 1.0 / n3;
    point_t x;
    int l = 0;
    for (int i = i1; i <= i2; ++i)
    {
      x.x = i * dx1;
      for (int j = j1; j <= j2; ++j)
      {
        x.y = j * dx2;
        for (int k = k1; k <= k2; ++k, ++l)
        {
          x.z = k * dx3;
          real_t y[3];
          sp_func_eval(mapping, &x, y);
          x1_node[l] = y[0];
          x2_node[l] = y[1];
          x3_node[l] = y[2];
        }
      }
    }
  }
  else
  {
    real_t dx = 1.0 / n1;
    for (int i = i1; i <= i2; ++i)
      x1_node[i] = i*dx;
    real_t dy = 1.0 / n2;
    for (int j = j1; j <= j2; ++j)
      x2_node[j] = j*dy;
    real_t dz = 1.0 / n3;
    for (int k = k1; k <= k2; ++k)
      x3_node[k] = k*dz;
  }

  real_t* coords[3];
  coords[0] = x1_node;
  coords[1] = x2_node;
  coords[2] = x3_node;

  // Write the patch grid.
  DBPutQuadmesh(dbfile, patch_grid_name, coord_names, coords, dimensions, 3, 
                SILO_FLOAT_TYPE, DB_COLLINEAR, NULL);
}

static void write_amr_patch_data(silo_file_t* file,
                                 const char** field_component_names,
                                 const char* patch_grid_name,
                                 amr_patch_t* patch,
                                 silo_field_metadata_t** field_metadata)
{
  ASSERT(patch->i2 > patch->i1);
  ASSERT(patch->j2 > patch->j1);
  ASSERT(patch->k2 > patch->k1);

  DBfile* dbfile = silo_file_dbfile(file);

  // Write the data.
  int n1 = patch->i2 - patch->i1;
  int n2 = patch->j2 - patch->j1;
  int n3 = patch->k2 - patch->k1;
  real_t* data = polymec_malloc(sizeof(real_t) * n1*n2*n3);
  int cell_dims[3] = {n1, n2, n3};
  DECLARE_AMR_PATCH_ARRAY(a, patch);
  for (int c = 0; c < patch->nc; ++c)
  {
    int l = 0;
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int j = patch->j1; j < patch->j2; ++j)
        for (int k = patch->k1; k < patch->k2; ++k, ++l)
          for (int c = 0; c < patch->nc; ++c)
            data[l] = a[i][j][k][c];
    DBoptlist* optlist = (field_metadata != NULL) ? optlist_from_metadata(field_metadata[c]) : NULL;
    DBPutQuadvar1(dbfile, field_component_names[c], patch_grid_name, patch->data,
                  cell_dims, 3, NULL, 0, SILO_FLOAT_TYPE, DB_ZONECENT, optlist);
    optlist_free(optlist);
  }

  // Clean up.
  polymec_free(data);
}

void silo_file_write_mapped_amr_patch(silo_file_t* file, 
                                      const char** field_component_names,
                                      const char* patch_grid_name,
                                      amr_patch_t* patch,
                                      sp_func_t* mapping)
{
  write_amr_patch_grid(file, patch_grid_name, patch->i1, patch->i2, 
                       patch->j1, patch->j2, patch->k1, patch->k2, mapping);
  write_amr_patch_data(file, field_component_names, patch_grid_name, patch, NULL);
}

bool silo_file_contains_amr_patch(silo_file_t* file, 
                                  const char* patch_name)
{
  DBfile* dbfile = silo_file_dbfile(file);
  return (DBInqVarExists(dbfile, patch_name) && 
          (DBInqVarType(dbfile, patch_name) == DB_QUAD_RECT));
}

void silo_file_write_amr_grid(silo_file_t* file, 
                              const char* grid_name,
                              amr_grid_t* grid)
{
  DBfile* dbfile = silo_file_dbfile(file);

  // This grid is really just a grouping of patches, as far as SILO is 
  // concerned.
  DBSetDir(dbfile, "/");
  int num_local_patches = amr_grid_num_local_patches(grid);
  char** patch_names = polymec_malloc(sizeof(char*) * num_local_patches);
  int patch_types[num_local_patches];
  int pos = 0, i, j, k, l = 0;
  int npx, npy, npz, nx, ny, nz, ng;
  amr_grid_get_extents(grid, &npx, &npy, &npz);
  amr_grid_get_patch_size(grid, &nx, &ny, &nz, &ng);
  sp_func_t* mapping = amr_grid_mapping(grid);
  while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k)) 
  {
    // Write out the grid for the patch itself.
    int i1 = i*npx, i2 = i*npx + nx, j1 = j*npy, j2 = j*npy + ny, k1 = k*npz, k2 = k*npz + nz;
    patch_names[l] = polymec_malloc(sizeof(char) * FILENAME_MAX);
    snprintf(patch_names[l], FILENAME_MAX, "%s_patch_%d_%d_%d", grid_name, i, j, k);
    write_amr_patch_grid(file, patch_names[l], i1, i2, j1, j2, k1, k2, mapping);

    // Name the thing and set its type.
    patch_types[i] = DB_QUAD_RECT;
  }
  ASSERT(l == num_local_patches);

  // Group all the patches together.
  DBPutMultimesh(dbfile, grid_name, num_local_patches, (const char* const*)patch_names, patch_types, NULL);

  // Clean up.
  for (int p = 0; p < num_local_patches; ++p)
    polymec_free(patch_names[p]);
  polymec_free(patch_names);

  // Add subdomain information for this grid.
  silo_file_add_subdomain_mesh(file, grid_name, DB_QUAD_RECT, NULL);
}

bool silo_file_contains_amr_grid(silo_file_t* file, 
                                 const char* grid_name)
{
  DBfile* dbfile = silo_file_dbfile(file);
  return (DBInqVarExists(dbfile, grid_name) && 
          (DBInqVarType(dbfile, grid_name) == DB_MULTIMESH));
}

void silo_file_write_amr_grid_data(silo_file_t* file, 
                                   const char** field_component_names,
                                   const char* grid_name,
                                   amr_grid_data_t* grid_data,
                                   silo_field_metadata_t** field_metadata)
{
  int num_local_patches = amr_grid_data_num_local_patches(grid_data);
  int num_components = amr_grid_data_num_components(grid_data);
  char field_names[num_components][FILENAME_MAX];
  int field_types[num_local_patches];
  amr_patch_t* patch;
  int pos = 0, i, j, k, l = 0;
  while (amr_grid_data_next(grid_data, &pos, &i, &j, &k, &patch))
  {
    // Write out the patch data itself.
    char patch_grid_name[FILENAME_MAX];
    snprintf(patch_grid_name, FILENAME_MAX-1, "%s_%d_%d_%d", grid_name, i, j, k);
    for (int c = 0; c < num_components; ++c)
      snprintf(field_names[c], FILENAME_MAX-1, "%s_%d_%d_%d", field_component_names[c], i, j, k);
    write_amr_patch_data(file, (const char**)field_names, patch_grid_name, patch, field_metadata);
    field_types[l] = DB_QUADVAR;
    ++l;
  }
  ASSERT(l == num_local_patches);

  // Finally, place multi-* entries into the Silo file.
  DBfile* dbfile = silo_file_dbfile(file);
  for (int c = 0; c < num_components; ++c)
  {
    DBPutMultivar(dbfile, field_component_names[c], num_local_patches, (const char* const*)field_names[c], field_types, NULL);

    // Add subdomain information for this component.
    silo_file_add_subdomain_mesh(file, field_component_names[c], DB_QUAD_RECT, NULL);
  }
}

void silo_file_write_amr_grid_hierarchy(silo_file_t* file, 
                                        const char* hierarchy_name,
                                        amr_grid_hierarchy_t* hierarchy)
{
}

bool silo_file_contains_amr_grid_hierarchy(silo_file_t* file, 
                                           const char* grid_hierarchy_name)
{
}

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
  int num_local_patches = 0;
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
      num_level_patches[i] = amr_grid_data_num_local_patches(grid_data);
      num_local_patches += num_level_patches[i];
      ++i;
    }
    level_region_names[0] = "@level%d@n";
    DBAddRegionArray(tree, num_levels, (const char* const*)level_region_names, 
                     0, level_maps_name, 1, seg_ids, num_level_patches, seg_types, NULL);
  }
  DBSetCwr(tree, "..");

  // Now define the patches.
  DBAddRegion(tree, "patches", 0, num_local_patches, NULL, 0, NULL, NULL, NULL, NULL); 
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
  DBfile* dbfile = silo_file_dbfile(file);
  DBPutMrgtree(dbfile, "mrgTree", "amr_mesh", tree, options);
  DBFreeMrgtree(tree);
}

bool silo_file_contains_amr_data_hierarchy(silo_file_t* file, 
                                           const char* data_hierarchy_name,
                                           const char* grid_hierarchy_name)
{
}
