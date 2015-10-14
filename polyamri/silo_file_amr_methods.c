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
                                 int N1, int N2, int N3, // index dimensions of containing space, if any
                                 int i1, int i2, int j1, int j2, int k1, int k2, // bounds of patch grid
                                 coord_mapping_t* mapping,
                                 bool hide_from_gui)
{
  ASSERT(i2 > i1);
  ASSERT(j2 > j1);
  ASSERT(k2 > k1);

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
  int coord_type;
  if (mapping != NULL)
  {
    coord_type = DB_NONCOLLINEAR;
    real_t dx1 = 1.0 / N1, dx2 = 1.0 / N2, dx3 = 1.0 / N3;
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
          point_t y;
          coord_mapping_map_point(mapping, &x, &y);
          x1_node[l] = y.x;
          x2_node[l] = y.y;
          x3_node[l] = y.z;
        }
      }
    }
  }
  else
  {
    coord_type = DB_COLLINEAR;
    real_t dx = 1.0 / N1;
    for (int i = i1; i <= i2; ++i)
      x1_node[i-i1] = i*dx;
    real_t dy = 1.0 / N2;
    for (int j = j1; j <= j2; ++j)
      x2_node[j-j1] = j*dy;
    real_t dz = 1.0 / N3;
    for (int k = k1; k <= k2; ++k)
      x3_node[k-k1] = k*dz;
  }

  real_t* coords[3];
  coords[0] = x1_node;
  coords[1] = x2_node;
  coords[2] = x3_node;

  // Write the patch grid.
  DBoptlist* optlist = NULL;
  int one = 1;
  if (hide_from_gui)
  {
    optlist = DBMakeOptlist(1); 
    DBAddOption(optlist, DBOPT_HIDE_FROM_GUI, &one);
  }
  DBPutQuadmesh(dbfile, patch_grid_name, coord_names, coords, dimensions, 3, 
                SILO_FLOAT_TYPE, coord_type, optlist);
  if (optlist != NULL)
    DBFreeOptlist(optlist);
}

static void write_amr_patch_data(silo_file_t* file,
                                 const char** field_component_names,
                                 const char* patch_grid_name,
                                 int N1, int N2, int N3, // index dimensions of containing space, if any
                                 amr_patch_t* patch,
                                 coord_mapping_t* mapping,
                                 silo_field_metadata_t** field_metadata,
                                 bool hide_from_gui)
{
  ASSERT(patch->i2 > patch->i1);
  ASSERT(patch->j2 > patch->j1);
  ASSERT(patch->k2 > patch->k1);

  DBfile* dbfile = silo_file_dbfile(file);
  int one = 1;

  // If we're mapped and there are vector fields, we need to treat them specially.
  int num_vector_components = 0;
  bool is_vector_component[patch->nc];
  int first_vector_component = -1;
  if ((mapping != NULL) && (field_metadata != NULL))
  {
    for (int c = 0; c < patch->nc; ++c)
    {
      if ((field_metadata[c] != NULL) && silo_field_metadata_is_vector_component(field_metadata[c]))
      {
        if ((num_vector_components % 3) == 0)
        {
          ASSERT(patch->nc >= c + 3); // If you start a vector, you better finish your vector.
          first_vector_component = c;
        }
        is_vector_component[c] = true;
        ++num_vector_components;
      }
    }
  }
  else
    memset(is_vector_component, 0, sizeof(bool) * patch->nc);
  ASSERT((num_vector_components % 3) == 0); // We should have complete sets of vector triples.

  // Write the data.
  int n1 = patch->i2 - patch->i1;
  int n2 = patch->j2 - patch->j1;
  int n3 = patch->k2 - patch->k1;
  real_t* data = polymec_malloc(sizeof(real_t) * n1*n2*n3);
  int cell_dims[3] = {n1, n2, n3};
  DECLARE_AMR_PATCH_ARRAY(a, patch);
  for (int c = 0; c < patch->nc; ++c)
  {
    DBoptlist* optlist = (field_metadata != NULL) ? optlist_from_metadata(field_metadata[c]) : NULL;
    if (hide_from_gui)
    {
      if (optlist == NULL)
        optlist = DBMakeOptlist(1); 
      DBAddOption(optlist, DBOPT_HIDE_FROM_GUI, &one);
    }
    if ((mapping != NULL) && is_vector_component[c])
    {
      // We need to map this vector field before we write it out.
      int c1 = first_vector_component, c2 = first_vector_component+1, c3 = first_vector_component+2;
      int which_component = c - c1;
      int l = 0;
      point_t x;
      real_t dx = 1.0 / N1, dy = 1.0 / N2, dz = 1.0 / N3;
      for (int i = 0; i < n1; ++i)
      {
        x.x = (i+patch->i1) * dx;
        for (int j = 0; j < n2; ++j)
        {
          x.y = (j+patch->j1) * dy;
          for (int k = 0; k < n3; ++k, ++l)
          {
            x.z = (k+patch->k1) * dz;
            vector_t v = {.x = a[i][j][k][c1], a[i][j][k][c2], a[i][j][k][c3]};
            vector_t v1;
            coord_mapping_map_vector(mapping, &x, &v, &v1);
            switch (which_component)
            {
              case 0: data[l] = v1.x; break;
              case 1: data[l] = v1.y; break;
              default: data[l] = v1.z;
            }
          }
        }
      }
    }
    else
    {
      // Copy the field data verbatim.
      int l = 0;
      for (int i = patch->i1; i < patch->i2; ++i)
        for (int j = patch->j1; j < patch->j2; ++j)
          for (int k = patch->k1; k < patch->k2; ++k, ++l)
            data[l] = a[i][j][k][c];
    }
    DBPutQuadvar1(dbfile, field_component_names[c], patch_grid_name, data,
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
                                      coord_mapping_t* mapping)
{
  int N1 = patch->i2 - patch->i1, N2 = patch->j2 - patch->j1, N3 = patch->k2 - patch->k1; 
  write_amr_patch_grid(file, patch_grid_name, N1, N2, N3, patch->i1, patch->i2, 
                       patch->j1, patch->j2, patch->k1, patch->k2, mapping, false);
  write_amr_patch_data(file, field_component_names, patch_grid_name, N1, N2, N3,
                       patch, mapping, NULL, false);
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
  int npx, npy, npz, nx, ny, nz, ng;
  amr_grid_get_extents(grid, &npx, &npy, &npz);
  amr_grid_get_patch_size(grid, &nx, &ny, &nz, &ng);

  coord_mapping_t* mapping = amr_grid_mapping(grid);
  char* patch_grid_names[num_local_patches];
  int patch_grid_types[num_local_patches];
  int N1 = npx * nx, N2 = npy * ny, N3 = npz * nz;
  int pos = 0, i, j, k, l = 0;
  while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k)) 
  {
    // Write out the grid for the patch itself.
    int i1 = nx*i, i2 = nx*(i+1), j1 = ny*j, j2 = ny*(j+1), k1 = nz*k, k2 = nz*(k+1);
    char patch_grid_name[FILENAME_MAX];
    snprintf(patch_grid_name, FILENAME_MAX-1, "%s_%d_%d_%d", grid_name, i, j, k);
    patch_grid_names[l] = string_dup(patch_grid_name);
    patch_grid_types[l] = DB_QUAD_RECT;
    write_amr_patch_grid(file, patch_grid_names[l], N1, N2, N3, i1, i2, j1, j2, k1, k2, mapping, true);
    ++l;
  }
  ASSERT(l == num_local_patches);

  // Group all the patches together.
  DBPutMultimesh(dbfile, grid_name, num_local_patches, (const char* const*)patch_grid_names, patch_grid_types, NULL);

  // Clean up.
  for (int p = 0; p < num_local_patches; ++p)
    string_free(patch_grid_names[p]);

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

  amr_grid_t* grid = amr_grid_data_grid(grid_data);
  coord_mapping_t* mapping = amr_grid_mapping(grid);
  int npx, npy, npz, nx, ny, nz, ng;
  amr_grid_get_extents(grid, &npx, &npy, &npz);
  amr_grid_get_patch_size(grid, &nx, &ny, &nz, &ng);

  char* field_names[num_components];
  char* multi_field_names[num_components][num_local_patches];
  int multi_field_types[num_local_patches];
  int N1 = npx * nx, N2 = npy * ny, N3 = npz * nz;

  char* patch_grid_names[num_local_patches];
  int patch_grid_types[num_local_patches];
  amr_patch_t* patch;
  int pos = 0, i, j, k, l = 0;
  while (amr_grid_data_next_local_patch(grid_data, &pos, &i, &j, &k, &patch))
  {
    // Write out the patch data itself.
    for (int c = 0; c < num_components; ++c)
    {
      char field_name[FILENAME_MAX];
      snprintf(field_name, FILENAME_MAX-1, "%s_%d_%d_%d", field_component_names[c], i, j, k);
      field_names[c] = string_dup(field_name);
      multi_field_names[c][l] = string_dup(field_name);
    }
    char patch_grid_name[FILENAME_MAX];
    snprintf(patch_grid_name, FILENAME_MAX-1, "%s_%d_%d_%d", grid_name, i, j, k);
    write_amr_patch_data(file, (const char**)field_names, patch_grid_name, N1, N2, N3, 
                         patch, mapping, field_metadata, true);
    multi_field_types[l] = DB_QUADVAR;
    ++l;

    for (int c = 0; c < num_components; ++c)
      string_free(field_names[c]);
  }
  ASSERT(l == num_local_patches);

  // Finally, place multi-* entries into the Silo file.
  DBfile* dbfile = silo_file_dbfile(file);
  for (int c = 0; c < num_components; ++c)
  {
    // We need to associate this multi-variable with our multi-mesh.
    DBoptlist* optlist = DBMakeOptlist(4);
    DBAddOption(optlist, DBOPT_MMESH_NAME, (void*)grid_name);
    DBPutMultivar(dbfile, field_component_names[c], num_local_patches, (const char* const*)multi_field_names[c], multi_field_types, NULL);
    DBFreeOptlist(optlist);

    // Add subdomain information for this component.
    silo_file_add_subdomain_mesh(file, field_component_names[c], DB_QUAD_RECT, NULL);
  }

  // Clean up.
  for (int c = 0; c < num_components; ++c)
    for (int p = 0; p < num_local_patches; ++p)
      string_free(multi_field_names[c][p]);
}

void silo_file_write_amr_grid_hierarchy(silo_file_t* file, 
                                        const char* hierarchy_name,
                                        amr_grid_hierarchy_t* hierarchy)
{
  // First, write out all the refinement levels in the hierarchy.
  int pos = 0, l = 0;
  amr_grid_t* grid;
  while (amr_grid_hierarchy_next_coarsest(hierarchy, &pos, &grid))
  {
    char grid_name[FILENAME_MAX];
    snprintf(grid_name, FILENAME_MAX-1, "%s_level_%d", hierarchy_name, l);
    silo_file_write_amr_grid(file, grid_name, grid);
  }
}

bool silo_file_contains_amr_grid_hierarchy(silo_file_t* file, 
                                           const char* grid_hierarchy_name)
{
}

void silo_file_write_amr_data_hierarchy(silo_file_t* file, 
                                        const char** field_component_names,
                                        const char* grid_hierarchy_name,
                                        amr_data_hierarchy_t* data_hierarchy,
                                        silo_field_metadata_t** field_metadata)
{
  int num_components = amr_data_hierarchy_num_components(data_hierarchy);

  // First, write out all the refinement levels in the hierarchy.
  int pos = 0, l = 0;
  amr_grid_data_t* data;
  while (amr_data_hierarchy_next_coarsest(data_hierarchy, &pos, &data))
  {
    char* field_names[num_components];
    for (int c = 0; c < num_components; ++c)
    {
      char field_name[FILENAME_MAX];
      snprintf(field_name, FILENAME_MAX-1, "%s_level_%d", field_component_names[c], l);
      field_names[c] = string_dup(field_name);
    }
    char grid_name[FILENAME_MAX];
    snprintf(grid_name, FILENAME_MAX-1, "%s_level_%d", grid_hierarchy_name, l);
    silo_file_write_amr_grid_data(file, (const char**)field_names, grid_name, data, field_metadata);
    ++l;

    // Clean up.
    for (int c = 0; c < num_components; ++c)
      string_free(field_names[c]);
  }

#if 0
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
#endif
}

bool silo_file_contains_amr_data_hierarchy(silo_file_t* file, 
                                           const char* data_hierarchy_name,
                                           const char* grid_hierarchy_name)
{
}
