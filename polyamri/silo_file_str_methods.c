// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "silo.h"
#include "core/silo_file.h"
#include "polyamri/silo_file_str_methods.h"

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

static void write_str_patch_grid(silo_file_t* file,
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
  int N = (n1+1) * (n2+1) * (n3+1);
  real_t x1_node[N], x2_node[N], x3_node[N];
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

static void write_str_patch_data(silo_file_t* file,
                                 const char** field_component_names,
                                 const char* patch_grid_name,
                                 str_grid_patch_t* patch,
                                 silo_field_metadata_t** field_metadata,
                                 int N1, int N2, int N3,
                                 coord_mapping_t* mapping,
                                 bool hide_from_gui)
{
  ASSERT(patch->i2 > patch->i1);
  ASSERT(patch->j2 > patch->j1);
  ASSERT(patch->k2 > patch->k1);

  DBfile* dbfile = silo_file_dbfile(file);
  int one = 1;

  // If we're given a mapping and there are vector fields, we need to treat 
  // them specially.
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
  DECLARE_STR_GRID_PATCH_ARRAY(a, patch);
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

void silo_file_write_str_grid_patch(silo_file_t* file, 
                                    const char** field_component_names,
                                    const char* patch_grid_name,
                                    str_grid_patch_t* patch,
                                    coord_mapping_t* mapping,
                                    coord_mapping_t* data_mapping)
{
  int N1 = patch->i2 - patch->i1, N2 = patch->j2 - patch->j1, N3 = patch->k2 - patch->k1; 
  write_str_patch_grid(file, patch_grid_name, N1, N2, N3, patch->i1, patch->i2, 
                       patch->j1, patch->j2, patch->k1, patch->k2, mapping, false);
  write_str_patch_data(file, field_component_names, patch_grid_name, 
                       patch, NULL, N1, N2, N3, data_mapping, false);
}

bool silo_file_contains_str_grid_patch(silo_file_t* file, 
                                       const char* patch_name)
{
  DBfile* dbfile = silo_file_dbfile(file);
  return (DBInqVarExists(dbfile, patch_name) && 
          (DBInqVarType(dbfile, patch_name) == DB_QUAD_RECT));
}

void silo_file_write_str_grid(silo_file_t* file, 
                              const char* grid_name,
                              str_grid_t* grid,
                              coord_mapping_t* mapping)
{
  DBfile* dbfile = silo_file_dbfile(file);

  // This grid is really just a grouping of patches, as far as SILO is 
  // concerned.
  DBSetDir(dbfile, "/");
  int num_local_patches = str_grid_num_patches(grid);
  int npx, npy, npz, nx, ny, nz;
  str_grid_get_extents(grid, &npx, &npy, &npz);
  str_grid_get_patch_size(grid, &nx, &ny, &nz);

  char* patch_grid_names[num_local_patches];
  int patch_grid_types[num_local_patches];
  int N1 = npx * nx, N2 = npy * ny, N3 = npz * nz;
  int pos = 0, i, j, k, l = 0;
  while (str_grid_next_patch(grid, &pos, &i, &j, &k)) 
  {
    // Write out the grid for the patch itself.
    int i1 = nx*i, i2 = nx*(i+1), j1 = ny*j, j2 = ny*(j+1), k1 = nz*k, k2 = nz*(k+1);
    char patch_grid_name[FILENAME_MAX];
    snprintf(patch_grid_name, FILENAME_MAX-1, "%s_%d_%d_%d", grid_name, i, j, k);
    patch_grid_names[l] = string_dup(patch_grid_name);
    patch_grid_types[l] = DB_QUAD_RECT;
    write_str_patch_grid(file, patch_grid_names[l], N1, N2, N3, i1, i2, j1, j2, k1, k2, mapping, true);
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

bool silo_file_contains_str_grid(silo_file_t* file, 
                                 const char* grid_name)
{
  DBfile* dbfile = silo_file_dbfile(file);
  return (DBInqVarExists(dbfile, grid_name) && 
          (DBInqVarType(dbfile, grid_name) == DB_MULTIMESH));
}

void silo_file_write_str_grid_cell_data(silo_file_t* file, 
                                        const char** field_component_names,
                                        const char* grid_name,
                                        str_grid_cell_data_t* cell_data,
                                        silo_field_metadata_t** field_metadata,
                                        coord_mapping_t* mapping)
{
  int num_local_patches = str_grid_cell_data_num_patches(cell_data);
  int num_components = str_grid_cell_data_num_components(cell_data);

  str_grid_t* grid = str_grid_cell_data_grid(cell_data);
  int npx, npy, npz, nx, ny, nz;
  str_grid_get_extents(grid, &npx, &npy, &npz);
  str_grid_get_patch_size(grid, &nx, &ny, &nz);

  char* field_names[num_components];
  char* multi_field_names[num_components][num_local_patches];
  int multi_field_types[num_local_patches];

  str_grid_patch_t* patch;
  int pos = 0, i, j, k, l = 0;
  while (str_grid_cell_data_next_patch(cell_data, &pos, &i, &j, &k, &patch, NULL))
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
    int N1 = patch->i2 - patch->i1, N2 = patch->j2 - patch->j1, N3 = patch->k2 - patch->k1; 
    write_str_patch_data(file, (const char**)field_names, patch_grid_name,  
                         patch, field_metadata, N1, N2, N3, mapping, true);
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

bool silo_file_contains_str_grid_cell_data(silo_file_t* file, 
                                           const char* cell_data_name,
                                           const char* grid_name)
{
  // FIXME
  DBfile* dbfile = silo_file_dbfile(file);
  return (DBInqVarExists(dbfile, grid_name) && 
          (DBInqVarType(dbfile, grid_name) == DB_MULTIMESH));
}

