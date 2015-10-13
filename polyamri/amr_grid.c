// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/unordered_map.h"
#include "polyamri/amr_grid.h"
#include "polyamri/amr_grid_data.h"
#include "polyamri/amr_grid_interpolator.h"

// This is a descriptor for "types" of patches -- this determines the nature 
// of interacting patches during ghost cell filling.
typedef enum
{
  NO_PATCH,
  LOCAL_SAME_LEVEL,
  LOCAL_FINER_LEVEL,
  LOCAL_COARSER_LEVEL,
  REMOTE_SAME_LEVEL,
  REMOTE_FINER_LEVEL,
  REMOTE_COARSER_LEVEL
} patch_type_t;

// This stores a function and its arguments for starting the process of 
// filling ghost cells.
typedef struct
{
  int (*method)(amr_grid_t*, amr_patch_t*, int, int, int, amr_patch_t*, int, int, int);
  int i_dest, j_dest, k_dest;
  int i_src, j_src, k_src;
} amr_grid_ghost_filler_starter_t;

// This stores a function and its arguments for finishing the process of 
// filling ghost cells.
typedef struct
{
  void (*method)(amr_grid_t*, int, amr_patch_t*, int, int, int, amr_patch_t*, int, int, int);
  int token;
  amr_patch_t* dest_patch;
  int i_dest, j_dest, k_dest;
  amr_patch_t* src_patch;
  int i_src, j_src, k_src;
} amr_grid_ghost_filler_finisher_t;

struct amr_grid_t
{
  // Intrinsic metadata.
  int nx, ny, nz, px, py, pz, ng; 
  bool x_periodic, y_periodic, z_periodic;
  
  // Information about which patches are present.
  int num_local_patches;
  int* local_patch_indices;
  patch_type_t* patch_types;
  int* remote_owners;

  // Mapping function.
  sp_func_t* mapping;

  // Associated properties.
  string_ptr_unordered_map_t* properties;

  // Ratio of refinement between this grid and its "neighbors"
  int ref_ratio; 

  // Coarser grid and associated interpolator.
  amr_grid_t* coarser;
  amr_grid_interpolator_t* coarse_interpolator;

  // Finer grid and associated interpolator.
  amr_grid_t* finer;
  amr_grid_interpolator_t* fine_interpolator;

  // Neighboring grids and interpolators.
  amr_grid_t* neighbors[6];
  amr_grid_interpolator_t* neighbor_interpolators[6];

  // Ghost filling machinery -- asynchronous!
  ptr_array_t* ghost_filler_starters;
  ptr_array_t* ghost_filler_finishers;

  // This flag is set by amr_grid_finalize() after a grid has been assembled.
  bool finalized;
};

amr_grid_t* amr_grid_new(int nx, int ny, int nz, 
                         int px, int py, int pz, 
                         int num_ghosts, 
                         bool periodic_in_x, 
                         bool periodic_in_y, 
                         bool periodic_in_z)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  ASSERT(px > 0);
  ASSERT(py > 0);
  ASSERT(pz > 0);
  ASSERT(num_ghosts >= 0);
  ASSERT(num_ghosts <= px);
  ASSERT(num_ghosts <= py);
  ASSERT(num_ghosts <= pz);

  amr_grid_t* grid = polymec_malloc(sizeof(amr_grid_t));
  grid->nx = nx;
  grid->ny = ny;
  grid->nz = nz;
  grid->px = px;
  grid->py = py;
  grid->pz = pz;
  grid->ng = num_ghosts;
  grid->num_local_patches = 0;
  grid->local_patch_indices = NULL;
  grid->x_periodic = periodic_in_x;
  grid->y_periodic = periodic_in_y;
  grid->z_periodic = periodic_in_z;
  grid->patch_types = polymec_malloc(sizeof(patch_type_t) * nx * ny * nz);
  grid->remote_owners = polymec_malloc(sizeof(int) * nx * ny * nz);
  for (int i = 0; i < nx * ny * nz; ++i)
  {
    grid->patch_types[i] = NO_PATCH;
    grid->remote_owners[i] = -1;
  }

  grid->ref_ratio = 0;
  grid->coarser = NULL;
  grid->coarse_interpolator = NULL;
  grid->finer = NULL;
  grid->fine_interpolator = NULL;

  for (int n = 0; n < 6; ++n)
  {
    grid->neighbors[n] = NULL;
    grid->neighbor_interpolators[n] = NULL;
  }

  grid->properties = string_ptr_unordered_map_new();
  grid->ghost_filler_starters = ptr_array_new();
  grid->ghost_filler_finishers = ptr_array_new();
  grid->finalized = false;

  return grid;
}

void amr_grid_free(amr_grid_t* grid)
{
  for (int n = 0; n < 6; ++n)
  {
    if (grid->neighbor_interpolators[n] != NULL)
      amr_grid_interpolator_free(grid->neighbor_interpolators[n]);
  }
  if (grid->coarse_interpolator != NULL)
    amr_grid_interpolator_free(grid->coarse_interpolator);
  if (grid->fine_interpolator != NULL)
    amr_grid_interpolator_free(grid->fine_interpolator);
  polymec_free(grid->patch_types);
  polymec_free(grid->remote_owners);
  string_ptr_unordered_map_free(grid->properties);
  ptr_array_free(grid->ghost_filler_starters);
  ptr_array_free(grid->ghost_filler_finishers);
  if (grid->local_patch_indices != NULL)
    polymec_free(grid->local_patch_indices);
  polymec_free(grid);
}

void amr_grid_set_mapping(amr_grid_t* grid, sp_func_t* mapping)
{
  ASSERT((mapping == NULL) || (sp_func_num_comp(mapping) == 3));
  grid->mapping = mapping;
}

sp_func_t* amr_grid_mapping(amr_grid_t* grid)
{
  return grid->mapping;
}

void amr_grid_set_neighbor(amr_grid_t* grid, 
                           amr_grid_neighbor_slot_t neighbor_slot,
                           amr_grid_t* neighbor,
                           amr_grid_interpolator_t* interpolator)
{
  ASSERT(!(grid->x_periodic && (neighbor_slot == AMR_GRID_X1_NEIGHBOR)));
  ASSERT(!(grid->x_periodic && (neighbor_slot == AMR_GRID_X2_NEIGHBOR)));
  ASSERT(!(grid->y_periodic && (neighbor_slot == AMR_GRID_Y1_NEIGHBOR)));
  ASSERT(!(grid->y_periodic && (neighbor_slot == AMR_GRID_Y2_NEIGHBOR)));
  ASSERT(!(grid->z_periodic && (neighbor_slot == AMR_GRID_Z1_NEIGHBOR)));
  ASSERT(!(grid->z_periodic && (neighbor_slot == AMR_GRID_Z2_NEIGHBOR)));
  grid->neighbors[neighbor_slot] = neighbor;
  grid->neighbor_interpolators[neighbor_slot] = interpolator;
}

void amr_grid_set_property(amr_grid_t* grid,
                           const char* name,
                           void* property,
                           void (*property_dtor)(void* property))
{
  string_ptr_unordered_map_insert_with_kv_dtors(grid->properties, string_dup(name), property, string_free, property_dtor);
}

void* amr_grid_property(amr_grid_t* grid,
                        const char* name)
{
  void** ptr = string_ptr_unordered_map_get(grid->properties, (char*)name);
  if (ptr != NULL)
    return *ptr;
  else
    return NULL;
}

void amr_grid_associate_finer(amr_grid_t* grid, 
                              amr_grid_t* finer_grid, 
                              int ref_ratio,
                              amr_grid_interpolator_t* interpolator)
{
  ASSERT(!grid->finalized);
  ASSERT(ref_ratio > 0);
  ASSERT((grid->ref_ratio == 0) || (grid->ref_ratio == ref_ratio));

  if (grid->ref_ratio == 0)
    grid->ref_ratio = ref_ratio;

  grid->finer = finer_grid;
  grid->fine_interpolator = interpolator;
}

void amr_grid_associate_coarser(amr_grid_t* grid, 
                                amr_grid_t* coarser_grid, 
                                int ref_ratio,
                                amr_grid_interpolator_t* interpolator)
{
  ASSERT(!grid->finalized);
  ASSERT(ref_ratio > 0);
  ASSERT((grid->ref_ratio == 0) || (grid->ref_ratio == ref_ratio));

  if (grid->ref_ratio == 0)
    grid->ref_ratio = ref_ratio;

  grid->coarser = coarser_grid;
  grid->coarse_interpolator = interpolator;
}

void amr_grid_add_local_patch(amr_grid_t* grid, int i, int j, int k)
{
  ASSERT(!grid->finalized);
  int patch_index = grid->nz*grid->ny*i + grid->nz*j + k;
  patch_type_t patch_type = grid->patch_types[patch_index];
  ASSERT(patch_type == NO_PATCH);
  grid->patch_types[patch_index] = LOCAL_SAME_LEVEL;
  ++(grid->num_local_patches);
}

void amr_grid_add_remote_patch(amr_grid_t* grid, int i, int j, int k, int remote_owner)
{
  ASSERT(!grid->finalized);
  int patch_index = grid->nz*grid->ny*i + grid->nz*j + k;
  patch_type_t patch_type = grid->patch_types[patch_index];
  ASSERT(patch_type == NO_PATCH);
  grid->patch_types[patch_index] = REMOTE_SAME_LEVEL;
  grid->remote_owners[patch_index] = remote_owner;
}

void amr_grid_finalize(amr_grid_t* grid)
{
  ASSERT(!grid->finalized);
  ASSERT(ptr_array_empty(grid->ghost_filler_starters));
  ASSERT(ptr_array_empty(grid->ghost_filler_finishers));

  // Make sure that we have accounted for all patches in our construction 
  // process.
  if (grid->num_local_patches < (grid->nx * grid->ny * grid->nz))
  {
    for (int i = 0; i < grid->nx; ++i)
    {
      for (int j = 0; j < grid->ny; ++j)
      {
        for (int k = 0; k < grid->nz; ++k)
        {
          int patch_index = grid->nz*grid->ny*i + grid->nz*j + k;
          if (grid->patch_types[patch_index] == NO_PATCH)
            polymec_error("amr_grid_finalize: no local or remote patch at (%d, %d, %d).", i, j, k);
        }
      }
    }
  }

  // Make an array of indices for locally-present patches.
  grid->local_patch_indices = polymec_malloc(sizeof(int) * 3 * grid->num_local_patches);
  int l = 0;
  for (int i = 0; i < grid->nx; ++i)
  {
    for (int j = 0; j < grid->ny; ++j)
    {
      for (int k = 0; k < grid->nz; ++k, ++l)
      {
        grid->local_patch_indices[3*l]   = i;
        grid->local_patch_indices[3*l+1] = j;
        grid->local_patch_indices[3*l+2] = k;
      }
    }
  }

#if 0
  // Assemble our ghost fillers.
  for (int i = 0; i < grid->nx; ++i)
  {
    for (int j = 0; j < grid->ny; ++j)
    {
      for (int k = 0; k < grid->nz; ++k)
      {
        patch_type_t patch_type = get_patch_type(grid, i, j, k);
        amr_grid_ghost_filler_t* ghost_filler = polymec_malloc(sizeof(amr_grid_ghost_filler_t));
        switch(patch_type)
        {
          case LOCAL_SAME_LEVEL:
            local_ghost_copy(grid, dest, i_dest, j_dest, k_dest, 
                src, i_src, j_src, k_src);
            break;
          case LOCAL_FINER_LEVEL:
            local_fine_to_coarse_interpolation(grid, dest, i_dest, j_dest, k_dest, 
                src, i_src, j_src, k_src);
            break;
          case LOCAL_COARSER_LEVEL:
            local_coarse_to_fine_interpolation(grid, dest, i_dest, j_dest, k_dest, 
                src, i_src, j_src, k_src);
            break;
          case REMOTE_SAME_LEVEL:
            token = start_remote_ghost_copy(grid, dest, i_dest, j_dest, k_dest, 
                src, i_src, j_src, k_src);
            break;
          case REMOTE_FINER_LEVEL:
            token = start_remote_fine_to_coarse_interpolation(grid, dest, i_dest, j_dest, k_dest, 
                src, i_src, j_src, k_src);
            break;
          case REMOTE_COARSER_LEVEL:
            token = start_remote_coarse_to_fine_interpolation(grid, dest, i_dest, j_dest, k_dest, 
                src, i_src, j_src, k_src);
          default:
            break;
        }
        ghost_filler->i_dest = i_dest;
        ptr_array_append_with_dtor(grid->ghost_fillers, ghost_filler, polymec_free);
      }
    }
  }
#endif

  grid->finalized = true;
}

bool amr_grid_next_local_patch(amr_grid_t* grid, int* pos, int* i, int* j, int* k)
{
  ASSERT(grid->finalized);
  ASSERT(*pos >= 0);
  bool result = (*pos < grid->num_local_patches);
  if (result)
  {
    int l = *pos;
    *i = grid->local_patch_indices[3*l];
    *j = grid->local_patch_indices[3*l+1];
    *k = grid->local_patch_indices[3*l+2];
    ++(*pos);
  }
  return result;
}

void amr_grid_get_periodicity(amr_grid_t* grid, bool* periodicity)
{
  ASSERT(periodicity != NULL);
  periodicity[0] = grid->x_periodic;
  periodicity[1] = grid->y_periodic;
  periodicity[2] = grid->z_periodic;
}

void amr_grid_get_extents(amr_grid_t* grid, int* nx, int* ny, int* nz)
{
  *nx = grid->nx;
  *ny = grid->ny;
  *nz = grid->nz;
}

void amr_grid_get_patch_size(amr_grid_t* grid, int* pnx, int* pny, int* pnz, int* png)
{
  *pnx = grid->px;
  *pny = grid->py;
  *pnz = grid->pz;
  *png = grid->ng;
}

int amr_grid_num_local_patches(amr_grid_t* grid)
{
  return grid->num_local_patches;
}

static inline patch_type_t get_patch_type(amr_grid_t* grid, int i, int j, int k)
{
  return grid->patch_types[grid->nz*grid->ny*i + grid->nz*j + k];
}

bool amr_grid_has_patch(amr_grid_t* grid, int i, int j, int k)
{
  return (get_patch_type(grid, i, j, k) == LOCAL_SAME_LEVEL);
}

static int local_copy_west_to_east(amr_grid_t* grid, 
                                   amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                   amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT((i_src == i_dest - 1) || ((i_dest == 0) && (i_src == grid->nx-1)));
  ASSERT(j_src == j_dest);
  ASSERT(k_src == k_dest);
  ASSERT(dest->nc == src->nc);
  ASSERT(dest->i1 == src->i1);
  ASSERT(dest->i2 == src->i2);
  ASSERT(dest->j1 == src->j1);
  ASSERT(dest->j2 == src->j2);
  ASSERT(dest->k1 == src->k1);
  ASSERT(dest->k2 == src->k2);

  // Get multi-arrays for the source and destination patches.
  DECLARE_AMR_PATCH_ARRAY(s, src);
  DECLARE_AMR_PATCH_ARRAY(d, dest);

  int nc = src->nc, ng = grid->ng;
  int i1 = src->i1, i2 = src->i2, 
      j1 = src->j1, j2 = src->j2,
      k1 = src->k1, k2 = src->k2;
  for (int jj = j1; jj < j2; ++jj)
    for (int kk = k1; kk < k2; ++kk)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[i1-g-1][jj][kk][c] = s[i2-g-1][jj][kk][c];
  return 0;
}

static int local_copy_east_to_west(amr_grid_t* grid, 
                                   amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                   amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT((i_src == i_dest + 1) || ((i_dest == grid->nx-1) && (i_src == 0)));
  ASSERT(j_src == j_dest);
  ASSERT(k_src == k_dest);
  ASSERT(dest->nc == src->nc);
  ASSERT(dest->i1 == src->i1);
  ASSERT(dest->i2 == src->i2);
  ASSERT(dest->j1 == src->j1);
  ASSERT(dest->j2 == src->j2);
  ASSERT(dest->k1 == src->k1);
  ASSERT(dest->k2 == src->k2);

  // Get multi-arrays for the source and destination patches.
  DECLARE_AMR_PATCH_ARRAY(s, src);
  DECLARE_AMR_PATCH_ARRAY(d, dest);

  int nc = src->nc, ng = grid->ng;
  int i1 = src->i1, i2 = src->i2, 
      j1 = src->j1, j2 = src->j2,
      k1 = src->k1, k2 = src->k2;
  for (int jj = j1; jj < j2; ++jj)
    for (int kk = k1; kk < k2; ++kk)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[i2+g][jj][kk][c] = s[i1+ng-g-1][jj][kk][c];
  return 0;
}

static int local_copy_south_to_north(amr_grid_t* grid, 
                                     amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                     amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT(i_src == i_dest);
  ASSERT((j_src == j_dest - 1) || ((j_dest == 0) && (j_src == grid->ny-1)));
  ASSERT(k_src == k_dest);
  ASSERT(dest->nc == src->nc);
  ASSERT(dest->i1 == src->i1);
  ASSERT(dest->i2 == src->i2);
  ASSERT(dest->j1 == src->j1);
  ASSERT(dest->j2 == src->j2);
  ASSERT(dest->k1 == src->k1);
  ASSERT(dest->k2 == src->k2);

  // Get multi-arrays for the source and destination patches.
  DECLARE_AMR_PATCH_ARRAY(s, src);
  DECLARE_AMR_PATCH_ARRAY(d, dest);

  int nc = src->nc, ng = grid->ng;
  int i1 = src->i1, i2 = src->i2, 
      j1 = src->j1, j2 = src->j2,
      k1 = src->k1, k2 = src->k2;
  for (int ii = i1; ii < i2; ++ii)
    for (int kk = k1; kk < k2; ++kk)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[ii][j1-g-1][kk][c] = s[ii][j2-g-1][kk][c];

  return 0;
}

static int local_copy_north_to_south(amr_grid_t* grid, 
                                     amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                     amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT(i_src == i_dest);
  ASSERT((j_src == j_dest + 1) || ((j_dest == grid->ny-1) && (j_src == 0)));
  ASSERT(k_src == k_dest);
  ASSERT(dest->nc == src->nc);
  ASSERT(dest->i1 == src->i1);
  ASSERT(dest->i2 == src->i2);
  ASSERT(dest->j1 == src->j1);
  ASSERT(dest->j2 == src->j2);
  ASSERT(dest->k1 == src->k1);
  ASSERT(dest->k2 == src->k2);

  // Get multi-arrays for the source and destination patches.
  DECLARE_AMR_PATCH_ARRAY(s, src);
  DECLARE_AMR_PATCH_ARRAY(d, dest);

  int nc = src->nc, ng = grid->ng;
  int i1 = src->i1, i2 = src->i2, 
      j1 = src->j1, j2 = src->j2,
      k1 = src->k1, k2 = src->k2;
  for (int ii = i1; ii < i2; ++ii)
    for (int kk = k1; kk < k2; ++kk)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[ii][j2+g][kk][c] = s[ii][j1+ng-g-1][kk][c];

  return 0;
}

static int local_copy_below_to_above(amr_grid_t* grid, 
                                     amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                     amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT(i_src == i_dest);
  ASSERT(j_src == j_dest);
  ASSERT((k_src == k_dest - 1) || ((k_dest == 0) && (k_src == grid->nz-1)));
  ASSERT(dest->nc == src->nc);
  ASSERT(dest->i1 == src->i1);
  ASSERT(dest->i2 == src->i2);
  ASSERT(dest->j1 == src->j1);
  ASSERT(dest->j2 == src->j2);
  ASSERT(dest->k1 == src->k1);
  ASSERT(dest->k2 == src->k2);

  // Get multi-arrays for the source and destination patches.
  DECLARE_AMR_PATCH_ARRAY(s, src);
  DECLARE_AMR_PATCH_ARRAY(d, dest);

  int nc = src->nc, ng = grid->ng;
  int i1 = src->i1, i2 = src->i2, 
      j1 = src->j1, j2 = src->j2,
      k1 = src->k1, k2 = src->k2;
  for (int ii = i1; ii < i2; ++ii)
    for (int jj = j1; jj < j2; ++jj)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[ii][jj][k1-g-1][c] = s[ii][jj][k2-g-1][c];

  return 0;
}

static int local_copy_above_to_below(amr_grid_t* grid, 
                                     amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                     amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT(i_src == i_dest);
  ASSERT(j_src == j_dest);
  ASSERT((k_src == k_dest + 1) || ((k_dest == grid->nz-1) && (k_src == 0)));
  ASSERT(dest->nc == src->nc);
  ASSERT(dest->i1 == src->i1);
  ASSERT(dest->i2 == src->i2);
  ASSERT(dest->j1 == src->j1);
  ASSERT(dest->j2 == src->j2);
  ASSERT(dest->k1 == src->k1);
  ASSERT(dest->k2 == src->k2);

  // Get multi-arrays for the source and destination patches.
  DECLARE_AMR_PATCH_ARRAY(s, src);
  DECLARE_AMR_PATCH_ARRAY(d, dest);

  int nc = src->nc, ng = grid->ng;
  int i1 = src->i1, i2 = src->i2, 
      j1 = src->j1, j2 = src->j2,
      k1 = src->k1, k2 = src->k2;
  for (int ii = i1; ii < i2; ++ii)
    for (int jj = j1; jj < j2; ++jj)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[ii][jj][k2+g][c] = s[ii][jj][k1+ng-g-1][c];

  return 0;
}

static int local_fine_to_coarse_west_to_east(amr_grid_t* grid, 
                                             amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                             amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT((i_src == i_dest - 1) || ((i_dest == 0) && (i_src == grid->nx-1)));
  ASSERT(j_src == j_dest);
  ASSERT(k_src == k_dest);
  ASSERT(grid->finer != NULL);
  ASSERT(grid->ref_ratio > 0);
  ASSERT(dest->nc == src->nc);

  int nc = dest->nc, ng = grid->ng;
  int ref_ratio = grid->ref_ratio;

  // Loop over all the fine patches that correspond to this coarser one, and 
  // perform the interpolation.
  for (int I = 0; I < ref_ratio; ++I)
  {
    int i_fine = ref_ratio * i_src + I; 
    for (int J = 0; J < ref_ratio; ++J)
    {
      int j_fine = ref_ratio * j_src + J;
      for (int K = 0; K < ref_ratio; ++K)
      {
        int k_fine = ref_ratio * k_src;
        amr_grid_interpolator_t* interp; // = get_interpolator(grid->finer, i_fine, j_fine, k_fine);
        if (interp == NULL) // no finer patch beneath us
          continue;

        // Now that we have the correct interpolator, figure out the 
        // source and destination (rectangular) regions.
        int i1, i2, j1, j2, k1, k2; // FIXME
        int ii1 = 0, ii2 = ng;
        int jj1 = J * grid->py/ref_ratio, jj2 = (J+1) * grid->py/ref_ratio - 1;
        int kk1 = K * grid->pz/ref_ratio, kk2 = (K+1) * grid->pz/ref_ratio - 1;

        // Interpolate.
        amr_grid_interpolator_interpolate(interp, i1, i2, j1, j2, k1, k2,
                                          ii1, ii2, jj1, jj2, kk1, kk2, dest);
      }
    }
  }
  return 0;
}

static int local_fine_to_coarse_east_to_west(amr_grid_t* grid, 
                                             amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                             amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT((i_src == i_dest + 1) || ((i_dest == grid->nx-1) && (i_src == 0)));
  ASSERT(j_src == j_dest);
  ASSERT(k_src == k_dest);
  ASSERT(grid->finer != NULL);
  ASSERT(grid->ref_ratio > 0);
  ASSERT(dest->nc == src->nc);

  int nc = dest->nc, ng = grid->ng;
  int ref_ratio = grid->ref_ratio;

  // Loop over all the fine patches that correspond to this coarser one, and 
  // perform the interpolation.
  for (int I = 0; I < ref_ratio; ++I)
  {
    int i_fine = ref_ratio * i_src + I; 
    for (int J = 0; J < ref_ratio; ++J)
    {
      int j_fine = ref_ratio * j_src + J;
      for (int K = 0; K < ref_ratio; ++K)
      {
        int k_fine = ref_ratio * k_src;
        amr_grid_interpolator_t* interp;// = get_interpolator(grid->finer, i_fine, j_fine, k_fine);
        if (interp == NULL) // no finer patch beneath us
          continue;

        // Now that we have the correct interpolator, figure out the 
        // source and destination (rectangular) regions.
        int i1, i2, j1, j2, k1, k2; // FIXME
        int ii1 = dest->i2, ii2 = dest->i2 + ng;
        int jj1 = J * grid->py/ref_ratio, jj2 = (J+1) * grid->py/ref_ratio - 1;
        int kk1 = K * grid->pz/ref_ratio, kk2 = (K+1) * grid->pz/ref_ratio - 1;

        // Interpolate.
        amr_grid_interpolator_interpolate(interp, i1, i2, j1, j2, k1, k2,
                                          ii1, ii2, jj1, jj2, kk1, kk2, dest);
      }
    }
  }
  return 0;
}

static void local_fine_to_coarse_interpolation(amr_grid_t* grid, 
                                               amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                               amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT(grid->finer != NULL);
  ASSERT(grid->ref_ratio > 0);
  ASSERT(dest->nc == src->nc);

  int nc = dest->nc, ng = grid->ng;
  int ref_ratio = grid->ref_ratio;

  // Loop over all the fine patches that correspond to this coarser one, and 
  // perform the interpolation.
  for (int I = 0; I < ref_ratio; ++I)
  {
    int i_fine = ref_ratio * i_src + I; 
    for (int J = 0; J < ref_ratio; ++J)
    {
      int j_fine = ref_ratio * j_src + J;
      for (int K = 0; K < ref_ratio; ++K)
      {
        int k_fine = ref_ratio * k_src;
        amr_grid_interpolator_t* interp;// = get_interpolator(grid->finer, i_fine, j_fine, k_fine);
        if (interp == NULL) // no finer patch beneath us
          continue;

        // Now that we have the correct interpolator, figure out the 
        // source and destination (rectangular) regions.
        int i1, i2, j1, j2, k1, k2; // source cells
        int ii1, ii2, jj1, jj2, kk1, kk2; // destination cells
        if ((i_src == i_dest - 1) || ((i_dest == 0) && (i_src == grid->nx-1)))
        {
          ii1 = 0, ii2 = ng;
          jj1 = J * grid->py/ref_ratio, jj2 = (J+1) * grid->py/ref_ratio - 1;
          kk1 = K * grid->pz/ref_ratio, kk2 = (K+1) * grid->pz/ref_ratio - 1;
        }
        else if ((i_src == i_dest + 1) || ((i_dest == grid->nx-1) && (i_src == 0)))
        {
          ii1 = dest->i2, ii2 = dest->i2 + ng;
          jj1 = J * grid->py/ref_ratio, jj2 = (J+1) * grid->py/ref_ratio - 1;
          kk1 = K * grid->pz/ref_ratio, kk2 = (K+1) * grid->pz/ref_ratio - 1;
        }
        else if ((j_src == j_dest - 1) || ((j_dest == 0) && (j_src == grid->ny-1)))
        {
          ii1 = I * grid->px/ref_ratio, ii2 = (I+1) * grid->px/ref_ratio - 1;
          jj1 = 0, jj2 = ng;
          kk1 = K * grid->pz/ref_ratio, kk2 = (K+1) * grid->pz/ref_ratio - 1;
        }
        else if ((j_src == j_dest + 1) || ((j_dest == grid->ny-1) && (j_src == 0)))
        {
          ii1 = I * grid->px/ref_ratio, ii2 = (I+1) * grid->px/ref_ratio - 1;
          jj1 = dest->j2, jj2 = dest->j2 + ng;
          kk1 = K * grid->pz/ref_ratio, kk2 = (K+1) * grid->pz/ref_ratio - 1;
        }
        else if ((k_src == k_dest - 1) || ((k_dest == 0) && (k_src == grid->nz-1)))
        {
          ii1 = I * grid->px/ref_ratio, ii2 = (I+1) * grid->px/ref_ratio - 1;
          jj1 = J * grid->py/ref_ratio, jj2 = (J+1) * grid->py/ref_ratio - 1;
          kk1 = 0, kk2 = ng;
        }
        else
        {
          ii1 = I * grid->px/ref_ratio, ii2 = (I+1) * grid->px/ref_ratio - 1;
          jj1 = J * grid->py/ref_ratio, jj2 = (J+1) * grid->py/ref_ratio - 1;
          kk1 = dest->k2, kk2 = dest->k2 + ng;
        }

        // Interpolate.
        amr_grid_interpolator_interpolate(interp, i1, i2, j1, j2, k1, k2,
                                          ii1, ii2, jj1, jj2, kk1, kk2, dest);
      }
    }
  }
}

static void local_coarse_to_fine_interpolation(amr_grid_t* grid, 
                                               amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                               amr_patch_t* src,  int i_src, int j_src, int k_src)
{
}

static int start_remote_ghost_copy(amr_grid_t* grid, 
                                   amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                   amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  return -1;
}

static void finish_remote_ghost_copy(amr_grid_t* grid, int token,
                                     amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                     amr_patch_t* src,  int i_src, int j_src, int k_src)
{
}

static int start_remote_fine_to_coarse_interpolation(amr_grid_t* grid, 
                                                     amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                                     amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  return -1;
}

static void finish_remote_fine_to_coarse_interpolation(amr_grid_t* grid, int token)
{
}

static int start_remote_coarse_to_fine_interpolation(amr_grid_t* grid, 
                                                     amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                                     amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  return -1;
}

static void finish_remote_coarse_to_fine_interpolation(amr_grid_t* grid, int token)
{
}

void amr_grid_start_filling_ghosts(amr_grid_t* grid, amr_grid_data_t* data)
{
  ASSERT(amr_grid_data_grid(data) == grid);
  ASSERT(ptr_array_empty(grid->ghost_filler_finishers));

  // Go through our list of starters and tick them off.
  for (int i = 0; i < grid->ghost_filler_starters->size; ++i)
  {
    // Start the process of filling the ghost cells.
    amr_grid_ghost_filler_starter_t* starter = grid->ghost_filler_starters->data[i];
    amr_patch_t* dest_patch = amr_grid_data_patch(data, starter->i_dest, starter->j_dest, starter->k_dest);
    ASSERT(dest_patch != NULL);
    amr_patch_t* src_patch = amr_grid_data_patch(data, starter->i_src, starter->j_src, starter->k_src);
    int token = starter->method(grid, 
                                dest_patch, 
                                starter->i_dest, starter->j_dest, starter->k_dest,
                                src_patch, 
                                starter->i_src, starter->j_src, starter->k_src);
    if (token != -1)
    {
      // If this process has an asynchronous part, create a finisher for it.
      amr_grid_ghost_filler_finisher_t* finisher = polymec_malloc(sizeof(amr_grid_ghost_filler_finisher_t));
      finisher->token = token;
      finisher->dest_patch = dest_patch;
      finisher->i_dest = starter->i_dest;
      finisher->j_dest = starter->j_dest;
      finisher->k_dest = starter->k_dest;
      finisher->src_patch = src_patch;
      finisher->i_src = starter->i_src;
      finisher->j_src = starter->j_src;
      finisher->k_src = starter->k_src;
      ptr_array_append_with_dtor(grid->ghost_filler_finishers, finisher, polymec_free);
    }
  }
}

void amr_grid_finish_filling_ghosts(amr_grid_t* grid, amr_grid_data_t* data)
{
  // Go through our existing list of "finishers."
  for (int i = 0; i < grid->ghost_filler_finishers->size; ++i)
  {
    amr_grid_ghost_filler_finisher_t* finisher = grid->ghost_filler_finishers->data[i];
    finisher->method(grid, finisher->token, 
                     finisher->dest_patch, 
                     finisher->i_dest, finisher->j_dest, finisher->k_dest,
                     finisher->src_patch, 
                     finisher->i_src, finisher->j_src, finisher->k_src);
  }

  // Clear the list.
  ptr_array_clear(grid->ghost_filler_finishers);
}

void amr_grid_fill_ghosts(amr_grid_t* grid, amr_grid_data_t* data)
{
  amr_grid_start_filling_ghosts(grid, data);
  amr_grid_finish_filling_ghosts(grid, data);
}
