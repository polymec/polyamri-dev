// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "polyamri/amr_grid.h"
#include "polyamri/amr_grid_interpolator.h"

// This is a descriptor for "types" of patches -- this determines the nature 
// of interacting patches during ghost cell filling.
typedef enum
{
  NONE,
  LOCAL_SAME_LEVEL,
  LOCAL_FINER_LEVEL,
  LOCAL_COARSER_LEVEL,
  REMOTE_SAME_LEVEL,
  REMOTE_FINER_LEVEL,
  REMOTE_COARSER_LEVEL
} patch_type_t;

struct amr_grid_t
{
  bbox_t domain;
  int nx, ny, nz, px, py, pz, ng, num_patches;
  patch_type_t* patch_types;
  int* remote_owners;
  amr_grid_interpolator_t** interpolators;
  bool x_periodic, y_periodic, z_periodic;

  // Relationships with other grids.
  int ref_ratio;
  amr_grid_t* coarser;
  amr_grid_t* finer;

  // Neighboring grids.
  amr_grid_t* neighbors[6];

  // Associated data.
  int_ptr_unordered_map_t* data;
};

amr_grid_t* amr_grid_new(bbox_t* domain, 
                         int nx, int ny, int nz, 
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
  grid->domain = *domain;
  grid->nx = nx;
  grid->ny = ny;
  grid->nz = nz;
  grid->px = px;
  grid->py = py;
  grid->pz = pz;
  grid->ng = num_ghosts;
  grid->num_patches = 0;
  grid->x_periodic = periodic_in_x;
  grid->y_periodic = periodic_in_y;
  grid->z_periodic = periodic_in_z;
  grid->patch_types = polymec_malloc(sizeof(patch_type_t) * nx * ny * nz);
  grid->remote_owners = polymec_malloc(sizeof(int) * nx * ny * nz);
  grid->interpolators = polymec_malloc(sizeof(amr_grid_interpolator_t*) * nx * ny * nz);
  for (int i = 0; i < nx * ny * nz; ++i)
  {
    grid->patch_types[i] = NONE;
    grid->remote_owners[i] = -1;
    grid->interpolators[i] = NULL;
  }

  grid->ref_ratio = 0;
  grid->coarser = grid->finer = NULL;

  grid->data = int_ptr_unordered_map_new();

  return grid;
}

void amr_grid_free(amr_grid_t* grid)
{
  for (int i = 0; i < grid->num_patches; ++i)
  {
    if (grid->interpolators[i] != NULL)
      amr_grid_interpolator_free(grid->interpolators[i]);
  }
  polymec_free(grid->patch_types);
  polymec_free(grid->interpolators);
  polymec_free(grid->remote_owners);
  polymec_free(grid->data);
  polymec_free(grid);
}

void amr_grid_set_neighbor(amr_grid_t* grid, 
                           amr_grid_neighbor_slot_t neighbor_slot,
                           amr_grid_t* neighbor)
{
  grid->neighbors[neighbor_slot] = neighbor;
}

void amr_grid_set_data(amr_grid_t* grid,
                       int data_index,
                       void* data,
                       void (*data_dtor)(void* data))
{
  int_ptr_unordered_map_insert_with_v_dtor(grid->data, data_index, data, data_dtor);
}

void* amr_grid_data(amr_grid_t* grid,
                    int data_index)
{
  void** ptr = int_ptr_unordered_map_get(grid->data, data_index);
  if (ptr != NULL)
    return *ptr;
  else
    return NULL;
}

void amr_grid_associate_finer(amr_grid_t* grid, amr_grid_t* finer_grid, int ref_ratio)
{
  ASSERT(ref_ratio > 0);
  ASSERT((grid->ref_ratio == 0) || (grid->ref_ratio == ref_ratio));

  if (grid->ref_ratio == 0)
    grid->ref_ratio = ref_ratio;

  grid->finer = finer_grid;
}

void amr_grid_associate_coarser(amr_grid_t* grid, amr_grid_t* coarser_grid, int ref_ratio)
{
  ASSERT(ref_ratio > 0);
  ASSERT((grid->ref_ratio == 0) || (grid->ref_ratio == ref_ratio));

  if (grid->ref_ratio == 0)
    grid->ref_ratio = ref_ratio;

  grid->coarser = coarser_grid;
}

bbox_t* amr_grid_domain(amr_grid_t* grid)
{
  return &grid->domain;
}

void amr_grid_get_periodicity(amr_grid_t* grid, bool* periodicity)
{
  ASSERT(periodicity != NULL);
  periodicity[0] = grid->x_periodic;
  periodicity[1] = grid->y_periodic;
  periodicity[2] = grid->z_periodic;
}

void amr_grid_add_local_patch(amr_grid_t* grid, int i, int j, int k)
{
  int patch_index = grid->nz*grid->ny*i + grid->nz*j + k;
  patch_type_t patch_type = grid->patch_types[patch_index];
  ASSERT(patch_type == NONE);
  grid->patch_types[patch_index] = LOCAL_SAME_LEVEL;
  ++(grid->num_patches);
}

void amr_grid_add_remote_patch(amr_grid_t* grid, int i, int j, int k, int remote_owner)
{
  int patch_index = grid->nz*grid->ny*i + grid->nz*j + k;
  patch_type_t patch_type = grid->patch_types[patch_index];
  ASSERT(patch_type == NONE);
  grid->patch_types[patch_index] = REMOTE_SAME_LEVEL;
  grid->remote_owners[patch_index] = remote_owner;
  ++(grid->num_patches);
}

int amr_grid_num_patches(amr_grid_t* grid)
{
  return grid->num_patches;
}

static inline patch_type_t get_patch_type(amr_grid_t* grid, int i, int j, int k)
{
  return grid->patch_types[grid->nz*grid->ny*i + grid->nz*j + k];
}

static inline amr_grid_interpolator_t* get_interpolator(amr_grid_t* grid, int i, int j, int k)
{
  return grid->interpolators[grid->nz*grid->ny*i + grid->nz*j + k];
}

static inline bool local_patch_is_present(amr_grid_t* grid, int i, int j, int k)
{
  return (get_patch_type(grid, i, j, k) == LOCAL_SAME_LEVEL);
}

amr_patch_set_t* amr_grid_create_patches(amr_grid_t* grid, int num_components)
{
  ASSERT(num_components > 0);

  amr_patch_set_t* patches = amr_patch_set_new();
  real_t dx = (grid->domain.x2 - grid->domain.x1) / grid->nx;
  real_t dy = (grid->domain.y2 - grid->domain.y1) / grid->ny;
  real_t dz = (grid->domain.z2 - grid->domain.z1) / grid->nz;
  for (int i = 0; i < grid->nx; ++i)
  {
    for (int j = 0; j < grid->ny; ++j)
    {
      for (int k = 0; k < grid->nz; ++k)
      {
        if (local_patch_is_present(grid, i, j, k))
        {
          amr_patch_t* patch = amr_patch_new(grid->px, grid->py, grid->pz, num_components, grid->ng);
          bbox_t domain = {.x1 = grid->domain.x1 + i*dx, 
                           .x2 = grid->domain.x1 + (i+1)*dx,
                           .y1 = grid->domain.y1 + j*dx, 
                           .y2 = grid->domain.y1 + (j+1)*dy,
                           .z1 = grid->domain.z1 + k*dz, 
                           .z2 = grid->domain.z1 + (k+1)*dz};
          amr_patch_set_add(patches, patch, &domain);
        }
      }
    }
  }
  return patches;
}

static void local_ghost_copy(amr_grid_t* grid, 
                             amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                             amr_patch_t* src,  int i_src, int j_src, int k_src)
{
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
  if ((i_src == i_dest - 1) || ((i_dest == 0) && (i_src == grid->nx-1)))
  {
    for (int jj = j1; jj < j2; ++jj)
      for (int kk = k1; kk < k2; ++kk)
        for (int g = 0; g < ng; ++g)
          for (int c = 0; c < nc; ++c)
            d[i1-g-1][jj][kk][c] = s[i2-g-1][jj][kk][c];
  }
  else if ((i_src == i_dest + 1) || ((i_dest == grid->nx-1) && (i_src == 0)))
  {
    for (int jj = j1; jj < j2; ++jj)
      for (int kk = k1; kk < k2; ++kk)
        for (int g = 0; g < ng; ++g)
          for (int c = 0; c < nc; ++c)
            d[i2+g][jj][kk][c] = s[i1+ng-g-1][jj][kk][c];
  }
  else if ((j_src == j_dest - 1) || ((j_dest == 0) && (j_src == grid->ny-1)))
  {
    for (int ii = i1; ii < i2; ++ii)
      for (int kk = k1; kk < k2; ++kk)
        for (int g = 0; g < ng; ++g)
          for (int c = 0; c < nc; ++c)
            d[ii][j1-g-1][kk][c] = s[ii][j2-g-1][kk][c];
  }
  else if ((j_src == j_dest + 1) || ((j_dest == grid->ny-1) && (j_src == 0)))
  {
    for (int ii = i1; ii < i2; ++ii)
      for (int kk = k1; kk < k2; ++kk)
        for (int g = 0; g < ng; ++g)
          for (int c = 0; c < nc; ++c)
            d[ii][j2+g][kk][c] = s[ii][j1+ng-g-1][kk][c];
  }
  else if ((k_src == k_dest - 1) || ((k_dest == 0) && (k_src == grid->nz-1)))
  {
    for (int ii = i1; ii < i2; ++ii)
      for (int jj = j1; jj < j2; ++jj)
        for (int g = 0; g < ng; ++g)
          for (int c = 0; c < nc; ++c)
            d[ii][jj][k1-g-1][c] = s[ii][jj][k2-g-1][c];
  }
  else
  {
    for (int ii = i1; ii < i2; ++ii)
      for (int jj = j1; jj < j2; ++jj)
        for (int g = 0; g < ng; ++g)
          for (int c = 0; c < nc; ++c)
            d[ii][jj][k2+g][c] = s[ii][jj][k1+ng-g-1][c];
  }
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
        amr_grid_interpolator_t* interp = get_interpolator(grid->finer, i_fine, j_fine, k_fine);
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

static void finish_remote_ghost_copy(amr_grid_t* grid, int token)
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

static int start_filling_ghosts(amr_grid_t* grid, 
                                amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT(dest != NULL);
  ASSERT(dest != src);
  ASSERT(((i_dest == i_src-1) || (i_dest == i_src+1) || 
          (j_dest == j_src-1) || (j_dest == j_src+1) || 
          (k_dest == k_src-1) || (k_dest == k_src+1)) ||
          (grid->x_periodic && (i_src == 0) && (i_dest == grid->nx-1)) ||
          (grid->x_periodic && (i_src == grid->nx-1) && (i_dest == 0)) ||
          (grid->y_periodic && (j_src == 0) && (j_dest == grid->ny-1)) ||
          (grid->y_periodic && (j_src == grid->ny-1) && (j_dest == 0)) ||
          (grid->z_periodic && (k_src == 0) && (k_dest == grid->nz-1)) ||
          (grid->z_periodic && (k_src == grid->nz-1) && (k_dest == 0)));

  int token = -1;
  patch_type_t src_patch_type = get_patch_type(grid, i_src, j_src, k_src);
  switch(src_patch_type)
  {
    case LOCAL_SAME_LEVEL:
      ASSERT(src != NULL);
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

  // FIXME: Stash the token somewhere we can get to this transaction.

  return token;
}

void amr_grid_start_filling_ghosts(amr_grid_t* grid, amr_patch_set_t* patches)
{
  if (grid->ng == 0)
    return;

  // If we made this patch set, we know how it's indexed. Get pointers 
  // to each of the local patches within.
  int pos = 0;
  amr_patch_t* patch_ptrs[grid->nx][grid->ny][grid->nz];
  for (int i = 0; i < grid->nx; ++i)
  {
    for (int j = 0; j < grid->ny; ++j)
    {
      for (int k = 0; k < grid->nz; ++k)
      {
        bbox_t* domain = NULL;
        patch_ptrs[i][j][k] = NULL;
        if (local_patch_is_present(grid, i, j, k))
          amr_patch_set_next(patches, &pos, &patch_ptrs[i][j][k], &domain);
      }
    }
  }

  // Now perform all local copies.
  for (int i = 0; i < grid->nx; ++i)
  {
    for (int j = 0; j < grid->ny; ++j)
    {
      for (int k = 0; k < grid->nz; ++k)
      {
        amr_patch_t* dest = patch_ptrs[i][j][k];
        if (dest == NULL) continue;

        if (i > 0)
        {
          amr_patch_t* src = patch_ptrs[i-1][j][k];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, i-1, j, k);
        }
        else if (grid->x_periodic)
        {
          amr_patch_t* src = patch_ptrs[grid->nx-1][j][k];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, grid->nx-1, j, k);
        }

        if (i < grid->nx-1)
        {
          amr_patch_t* src = patch_ptrs[i+1][j][k];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, i+1, j, k);
        }
        else if (grid->x_periodic)
        {
          amr_patch_t* src = patch_ptrs[0][j][k];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, 0, j, k);
        }

        if (j > 0)
        {
          amr_patch_t* src = patch_ptrs[i][j-1][k];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, i, j-1, k);
        }
        else if (grid->y_periodic)
        {
          amr_patch_t* src = patch_ptrs[i][grid->ny-1][k];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, i, grid->ny-1, k);
        }

        if (j < grid->ny-1)
        {
          amr_patch_t* src = patch_ptrs[i][j+1][k];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, i, j+1, k);
        }
        else if (grid->y_periodic)
        {
          amr_patch_t* src = patch_ptrs[i][0][k];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, i, 0, k);
        }

        if (k > 0)
        {
          amr_patch_t* src = patch_ptrs[i][j][k-1];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, i, j, k-1);
        }
        else if (grid->z_periodic)
        {
          amr_patch_t* src = patch_ptrs[i][j][grid->nz-1];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, i, j, grid->nz-1);
        }

        if (k < grid->ny-1)
        {
          amr_patch_t* src = patch_ptrs[i][j][k+1];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, i, j, k+1);
        }
        else if (grid->z_periodic)
        {
          amr_patch_t* src = patch_ptrs[i][j][0];
          int token = start_filling_ghosts(grid, dest, i, j, k, src, i, j, 0);
        }
      }
    }
  }
}

void amr_grid_finish_filling_ghosts(amr_grid_t* grid, amr_patch_set_t* patches)
{
#if 0
  switch(src_patch_type)
  {
    case REMOTE_SAME_LEVEL:
      finish_remote_ghost_copy(grid, token);
      break;
    case REMOTE_FINER_LEVEL:
      finish_remote_fine_to_coarse_interpolation(grid, token);
      break;
    case REMOTE_COARSER_LEVEL:
      finish_remote_coarse_to_fine_interpolation(grid, token);
    default:
      break;
  }
#endif
}

void amr_grid_fill_ghosts(amr_grid_t* grid, amr_patch_set_t* patches)
{
  amr_grid_start_filling_ghosts(grid, patches);
  amr_grid_finish_filling_ghosts(grid, patches);
}
