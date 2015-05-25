// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polyamri/amr_grid_level.h"
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

struct amr_grid_level_t
{
  bbox_t domain;
  int nx, ny, nz, tx, ty, tz, ng;
  patch_type_t* patch_types, num_patches;
  amr_grid_interpolator_t** interpolators;
  bool x_periodic, y_periodic, z_periodic;

  int ref_ratio;
  amr_grid_level_t* coarser;
  amr_grid_level_t* finer;
};

amr_grid_level_t* amr_grid_level_new(bbox_t* domain, 
                             int nx, int ny, int nz, 
                             int tx, int ty, int tz, 
                             int num_ghosts, 
                             bool periodic_in_x, bool periodic_in_y, bool periodic_in_z)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  ASSERT(tx > 0);
  ASSERT(ty > 0);
  ASSERT(tz > 0);
  ASSERT(num_ghosts >= 0);
  ASSERT(num_ghosts <= tx);
  ASSERT(num_ghosts <= ty);
  ASSERT(num_ghosts <= tz);

  amr_grid_level_t* level = polymec_malloc(sizeof(amr_grid_level_t));
  level->domain = *domain;
  level->nx = nx;
  level->ny = ny;
  level->nz = nz;
  level->tx = tx;
  level->ty = ty;
  level->tz = tz;
  level->ng = num_ghosts;
  level->x_periodic = periodic_in_x;
  level->y_periodic = periodic_in_y;
  level->z_periodic = periodic_in_z;
  level->patch_types = polymec_malloc(sizeof(patch_type_t) * nx * ny * nz);
  level->interpolators = polymec_malloc(sizeof(amr_grid_interpolator_t*) * nx * ny * nz);
  for (int i = 0; i < nx * ny * nz; ++i)
  {
    level->patch_types[i] = NONE;
    level->interpolators[i] = NULL;
  }

  level->ref_ratio = 0;
  level->coarser = level->finer = NULL;

  return level;
}

void amr_grid_level_free(amr_grid_level_t* level)
{
  for (int i = 0; i < level->num_patches; ++i)
  {
    if (level->interpolators[i] != NULL)
      amr_grid_interpolator_free(level->interpolators[i]);
  }
  polymec_free(level->patch_types);
  polymec_free(level->interpolators);
  polymec_free(level);
}

void amr_grid_level_associate_finer_level(amr_grid_level_t* level, amr_grid_level_t* finer_level, int ref_ratio)
{
  ASSERT(ref_ratio > 0);
  ASSERT((level->ref_ratio == 0) || (level->ref_ratio == ref_ratio));

  if (level->ref_ratio == 0)
    level->ref_ratio = ref_ratio;

  level->finer = finer_level;
}

void amr_grid_level_associate_coarser_level(amr_grid_level_t* level, amr_grid_level_t* coarser_level, int ref_ratio)
{
  ASSERT(ref_ratio > 0);
  ASSERT((level->ref_ratio == 0) || (level->ref_ratio == ref_ratio));

  if (level->ref_ratio == 0)
    level->ref_ratio = ref_ratio;

  level->coarser = coarser_level;
}

bbox_t* amr_grid_level_domain(amr_grid_level_t* level)
{
  return &level->domain;
}

void amr_grid_level_get_periodicity(amr_grid_level_t* level, bool* periodicity)
{
  ASSERT(periodicity != NULL);
  periodicity[0] = level->x_periodic;
  periodicity[1] = level->y_periodic;
  periodicity[2] = level->z_periodic;
}

void amr_grid_level_add_patch(amr_grid_level_t* level, int i, int j, int k)
{
  patch_type_t patch_type = level->patch_types[level->nz*level->ny*i + level->nz*j + k];
  ASSERT(patch_type == NONE);
  level->patch_types[level->nz*level->ny*i + level->nz*j + k] = LOCAL_SAME_LEVEL;
  ++(level->num_patches);
}

int amr_grid_level_num_patches(amr_grid_level_t* level)
{
  return level->num_patches;
}

static inline patch_type_t get_patch_type(amr_grid_level_t* level, int i, int j, int k)
{
  return level->patch_types[level->nz*level->ny*i + level->nz*j + k];
}

static inline amr_grid_interpolator_t* get_interpolator(amr_grid_level_t* level, int i, int j, int k)
{
  return level->interpolators[level->nz*level->ny*i + level->nz*j + k];
}

static inline bool local_patch_is_present(amr_grid_level_t* level, int i, int j, int k)
{
  return (get_patch_type(level, i, j, k) == LOCAL_SAME_LEVEL);
}

amr_patch_set_t* amr_grid_level_patch_set(amr_grid_level_t* level, int num_components)
{
  ASSERT(num_components > 0);

  amr_patch_set_t* patches = amr_patch_set_new();
  real_t dx = (level->domain.x2 - level->domain.x1) / level->nx;
  real_t dy = (level->domain.y2 - level->domain.y1) / level->ny;
  real_t dz = (level->domain.z2 - level->domain.z1) / level->nz;
  for (int i = 0; i < level->nx; ++i)
  {
    for (int j = 0; j < level->ny; ++j)
    {
      for (int k = 0; k < level->nz; ++k)
      {
        if (local_patch_is_present(level, i, j, k))
        {
          amr_patch_t* patch = amr_patch_new(level->tx, level->ty, level->tz, num_components, level->ng);
          bbox_t domain = {.x1 = level->domain.x1 + i*dx, 
                           .x2 = level->domain.x1 + (i+1)*dx,
                           .y1 = level->domain.y1 + j*dx, 
                           .y2 = level->domain.y1 + (j+1)*dy,
                           .z1 = level->domain.z1 + k*dz, 
                           .z2 = level->domain.z1 + (k+1)*dz};
          amr_patch_set_add(patches, patch, &domain);
        }
      }
    }
  }
  return patches;
}

static void local_ghost_copy(amr_grid_level_t* level, 
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

  int nc = src->nc, ng = level->ng;
  int i1 = src->i1, i2 = src->i2, 
      j1 = src->j1, j2 = src->j2,
      k1 = src->k1, k2 = src->k2;
  if ((i_src == i_dest - 1) || ((i_dest == 0) && (i_src == level->nx-1)))
  {
    for (int jj = j1; jj < j2; ++jj)
      for (int kk = k1; kk < k2; ++kk)
        for (int g = 0; g < ng; ++g)
          for (int c = 0; c < nc; ++c)
            d[i1-g-1][jj][kk][c] = s[i2-g-1][jj][kk][c];
  }
  else if ((i_src == i_dest + 1) || ((i_dest == level->nx-1) && (i_src == 0)))
  {
    for (int jj = j1; jj < j2; ++jj)
      for (int kk = k1; kk < k2; ++kk)
        for (int g = 0; g < ng; ++g)
          for (int c = 0; c < nc; ++c)
            d[i2+g][jj][kk][c] = s[i1+ng-g-1][jj][kk][c];
  }
  else if ((j_src == j_dest - 1) || ((j_dest == 0) && (j_src == level->ny-1)))
  {
    for (int ii = i1; ii < i2; ++ii)
      for (int kk = k1; kk < k2; ++kk)
        for (int g = 0; g < ng; ++g)
          for (int c = 0; c < nc; ++c)
            d[ii][j1-g-1][kk][c] = s[ii][j2-g-1][kk][c];
  }
  else if ((j_src == j_dest + 1) || ((j_dest == level->ny-1) && (j_src == 0)))
  {
    for (int ii = i1; ii < i2; ++ii)
      for (int kk = k1; kk < k2; ++kk)
        for (int g = 0; g < ng; ++g)
          for (int c = 0; c < nc; ++c)
            d[ii][j2+g][kk][c] = s[ii][j1+ng-g-1][kk][c];
  }
  else if ((k_src == k_dest - 1) || ((k_dest == 0) && (k_src == level->nz-1)))
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

static void local_fine_to_coarse_interpolation(amr_grid_level_t* level, 
                                               amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                               amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT(level->finer != NULL);
  ASSERT(level->ref_ratio > 0);
  ASSERT(dest->nc == src->nc);

  int nc = dest->nc, ng = level->ng;
  int ref_ratio = level->ref_ratio;

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
        amr_grid_interpolator_t* interp = get_interpolator(level->finer, i_fine, j_fine, k_fine);
        if (interp == NULL) // no finer patch beneath us
          continue;

        // Now that we have the correct interpolator, figure out the 
        // source and destination (rectangular) regions.
        int i1, i2, j1, j2, k1, k2; // source cells
        int ii1, ii2, jj1, jj2, kk1, kk2; // destination cells
        if ((i_src == i_dest - 1) || ((i_dest == 0) && (i_src == level->nx-1)))
        {
          ii1 = 0, ii2 = ng;
          jj1 = J * level->ty/ref_ratio, jj2 = (J+1) * level->ty/ref_ratio - 1;
          kk1 = K * level->tz/ref_ratio, kk2 = (K+1) * level->tz/ref_ratio - 1;
        }
        else if ((i_src == i_dest + 1) || ((i_dest == level->nx-1) && (i_src == 0)))
        {
          ii1 = dest->i2, ii2 = dest->i2 + ng;
          jj1 = J * level->ty/ref_ratio, jj2 = (J+1) * level->ty/ref_ratio - 1;
          kk1 = K * level->tz/ref_ratio, kk2 = (K+1) * level->tz/ref_ratio - 1;
        }
        else if ((j_src == j_dest - 1) || ((j_dest == 0) && (j_src == level->ny-1)))
        {
          ii1 = I * level->tx/ref_ratio, ii2 = (I+1) * level->tx/ref_ratio - 1;
          jj1 = 0, jj2 = ng;
          kk1 = K * level->tz/ref_ratio, kk2 = (K+1) * level->tz/ref_ratio - 1;
        }
        else if ((j_src == j_dest + 1) || ((j_dest == level->ny-1) && (j_src == 0)))
        {
          ii1 = I * level->tx/ref_ratio, ii2 = (I+1) * level->tx/ref_ratio - 1;
          jj1 = dest->j2, jj2 = dest->j2 + ng;
          kk1 = K * level->tz/ref_ratio, kk2 = (K+1) * level->tz/ref_ratio - 1;
        }
        else if ((k_src == k_dest - 1) || ((k_dest == 0) && (k_src == level->nz-1)))
        {
          ii1 = I * level->tx/ref_ratio, ii2 = (I+1) * level->tx/ref_ratio - 1;
          jj1 = J * level->ty/ref_ratio, jj2 = (J+1) * level->ty/ref_ratio - 1;
          kk1 = 0, kk2 = ng;
        }
        else
        {
          ii1 = I * level->tx/ref_ratio, ii2 = (I+1) * level->tx/ref_ratio - 1;
          jj1 = J * level->ty/ref_ratio, jj2 = (J+1) * level->ty/ref_ratio - 1;
          kk1 = dest->k2, kk2 = dest->k2 + ng;
        }

        // Interpolate.
        amr_grid_interpolator_interpolate(interp, i1, i2, j1, j2, k1, k2,
                                          ii1, ii2, jj1, jj2, kk1, kk2, dest);
      }
    }
  }
}

static void local_coarse_to_fine_interpolation(amr_grid_level_t* level, 
                                               amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                               amr_patch_t* src,  int i_src, int j_src, int k_src)
{
}

static int start_remote_ghost_copy(amr_grid_level_t* level, 
                                   amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                   amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  return -1;
}

static void finish_remote_ghost_copy(amr_grid_level_t* level, int token)
{
}

static int start_remote_fine_to_coarse_interpolation(amr_grid_level_t* level, 
                                                     amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                                     amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  return -1;
}

static void finish_remote_fine_to_coarse_interpolation(amr_grid_level_t* level, int token)
{
}

static int start_remote_coarse_to_fine_interpolation(amr_grid_level_t* level, 
                                                     amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                                     amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  return -1;
}

static void finish_remote_coarse_to_fine_interpolation(amr_grid_level_t* level, int token)
{
}

static int start_filling_ghosts(amr_grid_level_t* level, 
                                amr_patch_t* dest, int i_dest, int j_dest, int k_dest, 
                                amr_patch_t* src,  int i_src, int j_src, int k_src)
{
  ASSERT(dest != NULL);
  ASSERT(dest != src);
  ASSERT(((i_dest == i_src-1) || (i_dest == i_src+1) || 
          (j_dest == j_src-1) || (j_dest == j_src+1) || 
          (k_dest == k_src-1) || (k_dest == k_src+1)) ||
          (level->x_periodic && (i_src == 0) && (i_dest == level->nx-1)) ||
          (level->x_periodic && (i_src == level->nx-1) && (i_dest == 0)) ||
          (level->y_periodic && (j_src == 0) && (j_dest == level->ny-1)) ||
          (level->y_periodic && (j_src == level->ny-1) && (j_dest == 0)) ||
          (level->z_periodic && (k_src == 0) && (k_dest == level->nz-1)) ||
          (level->z_periodic && (k_src == level->nz-1) && (k_dest == 0)));

  int token = -1;
  patch_type_t src_patch_type = get_patch_type(level, i_src, j_src, k_src);
  switch(src_patch_type)
  {
    case LOCAL_SAME_LEVEL:
      ASSERT(src != NULL);
      local_ghost_copy(level, dest, i_dest, j_dest, k_dest, 
                       src, i_src, j_src, k_src);
      break;
    case LOCAL_FINER_LEVEL:
      local_fine_to_coarse_interpolation(level, dest, i_dest, j_dest, k_dest, 
                                         src, i_src, j_src, k_src);
      break;
    case LOCAL_COARSER_LEVEL:
      local_coarse_to_fine_interpolation(level, dest, i_dest, j_dest, k_dest, 
                                         src, i_src, j_src, k_src);
      break;
    case REMOTE_SAME_LEVEL:
      token = start_remote_ghost_copy(level, dest, i_dest, j_dest, k_dest, 
                                      src, i_src, j_src, k_src);
      break;
    case REMOTE_FINER_LEVEL:
      token = start_remote_fine_to_coarse_interpolation(level, dest, i_dest, j_dest, k_dest, 
                                                        src, i_src, j_src, k_src);
      break;
    case REMOTE_COARSER_LEVEL:
      token = start_remote_coarse_to_fine_interpolation(level, dest, i_dest, j_dest, k_dest, 
                                                        src, i_src, j_src, k_src);
    default:
      break;
  }

  // FIXME: Stash the token somewhere we can get to this transaction.

  return token;
}

void amr_grid_level_start_filling_ghosts(amr_grid_level_t* level, amr_patch_set_t* patches)
{
  if (level->ng == 0)
    return;

  // If we made this patch set, we know how it's indexed. Get pointers 
  // to each of the local patches within.
  int pos = 0;
  amr_patch_t* patch_ptrs[level->nx][level->ny][level->nz];
  for (int i = 0; i < level->nx; ++i)
  {
    for (int j = 0; j < level->ny; ++j)
    {
      for (int k = 0; k < level->nz; ++k)
      {
        bbox_t* domain = NULL;
        patch_ptrs[i][j][k] = NULL;
        if (local_patch_is_present(level, i, j, k))
          amr_patch_set_next(patches, &pos, &patch_ptrs[i][j][k], &domain);
      }
    }
  }

  // Now perform all local copies.
  for (int i = 0; i < level->nx; ++i)
  {
    for (int j = 0; j < level->ny; ++j)
    {
      for (int k = 0; k < level->nz; ++k)
      {
        amr_patch_t* dest = patch_ptrs[i][j][k];
        if (dest == NULL) continue;

        if (i > 0)
        {
          amr_patch_t* src = patch_ptrs[i-1][j][k];
          int token = start_filling_ghosts(level, dest, i, j, k, src, i-1, j, k);
        }
        else if (level->x_periodic)
        {
          amr_patch_t* src = patch_ptrs[level->nx-1][j][k];
          int token = start_filling_ghosts(level, dest, i, j, k, src, level->nx-1, j, k);
        }

        if (i < level->nx-1)
        {
          amr_patch_t* src = patch_ptrs[i+1][j][k];
          int token = start_filling_ghosts(level, dest, i, j, k, src, i+1, j, k);
        }
        else if (level->x_periodic)
        {
          amr_patch_t* src = patch_ptrs[0][j][k];
          int token = start_filling_ghosts(level, dest, i, j, k, src, 0, j, k);
        }

        if (j > 0)
        {
          amr_patch_t* src = patch_ptrs[i][j-1][k];
          int token = start_filling_ghosts(level, dest, i, j, k, src, i, j-1, k);
        }
        else if (level->y_periodic)
        {
          amr_patch_t* src = patch_ptrs[i][level->ny-1][k];
          int token = start_filling_ghosts(level, dest, i, j, k, src, i, level->ny-1, k);
        }

        if (j < level->ny-1)
        {
          amr_patch_t* src = patch_ptrs[i][j+1][k];
          int token = start_filling_ghosts(level, dest, i, j, k, src, i, j+1, k);
        }
        else if (level->y_periodic)
        {
          amr_patch_t* src = patch_ptrs[i][0][k];
          int token = start_filling_ghosts(level, dest, i, j, k, src, i, 0, k);
        }

        if (k > 0)
        {
          amr_patch_t* src = patch_ptrs[i][j][k-1];
          int token = start_filling_ghosts(level, dest, i, j, k, src, i, j, k-1);
        }
        else if (level->z_periodic)
        {
          amr_patch_t* src = patch_ptrs[i][j][level->nz-1];
          int token = start_filling_ghosts(level, dest, i, j, k, src, i, j, level->nz-1);
        }

        if (k < level->ny-1)
        {
          amr_patch_t* src = patch_ptrs[i][j][k+1];
          int token = start_filling_ghosts(level, dest, i, j, k, src, i, j, k+1);
        }
        else if (level->z_periodic)
        {
          amr_patch_t* src = patch_ptrs[i][j][0];
          int token = start_filling_ghosts(level, dest, i, j, k, src, i, j, 0);
        }
      }
    }
  }
}

void amr_grid_level_finish_filling_ghosts(amr_grid_level_t* level, amr_patch_set_t* patches)
{
#if 0
  switch(src_patch_type)
  {
    case REMOTE_SAME_LEVEL:
      finish_remote_ghost_copy(level, token);
      break;
    case REMOTE_FINER_LEVEL:
      finish_remote_fine_to_coarse_interpolation(level, token);
      break;
    case REMOTE_COARSER_LEVEL:
      finish_remote_coarse_to_fine_interpolation(level, token);
    default:
      break;
  }
#endif
}

void amr_grid_level_fill_ghosts(amr_grid_level_t* level, amr_patch_set_t* patches)
{
  amr_grid_level_start_filling_ghosts(level, patches);
  amr_grid_level_finish_filling_ghosts(level, patches);
}
