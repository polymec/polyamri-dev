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

// This stores a function and its arguments for performing the work of filling ghost cells.
typedef struct
{
  void* data;
  int (*method)(amr_grid_t*, amr_patch_t*, int, int, int, amr_patch_t*, int, int, int);
  amr_patch_t* dest_patch;
  int i_dest, j_dest, k_dest;
  amr_patch_t* src_patch;
  int i_src, j_src, k_src;
} ghost_filler_t;

// This is a generic buffer for storing serialized patch data.
typedef struct
{
  int* counts; // numbers of elements sent to/received from other procs.
  int* displs; // index displacements for data sent/received.
  MPI_Datatype type; // type of data sent/received.
  void* buffer; // actual data sent/received.
} patch_buffer_t;

static void patch_buffer_free(patch_buffer_t* buffer)
{
}

static void patch_buffer_copy_in(patch_buffer_t* buffer, amr_grid_data_t* grid_data)
{
}

static void patch_buffer_transform(patch_buffer_t* buffer)
{
}

static void patch_buffer_copy_out(patch_buffer_t* buffer, amr_grid_data_t* grid_data)
{
}

// This stores buffers and pointers to grid data for pending remote operations.
typedef struct
{
  patch_buffer_t* buffer;
  amr_grid_data_t* grid_data;
} remote_data_t;

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

  // Ratio of refinement between this grid and its "neighbors"
  int ref_ratio; 

  // MPI communicator.
  MPI_Comm comm;

  // Coarser grid and associated interpolator.
  amr_grid_t* coarser;
  amr_grid_interpolator_t* coarse_interpolator;

  // Finer grid and associated interpolator.
  amr_grid_t* finer;
  amr_grid_interpolator_t* fine_interpolator;

  // Neighboring grids and interpolators.
  amr_grid_t* neighbors[6];
  amr_grid_interpolator_t* neighbor_interpolators[6];

  // This flag is set by amr_grid_finalize() after a grid has been assembled.
  bool finalized;

  // Local patch buffer for assisting in local copies.
  ptr_array_t* local_buffers;

  // Data that hasn't been filled in remotely yet.
  int_ptr_unordered_map_t* pending_data;
};

static patch_buffer_t* patch_buffer_new(amr_grid_t* grid, int num_components)
{
  ASSERT(num_components > 0);
  return NULL;
}

static remote_data_t* remote_data_new(amr_grid_t* grid, amr_grid_data_t* grid_data)
{
  remote_data_t* data = polymec_malloc(sizeof(remote_data_t));
  data->buffer = patch_buffer_new(grid, amr_grid_data_num_components(grid_data));
  data->grid_data = grid_data;
  return data;
}

static void remote_data_free(remote_data_t* remote_data)
{
  patch_buffer_free(remote_data->buffer);
  polymec_free(remote_data);
}

amr_grid_t* amr_grid_new(MPI_Comm comm,
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
  grid->nx = nx;
  grid->ny = ny;
  grid->nz = nz;
  grid->px = px;
  grid->py = py;
  grid->pz = pz;
  grid->ng = num_ghosts;
  grid->x_periodic = periodic_in_x;
  grid->y_periodic = periodic_in_y;
  grid->z_periodic = periodic_in_z;
  grid->num_local_patches = 0;
  grid->local_patch_indices = NULL;
  grid->patch_types = polymec_malloc(sizeof(patch_type_t) * nx * ny * nz);
  grid->remote_owners = polymec_malloc(sizeof(int) * nx * ny * nz);
  grid->comm = comm;
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

  grid->finalized = false;
  grid->local_buffers = ptr_array_new();
  grid->pending_data = int_ptr_unordered_map_new();

  return grid;
}

void amr_grid_free(amr_grid_t* grid)
{
  int_ptr_unordered_map_free(grid->pending_data);
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
  if (grid->local_patch_indices != NULL)
    polymec_free(grid->local_patch_indices);
  ptr_array_free(grid->local_buffers);
  polymec_free(grid);
}

MPI_Comm amr_grid_comm(amr_grid_t* grid)
{
  return grid->comm;
}

void amr_grid_set_comm(amr_grid_t* grid, MPI_Comm comm)
{
  grid->comm = comm;
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
  DECLARE_3D_ARRAY(patch_type_t, patch_types, grid->patch_types, grid->nx, grid->ny, grid->nz);
  patch_type_t patch_type = patch_types[i][j][k];
  ASSERT(patch_type == NO_PATCH);
  patch_types[i][j][k] = LOCAL_SAME_LEVEL;
  ++(grid->num_local_patches);
}

void amr_grid_add_remote_patch(amr_grid_t* grid, int i, int j, int k, int remote_owner)
{
  ASSERT(!grid->finalized);
  DECLARE_3D_ARRAY(patch_type_t, patch_types, grid->patch_types, grid->nx, grid->ny, grid->nz);
  patch_type_t patch_type = patch_types[i][j][k];
  ASSERT(patch_type == NO_PATCH);
  patch_types[i][j][k] = REMOTE_SAME_LEVEL;
  DECLARE_3D_ARRAY(int, remote_owners, grid->remote_owners, grid->nx, grid->ny, grid->nz);
  remote_owners[i][j][k] = remote_owner;
}

void amr_grid_finalize(amr_grid_t* grid)
{
  ASSERT(!grid->finalized);

  // Make sure that we have accounted for all patches in our construction 
  // process.
  DECLARE_3D_ARRAY(patch_type_t, patch_types, grid->patch_types, grid->nx, grid->ny, grid->nz);
  if (grid->num_local_patches < (grid->nx * grid->ny * grid->nz))
  {
    for (int i = 0; i < grid->nx; ++i)
    {
      for (int j = 0; j < grid->ny; ++j)
      {
        for (int k = 0; k < grid->nz; ++k)
        {
          if (patch_types[i][j][k] == NO_PATCH)
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

  // Establish our local and remote copying patterns.
#if 0
  for (int i = 0; i < grid->nx; ++i)
  {
    for (int j = 0; j < grid->ny; ++j)
    {
      for (int k = 0; k < grid->nz; ++k)
      {
        patch_type_t patch_type = patch_types[i][j][k];
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

bool amr_grid_has_local_patch(amr_grid_t* grid, int i, int j, int k)
{
  DECLARE_3D_ARRAY(patch_type_t, patch_types, grid->patch_types, grid->nx, grid->ny, grid->nz);
  return (patch_types[i][j][k] == LOCAL_SAME_LEVEL);
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
  for (int j = j1; j < j2; ++j)
    for (int k = k1; k < k2; ++k) 
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[i1-g-1][j][k][c] = s[i2-g-1][j][k][c];
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
  for (int j = j1; j < j2; ++j)
    for (int k = k1; k < k2; ++k) 
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[i2+g][j][k][c] = s[i1+ng-g-1][j][k][c];
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
  for (int i = i1; i < i2; ++i)
    for (int k = k1; k < k2; ++k)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[i][j1-g-1][k][c] = s[i][j2-g-1][k][c];

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
  for (int i = i1; i < i2; ++i)
    for (int k = k1; k < k2; ++k)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[i][j2+g][k][c] = s[i][j1+ng-g-1][k][c];

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
  for (int i = i1; i < i2; ++i)
    for (int j = j1; j < j2; ++j)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[i][j][k1-g-1][c] = s[i][j][k2-g-1][c];

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
  for (int i = i1; i < i2; ++i)
    for (int j = j1; j < j2; ++j)
      for (int g = 0; g < ng; ++g)
        for (int c = 0; c < nc; ++c)
          d[i][j][k2+g][c] = s[i][j][k1+ng-g-1][c];

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

static void do_local_copies(amr_grid_t* grid, amr_grid_data_t* data)
{
  // Make sure we have a local patch buffer to use with the right number of components.
  int num_comp = amr_grid_data_num_components(data);
  if (num_comp > grid->local_buffers->size)
  {
    int old_size = grid->local_buffers->size;
    ptr_array_resize(grid->local_buffers, num_comp);
    for (int i = old_size; i < num_comp; ++i)
      grid->local_buffers->data[i] = NULL;
  }
  if (grid->local_buffers->data[num_comp-1] == NULL)
  {
    patch_buffer_t* new_buffer = patch_buffer_new(grid, num_comp);
    ptr_array_assign_with_dtor(grid->local_buffers, num_comp-1, new_buffer, DTOR(patch_buffer_free));
  }
  patch_buffer_t* buffer = grid->local_buffers->data[num_comp-1];

  // Copy the data into our local patch buffer and then out again.
  patch_buffer_copy_in(buffer, data);

  // Copy and/or interpolate locally.
  patch_buffer_transform(buffer);

  // Copy out.
  patch_buffer_copy_out(buffer, data);
}

static int start_remote_copies(amr_grid_t* grid, amr_grid_data_t* data)
{
  // Select the next available (non-negative) token.
  int token = 0;
  while (int_ptr_unordered_map_contains(grid->pending_data, token))
    ++token;

  // Make a new remote data thingy.
  int_ptr_unordered_map_insert_with_v_dtor(grid->pending_data, token, data, DTOR(remote_data_free));
  remote_data_t* remote_data = remote_data_new(grid, data);

  // Copy the data into the thingy's patch buffer.
  patch_buffer_copy_in(remote_data->buffer, remote_data->grid_data);

  // Initiate the data transfer.
  // FIXME

  return token;
}

static void finish_remote_copies(amr_grid_t* grid, int token)
{
  remote_data_t** remote_data_p = (remote_data_t**)int_ptr_unordered_map_get(grid->pending_data, token);
  ASSERT(remote_data_p != NULL)
  remote_data_t* remote_data = *remote_data_p;

  // Finish up the data transfer.
  // FIXME
  
  // Interpolate / transform as necessary.
  patch_buffer_transform(remote_data->buffer);

  // Copy the data into place.
  patch_buffer_copy_out(remote_data->buffer, remote_data->grid_data);
  
  // Dispose of the remote data for this token.
  int_ptr_unordered_map_delete(grid->pending_data, token);
}

int amr_grid_start_filling_ghosts(amr_grid_t* grid, amr_grid_data_t* data)
{
  ASSERT(amr_grid_data_grid(data) == grid);

  // Begin the remote copies and stick the given data into our "pending" list.
  int token = start_remote_copies(grid, data);

  // Perform all necessary local copies.
  do_local_copies(grid, data);

  return token;

#if 0
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
      finisher->i_dest = starter->i_dest;
      finisher->j_dest = starter->j_dest;
      finisher->k_dest = starter->k_dest;
      finisher->i_src = starter->i_src;
      finisher->j_src = starter->j_src;
      finisher->k_src = starter->k_src;
      ptr_array_append_with_dtor(grid->ghost_filler_finishers, finisher, polymec_free);
    }
  }
#endif
}

void amr_grid_finish_filling_ghosts(amr_grid_t* grid, int token)
{
  finish_remote_copies(grid, token);
#if 0
  // Go through our existing list of "finishers."
  for (int i = 0; i < grid->ghost_filler_finishers->size; ++i)
  {
    amr_grid_ghost_filler_finisher_t* finisher = grid->ghost_filler_finishers->data[i];
    amr_patch_t* dest_patch = amr_grid_data_patch(data, finisher->i_dest, finisher->j_dest, finisher->k_dest);
    ASSERT(dest_patch != NULL);
    amr_patch_t* src_patch = amr_grid_data_patch(data, finisher->i_src, finisher->j_src, finisher->k_src);
    finisher->method(grid, finisher->token, 
                     dest_patch, 
                     finisher->i_dest, finisher->j_dest, finisher->k_dest,
                     src_patch, 
                     finisher->i_src, finisher->j_src, finisher->k_src);
  }

  // Clear the list.
  ptr_array_clear(grid->ghost_filler_finishers);
#endif
}

void amr_grid_fill_ghosts(amr_grid_t* grid, amr_grid_data_t* data)
{
  int token = amr_grid_start_filling_ghosts(grid, data);
  amr_grid_finish_filling_ghosts(grid, token);
}

adj_graph_t* graph_from_amr_grid_patches(amr_grid_t* grid)
{
  // Create a graph whose vertices are the grid's local patches.
  adj_graph_t* g = adj_graph_new(grid->comm, grid->num_local_patches);

  // Allocate space in the graph for the edges (connections between patches).
  for (int i = 0; i < grid->num_local_patches; ++i)
  {
//    adj_graph_set_num_edges(g, i, num_connected_patches);
  }
  
  // Now fill in the edges.
  for (int i = 0; i < grid->num_local_patches; ++i)
  {
  }

  return g;
}

static size_t amr_grid_byte_size(void* obj)
{
  amr_grid_t* grid = obj;
  
  size_t storage =
    // Intrinsic metadata.
    7*sizeof(int) + 3*sizeof(bool) + 

    // Information about which patches are present.
    (1 + 2 * (grid->nx * grid->ny + grid->nz)) * sizeof(int);
    
  // Everything else (grid associations, ghost fillers) cannot be 
  // serialized, so the deserialized grid is manifestly not finalized.
  return storage;
}

static void* amr_grid_byte_read(byte_array_t* bytes, size_t* offset)
{
  // Read the intrinsic metadata.
  int nx, ny, nz, px, py, pz, ng;
  bool x_periodic, y_periodic, z_periodic;
  byte_array_read_ints(bytes, 1, &nx, offset);
  byte_array_read_ints(bytes, 1, &ny, offset);
  byte_array_read_ints(bytes, 1, &nz, offset);
  byte_array_read_ints(bytes, 1, &px, offset);
  byte_array_read_ints(bytes, 1, &py, offset);
  byte_array_read_ints(bytes, 1, &pz, offset);
  byte_array_read_ints(bytes, 1, &ng, offset);
  byte_array_read_ints(bytes, 1, (int*)&x_periodic, offset);
  byte_array_read_ints(bytes, 1, (int*)&y_periodic, offset);
  byte_array_read_ints(bytes, 1, (int*)&z_periodic, offset);

  // Create the empty grid.
  amr_grid_t* grid = amr_grid_new(MPI_COMM_WORLD, nx, ny, nz, px, py, pz,
                                  ng, x_periodic, y_periodic, z_periodic);
                                      
  // Read all the patch metadata.
  byte_array_read_ints(bytes, 1, &grid->num_local_patches, offset);
  byte_array_read_ints(bytes, nx*ny*nz, (int*)grid->patch_types, offset);
  byte_array_read_ints(bytes, nx*ny*nz, grid->remote_owners, offset);

  return grid;
}

static void amr_grid_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  amr_grid_t* grid = obj;

  // Write the intrinsic metadata.
  byte_array_write_ints(bytes, 1, &grid->nx, offset);
  byte_array_write_ints(bytes, 1, &grid->ny, offset);
  byte_array_write_ints(bytes, 1, &grid->nz, offset);
  byte_array_write_ints(bytes, 1, &grid->px, offset);
  byte_array_write_ints(bytes, 1, &grid->py, offset);
  byte_array_write_ints(bytes, 1, &grid->pz, offset);
  byte_array_write_ints(bytes, 1, &grid->ng, offset);
  byte_array_write_ints(bytes, 1, (int*)&grid->x_periodic, offset);
  byte_array_write_ints(bytes, 1, (int*)&grid->y_periodic, offset);
  byte_array_write_ints(bytes, 1, (int*)&grid->z_periodic, offset);

  // Write all the patch metadata.
  byte_array_write_ints(bytes, 1, &grid->num_local_patches, offset);
  byte_array_write_ints(bytes, grid->nx*grid->ny*grid->nz, (int*)grid->patch_types, offset);
  byte_array_write_ints(bytes, grid->nx*grid->ny*grid->nz, grid->remote_owners, offset);
}

serializer_t* amr_grid_serializer()
{
  return serializer_new("amr_grid", amr_grid_byte_size, amr_grid_byte_read, amr_grid_byte_write, NULL);
}

