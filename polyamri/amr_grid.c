// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/unordered_map.h"
#include "core/exchanger.h"
#include "polyamri/amr_grid.h"
#include "polyamri/amr_grid_data.h"
#include "polyamri/amr_grid_interpolator.h"

// This is a descriptor for "types" of patches -- this determines the nature 
// of interacting patches during ghost cell filling.
typedef enum
{
  NO_PATCH = 1,
  LOCAL_SAME_LEVEL = 2,
  LOCAL_FINER_LEVEL = 4,
  LOCAL_COARSER_LEVEL = 8,
  REMOTE_SAME_LEVEL = 16,
  REMOTE_FINER_LEVEL = 32,
  REMOTE_COARSER_LEVEL = 64
} patch_type_t;

struct amr_grid_t
{
  // Intrinsic metadata.
  int nx, ny, nz, px, py, pz;
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

  // Grid exchangers.
  exchanger_t* cell_ex;
  exchanger_t* x_face_ex;
  exchanger_t* y_face_ex;
  exchanger_t* z_face_ex;
  exchanger_t* x_edge_ex;
  exchanger_t* y_edge_ex;
  exchanger_t* z_edge_ex;
  exchanger_t* node_ex;

  // Local patch buffer for assisting in local copies, one list for each
  // type of centering.
  ptr_array_t* local_buffers[8];

  // Pending remote data from other processes.
  int_ptr_unordered_map_t* pending_data;
};

//------------------------------------------------------------------------
//                              Patch buffer
//------------------------------------------------------------------------

// This is a buffer for storing serialized patch data.
typedef struct
{
  int num_patches;
  int* patch_send_offsets;
  int* patch_receive_offsets;
  real_t* data;
} patch_buffer_t;

// Returns the number of values that patch (i, j, k) transmits to patch (i1, j1, k1).
// The patch mask identifies those kinds of patches included in the total.
static int num_transmitted_patch_values(amr_grid_t* grid, 
                                        amr_grid_data_centering_t centering, int patch_mask, 
                                        int i, int j, int k, int i1, int j1, int k1)
{
  // Make sure the patches abut one another in index space.
  if ((ABS(i - i1) + ABS(j - j1) + ABS(k - k1)) > 1)
    return 0;

  // Figure out how many values are on a side.
  int x_padding, y_padding, z_padding;
  switch(centering)
  {
    case AMR_GRID_CELL: 
      x_padding = y_padding = z_padding = 0; break;
    case AMR_GRID_X_FACE: 
      x_padding = 1; y_padding = z_padding = 0; break;
    case AMR_GRID_Y_FACE: 
      y_padding = 1; x_padding = z_padding = 0; break;
    case AMR_GRID_Z_FACE: 
      z_padding = 1; x_padding = y_padding = 0; break;
    case AMR_GRID_X_EDGE: 
      x_padding = 0; y_padding = z_padding = 1; break;
    case AMR_GRID_Y_EDGE: 
      y_padding = 0; x_padding = z_padding = 1; break;
    case AMR_GRID_Z_EDGE: 
      z_padding = 0; x_padding = y_padding = 1; break;
    case AMR_GRID_NODE: 
      x_padding = y_padding = z_padding = 1;
  }
  int px, py, pz;
  amr_grid_get_patch_size(grid, &px, &py, &pz);
  int x_size = px + x_padding, y_size = py + y_padding, z_size = pz + z_padding;

  // Handle neighbor grid cases, if any.
  amr_grid_t* grid1 = grid;
  int cross_section_size = 0;
  if (i1 == -1)
  {
    grid1 = (grid->neighbors[0] != NULL) ? grid->neighbors[0] : NULL;
    i1 = grid1->nx - 1;
    cross_section_size = y_size * z_size;
  }
  else if (i1 == grid->px)
  {
    grid1 = (grid->neighbors[1] != NULL) ? grid->neighbors[1] : NULL;
    i1 = 0;
    cross_section_size = y_size * z_size;
  }
  else if (j1 == -1)
  {
    grid1 = (grid->neighbors[2] != NULL) ? grid->neighbors[2] : NULL;
    j1 = grid1->ny - 1;
    cross_section_size = x_size * z_size;
  }
  else if (j1 == grid->py)
  {
    grid1 = (grid->neighbors[3] != NULL) ? grid->neighbors[3] : NULL;
    j1 = 0;
    cross_section_size = x_size * z_size;
  }
  else if (k1 == -1)
  {
    grid1 = (grid->neighbors[4] != NULL) ? grid->neighbors[4] : NULL;
    k1 = grid1->nz - 1;
    cross_section_size = x_size * y_size;
  }
  else if (k1 == grid->px)
  {
    grid1 = (grid->neighbors[5] != NULL) ? grid->neighbors[5] : NULL;
    k1 = 0;
    cross_section_size = x_size * y_size;
  }

  if (grid1 == NULL)
    return 0;

  DECLARE_3D_ARRAY(patch_type_t, patch_types, grid1->patch_types, grid1->nx, grid1->ny, grid1->nz);
  patch_type_t patch_type = patch_types[i1][j1][k1];
  if ((patch_mask | patch_type) == 0)
    return 0;

  return cross_section_size;
}

// General patch buffer constructor.
static patch_buffer_t* patch_buffer_new(amr_grid_t* grid, 
                                        amr_grid_data_centering_t centering,
                                        int num_components, int num_ghosts,
                                        int patch_mask)
{
  ASSERT(num_components > 0);
  ASSERT(num_ghosts >= 0);

  patch_buffer_t* buffer = polymec_malloc(sizeof(patch_buffer_t));
  buffer->num_patches = grid->num_local_patches;
  buffer->patch_send_offsets = polymec_malloc(sizeof(int) * (buffer->num_patches + 1));
  buffer->patch_receive_offsets = polymec_malloc(sizeof(int) * (buffer->num_patches + 1));

  // The ordering of the patch buffer's data is as follows:
  //
  // 1. Sent data: all local patches abutting other local patches (on this grid, on a 
  //    neighboring grid, or on a coarser or finer grid) will send the appropriate number of 
  //    layers of interior values to these other local patches. The patches are traversed in the
  //    order that they are locally stored, and each patch will pack data from its -x, +x, -y, +y, 
  //    and -z, +z boundaries (in that order).
  // 2. Received data: these local patches then expect the same number of layers of ghost values
  //    from the same local patches, traversed in the same order.

  // Now count up the sent data for each patch.
  int local_mask = LOCAL_SAME_LEVEL | LOCAL_FINER_LEVEL | LOCAL_COARSER_LEVEL;
  int pos = 0, i, j, k, p = 0, l = 0;
  while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k))
  {
    // Mark our spot in the buffer.
    buffer->patch_send_offsets[p] = l;

    // Figure out what we need to send to our neighbors.
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i-1, j, k);
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i+1, j, k);
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j-1, k);
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j+1, k);
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j, k-1);
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j, k+1);

    ++p;
  }
  buffer->patch_send_offsets[p] = l;

  // Count up the receive data for each patch.
  pos = 0;
  p = 0;
  while (amr_grid_next_local_patch(grid, &pos, &i, &j, &k))
  {
    // Mark our spot in the buffer.
    buffer->patch_receive_offsets[p] = l;

    // Figure out what we need to send to our neighbors.
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i-1, j, k);
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i+1, j, k);
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j-1, k);
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j+1, k);
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j, k-1);
    l += num_ghosts * num_components * num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j, k+1);

    ++p;
  }
  buffer->patch_receive_offsets[p] = l;

  // Allocate storage for the send/receive buffer.
  buffer->data = polymec_malloc(sizeof(real_t) * l);

  return buffer;
}

// Creates a patch buffer that will perform copies between local patches
// in the given grid.
static patch_buffer_t* local_patch_buffer_new(amr_grid_t* grid, 
                                              amr_grid_data_centering_t centering,
                                              int num_components, int num_ghosts)
{
  int local_mask = LOCAL_SAME_LEVEL | LOCAL_FINER_LEVEL | LOCAL_COARSER_LEVEL;
  return patch_buffer_new(grid, centering, num_components, num_ghosts, local_mask);
}

// Creates a patch buffer that will perform copies between local and remote 
// patches in the grid underlying the given grid data object.
static patch_buffer_t* remote_patch_buffer_new(amr_grid_data_t* grid_data)
{
  amr_grid_t* grid = amr_grid_data_grid(grid_data);
  amr_grid_data_centering_t centering = amr_grid_data_centering(grid_data);
  int num_components = amr_grid_data_num_components(grid_data);
  int num_ghosts = amr_grid_data_num_ghosts(grid_data);
  int remote_mask = REMOTE_SAME_LEVEL | REMOTE_FINER_LEVEL | REMOTE_COARSER_LEVEL;
  return patch_buffer_new(grid, centering, num_components, num_ghosts, remote_mask);
}

// Destroys the given patch buffer.
static void patch_buffer_free(patch_buffer_t* buffer)
{
  polymec_free(buffer->patch_send_offsets);
  polymec_free(buffer->patch_receive_offsets);
  polymec_free(buffer->data);
  polymec_free(buffer);
}

// Copies data from the patches in the given grid to the patch buffer.
static void patch_buffer_copy_in(patch_buffer_t* buffer, amr_grid_data_t* grid_data)
{
  amr_grid_t* grid = amr_grid_data_grid(grid_data);
  amr_grid_data_centering_t centering = amr_grid_data_centering(grid_data);
  int num_components = amr_grid_data_num_components(grid_data);
  int num_ghosts = amr_grid_data_num_ghosts(grid_data);

  // Now count up the sent data for each patch.
  int local_mask = LOCAL_SAME_LEVEL | LOCAL_FINER_LEVEL | LOCAL_COARSER_LEVEL;
  int pos = 0, i, j, k, p = 0;
  amr_patch_t* patch;
  DECLARE_3D_ARRAY(patch_type_t, patch_types, grid->patch_types, grid->nx, grid->ny, grid->nz);
  while (amr_grid_data_next_local_patch(grid_data, &pos, &i, &j, &k, &patch))
  {
    patch_type_t patch_type = patch_types[i][j][k];
    if (patch_type == LOCAL_SAME_LEVEL)
      polymec_error("patch_buffer_copy_in: adaptive mesh refinement is not yet supported!");

    DECLARE_AMR_PATCH_ARRAY(array, patch);
    int offset = buffer->patch_send_offsets[p];

    // Copy in -x values.
    int num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i-1, j, k);
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            buffer->data[offset] = array[patch->i1+g][j][k][c];

    // Copy in +x values.
    num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i+1, j, k);
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            buffer->data[offset] = array[patch->i2-g-1][j][k][c];

    // Copy in -y values.
    num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j-1, k);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            buffer->data[offset] = array[i][patch->j1+g][k][c];

    // Copy in +y values.
    num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j+1, k);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            buffer->data[offset] = array[i][patch->j2-g-1][k][c];

    // Copy in -z values.
    num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j, k-1);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int j = patch->j1; j < patch->j2; ++j)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            buffer->data[offset] = array[i][j][patch->k1+g][c];

    // Copy in +z values.
    num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j, k+1);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int j = patch->j1; j < patch->j2; ++j)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            buffer->data[offset] = array[i][j][patch->k2-g-1][c];

    ASSERT(offset == buffer->patch_send_offsets[p+1]);
    ++p;
  }
}

// Copies data from the patch buffer to the patches in the given grid.
static void patch_buffer_copy_out(patch_buffer_t* buffer, amr_grid_data_t* grid_data)
{
  amr_grid_t* grid = amr_grid_data_grid(grid_data);
  amr_grid_data_centering_t centering = amr_grid_data_centering(grid_data);
  int num_components = amr_grid_data_num_components(grid_data);
  int num_ghosts = amr_grid_data_num_ghosts(grid_data);

  // Now count up the sent data for each patch.
  int local_mask = LOCAL_SAME_LEVEL | LOCAL_FINER_LEVEL | LOCAL_COARSER_LEVEL;
  int pos = 0, i, j, k, p = 0;
  amr_patch_t* patch;
  DECLARE_3D_ARRAY(patch_type_t, patch_types, grid->patch_types, grid->nx, grid->ny, grid->nz);
  while (amr_grid_data_next_local_patch(grid_data, &pos, &i, &j, &k, &patch))
  {
    patch_type_t patch_type = patch_types[i][j][k];
    if (patch_type == LOCAL_SAME_LEVEL)
      polymec_error("patch_buffer_copy_out: adaptive mesh refinement is not yet supported!");

    DECLARE_AMR_PATCH_ARRAY(array, patch);
    int offset = buffer->patch_receive_offsets[p];

    // Copy out -x values.
    int num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i-1, j, k);
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            array[g][j][k][c] = buffer->data[offset];

    // Copy out +x values.
    num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i+1, j, k);
    for (int j = patch->j1; j < patch->j2; ++j)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            array[patch->i2+num_ghosts-g-1][j][k][c] = buffer->data[offset];

    // Copy out -y values.
    num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j-1, k);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            array[i][g][k][c] = buffer->data[offset];

    // Copy out +y values.
    num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j+1, k);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int k = patch->k1; k < patch->k2; ++k)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            array[i][patch->j2+num_ghosts-g-1][k][c] = buffer->data[offset];

    // Copy out -z values.
    num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j, k-1);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int j = patch->j1; j < patch->j2; ++j)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            array[i][j][g][c] = buffer->data[offset];

    // Copy out +z values.
    num_values = num_transmitted_patch_values(grid, centering, local_mask, i, j, k, i, j, k+1);
    for (int i = patch->i1; i < patch->i2; ++i)
      for (int j = patch->j1; j < patch->j2; ++j)
        for (int g = 0; g < num_ghosts; ++g)
          for (int c = 0; c < num_components; ++c, ++offset)
            array[i][j][patch->k2+num_ghosts-g-1][c] = buffer->data[offset];

    ASSERT(offset == buffer->patch_receive_offsets[p+1]);
    ++p;
  }
}

//------------------------------------------------------------------------
//                              Remote data
//------------------------------------------------------------------------

// This type stores a buffer and a pointer to grid data for pending remote 
// operations. It's really just a container with a constructor and destructor.
typedef struct
{
  patch_buffer_t* buffer;
  amr_grid_data_t* grid_data;
  exchanger_t* ex;
} remote_data_t;

// Creates a remote data object for the given grid data.
static remote_data_t* remote_data_new(amr_grid_data_t* grid_data)
{
  remote_data_t* data = polymec_malloc(sizeof(remote_data_t));
  data->buffer = remote_patch_buffer_new(grid_data);
  data->grid_data = grid_data;

  // Figure out the right exchanger for this data.
  amr_grid_t* grid = amr_grid_data_grid(grid_data);
  amr_grid_data_centering_t centering = amr_grid_data_centering(grid_data);
  switch(centering)
  {
    case AMR_GRID_CELL: data->ex = grid->cell_ex; break;
    case AMR_GRID_X_FACE: data->ex = grid->x_face_ex; break;
    case AMR_GRID_Y_FACE: data->ex = grid->y_face_ex; break;
    case AMR_GRID_Z_FACE: data->ex = grid->z_face_ex; break;
    case AMR_GRID_X_EDGE: data->ex = grid->x_edge_ex; break;
    case AMR_GRID_Y_EDGE: data->ex = grid->y_edge_ex; break;
    case AMR_GRID_Z_EDGE: data->ex = grid->z_edge_ex; break;
    case AMR_GRID_NODE: data->ex = grid->node_ex;
  }

  return data;
}

static void remote_data_free(remote_data_t* remote_data)
{
  patch_buffer_free(remote_data->buffer);
  polymec_free(remote_data);
}

//------------------------------------------------------------------------
//                          Exchanger constructors
//------------------------------------------------------------------------

static exchanger_t* grid_cell_exchanger_new(amr_grid_t* grid)
{
  return NULL;
}

static exchanger_t* grid_x_face_exchanger_new(amr_grid_t* grid)
{
  return NULL;
}

static exchanger_t* grid_y_face_exchanger_new(amr_grid_t* grid)
{
  return NULL;
}

static exchanger_t* grid_z_face_exchanger_new(amr_grid_t* grid)
{
  return NULL;
}

static exchanger_t* grid_x_edge_exchanger_new(amr_grid_t* grid)
{
  return NULL;
}

static exchanger_t* grid_y_edge_exchanger_new(amr_grid_t* grid)
{
  return NULL;
}

static exchanger_t* grid_z_edge_exchanger_new(amr_grid_t* grid)
{
  return NULL;
}

static exchanger_t* grid_node_exchanger_new(amr_grid_t* grid)
{
  return NULL;
}

//------------------------------------------------------------------------
//                              AMR grid
//------------------------------------------------------------------------

amr_grid_t* amr_grid_new(MPI_Comm comm,
                         int nx, int ny, int nz, 
                         int px, int py, int pz, 
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

  amr_grid_t* grid = polymec_malloc(sizeof(amr_grid_t));
  grid->nx = nx;
  grid->ny = ny;
  grid->nz = nz;
  grid->px = px;
  grid->py = py;
  grid->pz = pz;
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
  grid->cell_ex = NULL;
  grid->x_face_ex = NULL;
  grid->y_face_ex = NULL;
  grid->z_face_ex = NULL;
  grid->x_edge_ex = NULL;
  grid->y_edge_ex = NULL;
  grid->z_edge_ex = NULL;
  grid->node_ex = NULL;

  for (int centering = 0; centering < 8; ++centering)
    grid->local_buffers[centering] = ptr_array_new();
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

  if (grid->cell_ex != NULL)
    exchanger_free(grid->cell_ex);
  if (grid->x_face_ex != NULL)
    exchanger_free(grid->x_face_ex);
  if (grid->y_face_ex != NULL)
    exchanger_free(grid->y_face_ex);
  if (grid->z_face_ex != NULL)
    exchanger_free(grid->z_face_ex);
  if (grid->x_edge_ex != NULL)
    exchanger_free(grid->x_edge_ex);
  if (grid->y_edge_ex != NULL)
    exchanger_free(grid->y_edge_ex);
  if (grid->z_edge_ex != NULL)
    exchanger_free(grid->z_edge_ex);
  if (grid->node_ex != NULL)
    exchanger_free(grid->node_ex);

  for (int centering = 0; centering < 8; ++centering)
    ptr_array_free(grid->local_buffers[centering]);

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

  // Establish our remote copying pattern.
  grid->cell_ex = grid_cell_exchanger_new(grid);
  grid->x_face_ex = grid_x_face_exchanger_new(grid);
  grid->y_face_ex = grid_y_face_exchanger_new(grid);
  grid->z_face_ex = grid_z_face_exchanger_new(grid);
  grid->x_edge_ex = grid_x_edge_exchanger_new(grid);
  grid->y_edge_ex = grid_y_edge_exchanger_new(grid);
  grid->z_edge_ex = grid_z_edge_exchanger_new(grid);
  grid->node_ex = grid_node_exchanger_new(grid);
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

void amr_grid_get_patch_size(amr_grid_t* grid, int* pnx, int* pny, int* pnz)
{
  *pnx = grid->px;
  *pny = grid->py;
  *pnz = grid->pz;
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

static void do_local_copies(amr_grid_t* grid, amr_grid_data_t* data)
{
  // Make sure we have a local patch buffer to use with the right number of components.
  amr_grid_data_centering_t centering = amr_grid_data_centering(data);
  int num_comp = amr_grid_data_num_components(data);
  int num_ghosts = amr_grid_data_num_ghosts(data);
  ptr_array_t* local_buffers = grid->local_buffers[centering];
  if (num_comp > local_buffers->size)
  {
    int old_size = local_buffers->size;
    ptr_array_resize(local_buffers, num_comp);
    for (int i = old_size; i < num_comp; ++i)
      local_buffers->data[i] = NULL;
  }
  if (local_buffers->data[num_comp-1] == NULL)
  {
    patch_buffer_t* new_buffer = local_patch_buffer_new(grid, centering, num_comp, num_ghosts);
    ptr_array_assign_with_dtor(local_buffers, num_comp-1, new_buffer, DTOR(patch_buffer_free));
  }
  patch_buffer_t* buffer = local_buffers->data[num_comp-1];

  // Copy the data into our local patch buffer and then out again.
  patch_buffer_copy_in(buffer, data);

  // Copy out.
  patch_buffer_copy_out(buffer, data);
}

static int start_remote_copies(amr_grid_t* grid, amr_grid_data_t* data)
{
  // Make a new remote data thingy.
  remote_data_t* remote_data = remote_data_new(data);

  // Copy the grid data into the thingy's patch buffer.
  patch_buffer_t* buffer = remote_data->buffer;
  patch_buffer_copy_in(buffer, data);

  // Initiate the data transfer and receive a token for the exchange.
  int num_comp = amr_grid_data_num_components(data);
  int token = exchanger_start_exchange(remote_data->ex, buffer->data, num_comp, 0, MPI_REAL_T);

  // Stash the remote data into our set of pending data.
  int_ptr_unordered_map_insert_with_v_dtor(grid->pending_data, token, remote_data, DTOR(remote_data_free));

  return token;
}

int amr_grid_start_filling_ghosts(amr_grid_t* grid, amr_grid_data_t* data)
{
  ASSERT(amr_grid_data_grid(data) == grid);

  // Begin the remote copies and stick the given data into our "pending" list.
  int token = start_remote_copies(grid, data);

  // Perform all necessary local copies.
  do_local_copies(grid, data);

  return token;
}

static void finish_remote_copies(amr_grid_t* grid, int token)
{
  remote_data_t** remote_data_p = (remote_data_t**)int_ptr_unordered_map_get(grid->pending_data, token);
  ASSERT(remote_data_p != NULL)
  remote_data_t* remote_data = *remote_data_p;

  // Finish up the data transfer.
  exchanger_finish_exchange(remote_data->ex, token);
  
  // Copy the data into place.
  patch_buffer_copy_out(remote_data->buffer, remote_data->grid_data);
  
  // Dispose of the remote data for this token.
  int_ptr_unordered_map_delete(grid->pending_data, token);
}

void amr_grid_finish_filling_ghosts(amr_grid_t* grid, int token)
{
  finish_remote_copies(grid, token);
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
  int nx, ny, nz, px, py, pz;
  bool x_periodic, y_periodic, z_periodic;
  byte_array_read_ints(bytes, 1, &nx, offset);
  byte_array_read_ints(bytes, 1, &ny, offset);
  byte_array_read_ints(bytes, 1, &nz, offset);
  byte_array_read_ints(bytes, 1, &px, offset);
  byte_array_read_ints(bytes, 1, &py, offset);
  byte_array_read_ints(bytes, 1, &pz, offset);
  byte_array_read_ints(bytes, 1, (int*)&x_periodic, offset);
  byte_array_read_ints(bytes, 1, (int*)&y_periodic, offset);
  byte_array_read_ints(bytes, 1, (int*)&z_periodic, offset);

  // Create the empty grid.
  amr_grid_t* grid = amr_grid_new(MPI_COMM_WORLD, nx, ny, nz, px, py, pz,
                                  x_periodic, y_periodic, z_periodic);
                                      
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

