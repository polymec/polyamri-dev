// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "polyamri/partition_amr_grid.h"

#if POLYMEC_HAVE_MPI

#include "core/adj_graph.h"

// We use an unpublished function in polymec to do our graph cutting.
extern int64_t* partition_graph(adj_graph_t* global_graph, 
                                MPI_Comm comm,
                                int* weights,
                                real_t imbalance_tol);

// Creates a grid from a subset of the patches of a global grid.
static amr_grid_t* create_subgrid(MPI_Comm comm,
                                  amr_grid_t* global_grid,
                                  int64_t* global_partition,
                                  index_t* vtx_dist,
                                  int* indices, int num_indices)
{
  // FIXME
  return NULL;
}

static void amr_grid_distribute(amr_grid_t** grid, 
                                MPI_Comm comm,
                                adj_graph_t* global_graph, 
                                int64_t* global_partition)
{
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Make sure we're all here.
  MPI_Barrier(comm);

  // We use the int_array serializer, so we need to make sure it's 
  // registered on all processes. See serializer.h for details.
  {
    serializer_t* s = int_array_serializer();
    s = NULL;
  }

  amr_grid_t* global_grid = *grid;
  amr_grid_t* local_grid = NULL;
  uint64_t vtx_dist[nprocs+1];
  if (rank == 0)
  {
    // Take stock of how many patches we'll have per process.
    int num_global_patches = amr_grid_num_local_patches(global_grid);
    int num_local_patches[nprocs];
    memset(num_local_patches, 0, sizeof(int) * nprocs);
    for (int i = 0; i < num_global_patches; ++i)
      num_local_patches[global_partition[i]]++;

    // Construct the distribution of vertices for the partitioning.
    vtx_dist[0] = 0;
    for (int p = 0; p < nprocs; ++p)
      vtx_dist[p+1] = vtx_dist[p] + num_local_patches[p];

    // Carve out the patches that will stick around on process 0.
    {
      int indices[num_local_patches[0]], k = 0;
      for (int i = 0; i < num_global_patches; ++i)
      {
        if (global_partition[i] == rank)
          indices[k++] = i;
      }
      local_grid = create_subgrid(comm, global_grid, global_partition, NULL, indices, num_local_patches[0]);
    }

    // Now do the other processes.
    serializer_t* ser = amr_grid_serializer();
    byte_array_t* bytes = byte_array_new();
    for (int p = 1; p < nprocs; ++p)
    {
      // Share the vertex distribution.
      MPI_Send(vtx_dist, nprocs+1, MPI_UINT64_T, p, p, comm);

      // Create the pth subgrid.
      int indices[num_local_patches[p]], k = 0;
      for (int i = 0; i < num_global_patches; ++i)
      {
        if (global_partition[i] == p)
          indices[k++] = i;
      }
      amr_grid_t* p_grid = create_subgrid(comm, global_grid, global_partition, NULL, indices, num_local_patches[p]);

      // Serialize it and send its size (and it) to process p.
      size_t offset = 0;
      serializer_write(ser, p_grid, bytes, &offset);
      MPI_Send(&bytes->size, 1, MPI_INT, p, p, comm);
      MPI_Send(bytes->data, bytes->size, MPI_BYTE, p, p, comm);

      // Clean up.
      byte_array_clear(bytes);
      amr_grid_free(p_grid);
    }
    ser = NULL;
    byte_array_free(bytes);
  }
  else
  {
    // Receive the vertex distribution of the incoming mesh.
    MPI_Status status;
    MPI_Recv(vtx_dist, nprocs+1, MPI_UINT64_T, 0, rank, comm, &status);

    // Receive the size of the incoming mesh.
    int mesh_size;
    MPI_Recv(&mesh_size, 1, MPI_INT, 0, rank, comm, &status);

    // Now receive the mesh.
    byte_array_t* bytes = byte_array_new();
    byte_array_resize(bytes, mesh_size);

    MPI_Recv(bytes->data, mesh_size, MPI_BYTE, 0, rank, comm, &status);
    serializer_t* ser = amr_grid_serializer();
    size_t offset = 0;
    local_grid = serializer_read(ser, bytes, &offset);
    
    byte_array_free(bytes);
    ser = NULL;
  }

  *grid = local_grid;

  // Clean up.
  if (global_grid != NULL)
    amr_grid_free(global_grid);

#if 0
  // Extract the boundary faces and the original (global) ghost cell indices
  // associated with them.
  ASSERT(mesh_has_tag(local_mesh->face_tags, "parallel_boundary_faces"));
  int num_pbfaces;
  int* pbfaces = mesh_tag(local_mesh->face_tags, "parallel_boundary_faces", &num_pbfaces);
  int_array_t* pbgcells = mesh_tag_property(local_mesh->face_tags, "parallel_boundary_faces", 
                                            "ghost_cell_indices");
  ASSERT(pbgcells != NULL);
  int_array_t* global_cell_indices = mesh_property(local_mesh, "global_cell_indices");
  ASSERT(global_cell_indices != NULL);

  // Here's an inverse map for global cell indices to the local ones we have.
  int_int_unordered_map_t* inverse_cell_map = int_int_unordered_map_new();
  for (int i = 0; i < local_mesh->num_cells; ++i)
    int_int_unordered_map_insert(inverse_cell_map, global_cell_indices->data[i], i);

  // Now we create pairwise global cell indices for the exchanger.
  int_ptr_unordered_map_t* ghost_cell_indices = int_ptr_unordered_map_new();
  int next_ghost_index = local_mesh->num_cells;
  for (int i = 0; i < num_pbfaces; ++i)
  {
    int f = pbfaces[i];
    ASSERT(local_mesh->face_cells[2*f+1] < -1);
    // Get the destination process for this ghost cell.
    int proc = -local_mesh->face_cells[2*f+1] - 2;
    ASSERT(proc >= 0);
    ASSERT(proc < nprocs);

    // Generate a mapping from a global index to its ghost cell.
    int global_ghost_index = pbgcells->data[i];
    int local_ghost_index;
    int* ghost_index_p = int_int_unordered_map_get(inverse_cell_map, global_ghost_index);
    if (ghost_index_p == NULL)
    {
      local_ghost_index = next_ghost_index++;
      int_int_unordered_map_insert(inverse_cell_map, global_ghost_index, local_ghost_index);
    }
    else
      local_ghost_index = *ghost_index_p;
    local_mesh->face_cells[2*f+1] = local_ghost_index;

    if (!int_ptr_unordered_map_contains(ghost_cell_indices, proc))
      int_ptr_unordered_map_insert_with_v_dtor(ghost_cell_indices, proc, int_array_new(), DTOR(int_array_free));
    int_array_t* indices = *int_ptr_unordered_map_get(ghost_cell_indices, proc);
    int local_cell = local_mesh->face_cells[2*f];
    int global_cell = global_cell_indices->data[local_cell];
    int_array_append(indices, global_cell);
    int_array_append(indices, global_ghost_index);
  }

  // Make sure everything lines up.
  ASSERT(next_ghost_index - local_mesh->num_cells == local_mesh->num_ghost_cells);

  exchanger_t* ex = mesh_exchanger(local_mesh);
  int pos = 0, proc;
  int_array_t* indices;
  while (int_ptr_unordered_map_next(ghost_cell_indices, &pos, &proc, (void**)&indices))
  {
    if (proc != rank)
    {
      ASSERT((indices->size > 0) && ((indices->size % 2) == 0));
      int num_pairs = indices->size/2;

      // Sort the indices array lexicographically by pairs so that all of the 
      // exchanger send/receive transactions have the same order across 
      // processes. This requires a specialized sort, since we have to 
      // arrange the integers within the pairs in ascending order, sort 
      // them, and then switch them back.
      sort_global_cell_pairs(indices->data, num_pairs);
      
      int send_indices[num_pairs], recv_indices[num_pairs];
      for (int i = 0; i < num_pairs; ++i)
      {
        send_indices[i] = *int_int_unordered_map_get(inverse_cell_map, indices->data[2*i]);
        ASSERT(send_indices[i] < local_mesh->num_cells);
        recv_indices[i] = *int_int_unordered_map_get(inverse_cell_map, indices->data[2*i+1]);
        ASSERT(recv_indices[i] >= local_mesh->num_cells);
      }
      exchanger_set_send(ex, proc, send_indices, num_pairs, true);
      exchanger_set_receive(ex, proc, recv_indices, num_pairs, true);
    }
  }

  // Clean up again.
  int_ptr_unordered_map_free(ghost_cell_indices);
  int_int_unordered_map_free(inverse_cell_map);

  // Destroy the tag and global cell index properties. We don't want to have 
  // to rely on them further.
  mesh_delete_tag(local_mesh->face_tags, "parallel_boundary_faces");
  mesh_delete_property(local_mesh, "global_cell_indices");
#endif
  STOP_FUNCTION_TIMER();
}
#endif

exchanger_t* partition_amr_grid(amr_grid_t** grid, 
                                MPI_Comm comm, 
                                int* weights, 
                                real_t imbalance_tol)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  ASSERT((*grid == NULL) || (amr_grid_comm(*grid) == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (*grid != NULL));

  // On a single process, partitioning has no meaning, but we do replace the communicator
  // if needed. NOTE: the exchanger will still have its original communicator, but this 
  // shouldn't matter in any practical sense.
  if (nprocs == 1)
  {
    if (comm != amr_grid_comm(*grid))
      amr_grid_set_comm(*grid, comm);
    STOP_FUNCTION_TIMER();
    return exchanger_new(comm);
  }

  log_debug("partition_amr_grid: Partitioning grid into %d subdomains.", nprocs);

  // If grids on rank != 0 are not NULL, we delete them.
  amr_grid_t* g = *grid;
  if ((rank != 0) && (g != NULL))
  {
    amr_grid_free(g);
    *grid = g = NULL; 
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* global_graph = (g != NULL) ? graph_from_amr_grid_patches(g) : NULL;

#ifndef NDEBUG
  // Make sure there are enough patches to go around for the processes we're given.
  if (rank == 0)
  {
    ASSERT(amr_grid_num_local_patches(g) >= nprocs);
  }
#endif

  // Map the graph to the different domains, producing a local partition vector.
  int64_t* global_partition = (rank == 0) ? partition_graph(global_graph, comm, weights, imbalance_tol): NULL;

  // Distribute the grid.
  log_debug("partition_amr_grid: Distributing grid to %d processes.", nprocs);
  amr_grid_distribute(grid, comm, global_graph, global_partition);

  // Set up an exchanger to distribute field data.
  int num_vertices = (g != NULL) ? adj_graph_num_vertices(global_graph) : 0;
  exchanger_t* distributor = create_distributor(comm, global_partition, num_vertices);

  // Clean up.
  if (global_graph != NULL)
    adj_graph_free(global_graph);
  if (global_partition != NULL)
    polymec_free(global_partition);

  // Return the migrator.
  STOP_FUNCTION_TIMER();
  return distributor;
#else
  // Replace the communicator if needed.
  if (comm != amr_grid_comm(*grid))
    amr_grid_set_comm(*grid, comm);
  return exchanger_new(comm);
#endif
}

//exchanger_t* repartition_amr_grid(amr_grid_t** grid, int* weights, real_t imbalance_tol)
//{
//}

int64_t* partition_vector_from_amr_grid(amr_grid_t* global_grid, 
                                        MPI_Comm comm, 
                                        int* weights, 
                                        real_t imbalance_tol)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (global_grid != NULL));

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
  {
    // Dumb, but correct.
    int num_local_patches = amr_grid_num_local_patches(global_grid);
    int64_t* global_partition = polymec_malloc(sizeof(int64_t) * num_local_patches);
    memset(global_partition, 0, sizeof(int64_t) * num_local_patches);
    return global_partition;
  }

  // Generate a global adjacency graph for the patches in the grid.
  adj_graph_t* global_graph = (rank == 0) ? graph_from_amr_grid_patches(global_grid) : NULL;

#ifndef NDEBUG
  // Make sure there are enough patches to go around for the processes we're given.
  if (rank == 0)
  {
    ASSERT(amr_grid_num_local_patches(global_grid) >= nprocs);
  }
#endif

  // Map the graph to the different domains, producing a local partition vector.
  int64_t* global_partition = (rank == 0) ? partition_graph(global_graph, comm, weights, imbalance_tol): NULL;

  // Get rid of the graph.
  if (global_graph != NULL)
    adj_graph_free(global_graph);

  STOP_FUNCTION_TIMER();
  return global_partition;

#else
  // This is dumb, but we were asked for it.
  int num_local_patches = amr_grid_num_local_patches(global_grid);
  int64_t* global_partition = polymec_malloc(sizeof(int64_t) * num_local_patches);
  memset(global_partition, 0, sizeof(int64_t) * num_local_patches);
  return global_partition;
#endif
}

exchanger_t* distribute_amr_grid(amr_grid_t** grid, MPI_Comm comm, int64_t* global_partition)
{
#if POLYMEC_HAVE_MPI
  ASSERT((*grid == NULL) || (amr_grid_comm(*grid) == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (global_partition != NULL));
  ASSERT((rank != 0) || (*grid != NULL));

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(comm);

  START_FUNCTION_TIMER();

  // If grids on rank != 0 are not NULL, we delete them.
  amr_grid_t* g = *grid;
  if ((rank != 0) && (g != NULL))
  {
    amr_grid_free(g);
    *grid = g = NULL; 
  }

  // Generate a global adjacency graph for the grid.
  adj_graph_t* global_graph = (g != NULL) ? graph_from_amr_grid_patches(g) : NULL;

  // Distribute the grid.
  amr_grid_distribute(grid, comm, global_graph, global_partition);

  // Set up an exchanger to distribute field data.
  int num_vertices = (g != NULL) ? adj_graph_num_vertices(global_graph) : 0;
  exchanger_t* distributor = create_distributor(comm, global_partition, num_vertices);

  // Get rid of the graph.
  if (global_graph != NULL)
    adj_graph_free(global_graph);

  STOP_FUNCTION_TIMER();
  return distributor;
#else
  return exchanger_new(comm);
#endif
}

