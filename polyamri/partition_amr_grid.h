// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_PARTITION_AMR_GRID_H
#define POLYAMRI_PARTITION_AMR_GRID_H

#include "core/exchanger.h"
#include "geometry/coord_mapping.h"
#include "polyamri/amr_grid.h"

// This function partitions the given AMR grid on rank 0 with the given 
// load weights, distributing the patches to parallel domains on the given 
// communicator to balance the load. If a non-NULL coordinate mapping is given,
// it is used to provide spatial information about how the patches If weights is NULL, the patches are assumed 
// all to have equal weights. The weights array should be interpreted as a 3-dimensional
// array that assigns an integer to each patch in the given grid. This function creates 
// and returns an exchanger object that can be used to distribute data from the rank 0 
// to the partition. The grid on rank 0 (as well as any non-NULL grid on rank != 0) is 
// consumed. In each case, the grid is replaced with a partitioned grid.
exchanger_t* partition_amr_grid(amr_grid_t** grid, MPI_Comm comm, int* weights, real_t imbalance_tol);

// This function repartitions the given mesh with the given load weights, 
// alloting the patches to parallel domains to balance their load. If weights is 
// NULL, the patches are all assumed to have equal weights. The weights array should 
// be interpreted as a 3-dimensional array that assigns an integer to each patch in the 
// given grid. This function creates and returns an exchanger object that can be used to 
// migrate data from the old partition to the new. The grid is consumed and replaced with a 
// repartitioned grid.
//exchanger_t* repartition_amr_grid(amr_grid_t** grid, int* weights, real_t imbalance_tol);

// While partition_amr_grid and repartition_amr_grid are all-in-one grid partitioners, the 
// following functions allow one to mix-n-match the pieces of the underlying algorithms.

// This function creates a newly-allocated global partition vector that can be used 
// to distribute a global AMR grid on rank 0 to all processes on the given communicator, 
// according to the given (patch) weights and the specified imbalance tolerance. The partition 
// vector assigns a process ID to each patch in the grid. Grid objects on non-zero ranks are 
// ignored.
int64_t* partition_vector_from_amr_grid(amr_grid_t* global_grid, 
                                        MPI_Comm comm, 
                                        int* weights, 
                                        real_t imbalance_tol);

// Given a global partition vector, distributes the patches in the AMR grid from rank 0 to 
// all processes in the given communicator, and returns a distributor object that can be used 
// to distribute its data. The grid on rank 0 is replaced with a partitioned grid, and grids are 
// written (or overwritten) on other ranks.
exchanger_t* distribute_amr_grid(amr_grid_t** grid, MPI_Comm comm, int64_t* global_partition);

#endif

