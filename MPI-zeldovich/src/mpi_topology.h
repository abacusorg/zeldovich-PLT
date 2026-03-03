#ifndef MPI_TOPOLOGY_H
#define MPI_TOPOLOGY_H

#include <mpi.h>

/* Dimensions: [grid_x, grid_z] & grid_x * grid_z = num_ranks
 * Periodic: X and Z 
 * Reorder: Enabled (reorder=1) for potential hardware-aware optimization
 * Rank mapping: Linear rank = rank_x * grid_z + rank_z (row-major order)
 * 
 * REQUIREMENTS:
 * -------------
 * - num_ranks must be factorable (grid_x * grid_z == num_ranks)
 * - Prime or non-factorable rank counts raises error in main.cpp
 * 
 * TO-DO:
 * ----------------------
 * Read NumZRanks from param file to exactly match Abacus grid?
 */

// Global 2D Cartesian communicator
extern MPI_Comm comm_2d;

#endif
