#ifndef HERMITIAN_DECOMPOSITION_H
#define HERMITIAN_DECOMPOSITION_H

// ====================================================================================
// GRID DECOMPOSITION UTILITIES
// ====================================================================================
// Functions for calculating grid bounds & decompositions for
// the 2D (X,Z) pencil decomposition used for MPI comm. buffers
//
// Depends on: config.h, types.h
// ====================================================================================

#include "config.h"
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

// Calculates the grid decomposition dims from num_ranks
void calculate_grid_factors(int num_ranks, int *grid_x_out, int *grid_z_out);

// Returns GridBounds struct  with x_start, x_end, z_start, z_end 
GridBounds get_grid_bounds(int dest, int N, int num_pencil_ranks);

// If USE_X_PADDING=0: padded region equals core region (no padding)
ExtendedGridBounds get_extended_grid_bounds(int rank, int N, int num_ranks, int grid_x, int grid_z);

// Unified function to get padded or core bounds depending on USE_X_PADDING
GridBounds get_padded_bounds_simple(int dest, int N, int num_ranks);

#ifdef __cplusplus
}
#endif

#endif 

