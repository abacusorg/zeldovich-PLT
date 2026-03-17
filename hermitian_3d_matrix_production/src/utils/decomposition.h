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

// CPD-aligned: grid_x * grid_z = num_ranks with grid_x|cpd and grid_z|cpd when possible.
// Uses N (PPD) and cpd from param file so rank x/z ranges fall on CPD slab boundaries.
void calculate_grid_factors_cpd_aligned(int num_ranks, int N, int cpd,
                                        int *grid_x_out, int *grid_z_out);

// Returns GridBounds struct  with x_start, x_end, z_start, z_end 
GridBounds get_grid_bounds(int dest, int N, int num_pencil_ranks);

// CPD-aligned bounds: works when grid_x/grid_z do not divide cpd (integer division).
// Uses integer division: slab_x_start(rx) = (rx*cpd)/grid_x, slab_x_end(rx) = ((rx+1)*cpd)/grid_x.
GridBounds get_grid_bounds_CPD_aligned(int dest, int N, int num_ranks,
                                       int grid_x, int grid_z, int cpd);

// Inverse: given CPD slab index s (x-direction), which rx owns it?
// rx = ((slab + 1) * grid_x - 1) / cpd
int slab_to_rx(int slab, int cpd, int grid_x);

// CPD-aligned extended bounds (core from get_grid_bounds_CPD_aligned, same padding logic).
ExtendedGridBounds get_extended_grid_bounds_CPD_aligned(int rank, int N, int num_ranks,
                                                        int grid_x, int grid_z, int cpd);

// If USE_X_PADDING=0: padded region equals core region (no padding)
ExtendedGridBounds get_extended_grid_bounds(int rank, int N, int num_ranks, int grid_x, int grid_z);

// Unified function to get padded or core bounds depending on USE_X_PADDING
GridBounds get_padded_bounds_simple(int dest, int N, int num_ranks);

#ifdef __cplusplus
}
#endif

#endif 

