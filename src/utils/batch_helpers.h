#ifndef HERMITIAN_BATCH_HELPERS_H
#define HERMITIAN_BATCH_HELPERS_H

#include <stdint.h>
#include "config.h"
#include "types.h"
#include <mpi.h>

// ====================================================================================
// Calculate how many Y-slices a given rank processes in a given batch
// Returns: 0 (no slices), 1 (self-conjugate), or 2 (conjugate pair)
// This allows each rank to independently calculate send/recv counts per batch

int get_rank_batch_slice_count(int target_rank, int batch_idx, int N, int num_ranks);

// ====================================================================================

// Get the specific Y-values that a rank is processing in a given batch
// Fills out_y_values array (must be pre-allocated with size >= 2)
// Sets out_count to number of Y-values (0, 1, or 2)

void get_rank_batch_y_values(int target_rank, int batch_idx, int N, int num_ranks,
                              int *out_y_values, int *out_count);
// ====================================================================================

// Calculate sendcounts, recvcounts, total send and recv counts, and send displacements for a given batch
// Allocates new arrays (caller must free them)
// Uses int64_t to avoid overflow: dest_region_size * my_batch_slice_count * narray can exceed INT_MAX for large N
// Note: recv displacements (rdispls) are computed in main.cpp using persistent recv_buffer cursors
// grid_x/grid_z are used for CPD-aligned dest bounds so send counts match main's decomposition.

void calculate_batch_send_recv_counts(
    int rank, int num_ranks, int N, int narray, int batch_idx,
    int my_batch_slice_count, int my_pencils,
    int grid_x, int grid_z, int cpd,
    int64_t **out_sendcounts, int64_t **out_sdispls,
    int64_t **out_recvcounts,
    int64_t *out_total_send, int64_t *out_total_recv);
// ====================================================================================

#endif 

