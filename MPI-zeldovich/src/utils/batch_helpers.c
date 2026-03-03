// ====================================================================================
// BATCH PROCESSING HELPERS
// ====================================================================================

#include "utils/batch_helpers.h"
#include "utils/decomposition.h"
#include <stdlib.h>
#include <stdbool.h>

// ====================================================================================
// Maps from rank's local batch idx to the global y pair idx, based on rank division
// Each rank compute which Y-slices it's processing in a batch, w/o communication
// ====================================================================================

int get_rank_batch_slice_count(int target_rank, int batch_idx, int N, int num_ranks)
{
    // Calculate total pairs and distribution (same logic as main)
    int total_pairs = N/2 + 1;
    int pairs_per_rank = total_pairs / num_ranks;
    int remainder = total_pairs % num_ranks;
    
    // Calculate how many pairs this target_rank has
    int rank_num_pairs;
    int rank_pair_start;
    
    if (target_rank < remainder) {
        rank_num_pairs = pairs_per_rank + 1; // extra pair for remainder ranks
        rank_pair_start = target_rank * (pairs_per_rank + 1);
    } else {
        rank_num_pairs = pairs_per_rank;     // regular pairs for non-remainder ranks
        rank_pair_start = remainder * (pairs_per_rank + 1) + (target_rank - remainder) * pairs_per_rank;
    }
    
    // Check if batch_idx is within this rank's range
    if (batch_idx >= rank_num_pairs) {
        return 0;  // This rank has no data for this batch
    }
    
    // Calculate which pair this batch corresponds to
    int pair_idx = rank_pair_start + batch_idx;
    int y_primary = pair_idx;
    int y_mirror = (y_primary == 0 || y_primary == N/2) ? y_primary : N - y_primary;
    bool is_self_conjugate = (y_primary == y_mirror);
    
    return is_self_conjugate ? 1 : 2;
}

// ====================================================================================
// What y values is the rank processing in a given batch?
// Example 
// N=8, num_ranks=3:
// Rank 0, batch 0: y_values = [0], count = 1 (self-conjugate Y=0)
// Rank 0, batch 1: y_values = [1, 7], count = 2 (conjugate pair)
// Rank 1, batch 0: y_values = [2, 6], count = 2 (conjugate pair)
// ...
// ====================================================================================

void get_rank_batch_y_values(int target_rank, int batch_idx, int N, int num_ranks,
                              int *out_y_values, int *out_count)
{
    // Calculate total pairs and distribution (same logic as main)
    int total_pairs = N/2 + 1;
    int pairs_per_rank = total_pairs / num_ranks;
    int remainder = total_pairs % num_ranks;
    
    // Calculate how many pairs this target_rank has
    int rank_num_pairs;
    int rank_pair_start;
    
    if (target_rank < remainder) {
        rank_num_pairs = pairs_per_rank + 1;
        rank_pair_start = target_rank * (pairs_per_rank + 1);
    } else {
        rank_num_pairs = pairs_per_rank;
        rank_pair_start = remainder * (pairs_per_rank + 1) + (target_rank - remainder) * pairs_per_rank;
    }
    
    // Check if batch_idx is within this rank's range
    if (batch_idx >= rank_num_pairs) {
        *out_count = 0;  // This rank has no data for this batch
        return;
    }
    
    // Calculate which pair this batch corresponds to
    int pair_idx = rank_pair_start + batch_idx;
    int y_primary = pair_idx;
    int y_mirror = (y_primary == 0 || y_primary == N/2) ? y_primary : N - y_primary;
    bool is_self_conjugate = (y_primary == y_mirror);
    
    // Fill output
    out_y_values[0] = y_primary;
    if (is_self_conjugate) {
        *out_count = 1;
    } else {
        out_y_values[1] = y_mirror;
        *out_count = 2;
    }
}

// ====================================================================================
// Calculates MPI send/recv arrays for MPI_Alltoallv_c
// Allocates and returns: sendcounts, sdispls, recvcounts
// Note: rdispls are computed in main.cpp using persistent recv_buffer cursors
// Uses padded bounds for destination regions (includes overlap regions)
// ====================================================================================

void calculate_batch_send_recv_counts(
    int rank, int num_ranks, int N, int narray, int batch_idx,
    int my_batch_slice_count, int my_pencils,
    int64_t **out_sendcounts, int64_t **out_sdispls,
    int64_t **out_recvcounts,
    int64_t *out_total_send, int64_t *out_total_recv)
{
    (void)rank;  // Unused but kept for compatibility
    
    // Allocate arrays, each array's size = num_ranks
    int64_t *sendcounts = (int64_t*)malloc(sizeof(int64_t) * num_ranks);
    int64_t *sdispls = (int64_t*)malloc(sizeof(int64_t) * num_ranks);
    int64_t *recvcounts = (int64_t*)malloc(sizeof(int64_t) * num_ranks);
    
    // SEND COUNTS: I send my batch's Y-slices to all ranks (each gets their (X,Z) region)
    // Use int64_t to avoid overflow: dest_region_size * my_batch_slice_count * narray can exceed INT_MAX
    int64_t total_send = 0;
    for (int dest = 0; dest < num_ranks; dest++) {
        // V13: Use (optional) PADDED bounds for destination (includes overlap regions)
        GridBounds dest_bounds = get_padded_bounds_simple(dest, N, num_ranks);
        int64_t dest_region_size = (int64_t)(dest_bounds.x_end - dest_bounds.x_start) *
                                  (int64_t)(dest_bounds.z_end - dest_bounds.z_start);
        
        // Send: my Y-slices * dest's (X,Z) region * narray
        sendcounts[dest] = dest_region_size * (int64_t)my_batch_slice_count * (int64_t)narray;
        sdispls[dest] = total_send;
        total_send += sendcounts[dest];
    }
    
    // RECV COUNTS: I receive all ranks' batch Y-slices for my (X,Z) region
    int64_t total_recv = 0;
    for (int src = 0; src < num_ranks; src++) {
        // How many Y-slices is src sending in this batch? (0, 1, 2)
        int src_batch_slice_count = get_rank_batch_slice_count(src, batch_idx, N, num_ranks);
        
        // Receive: my pencils * src's Y-slices * narray
        recvcounts[src] = (int64_t)my_pencils * (int64_t)src_batch_slice_count * (int64_t)narray;
        total_recv += recvcounts[src];
    }
    
    // Return results
    *out_sendcounts = sendcounts;
    *out_sdispls = sdispls;
    *out_recvcounts = recvcounts;
    *out_total_send = total_send;
    *out_total_recv = total_recv;
}

