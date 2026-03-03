#ifndef MPI_EXCHANGE_H
#define MPI_EXCHANGE_H

// ====================================================================================
// MPI COMMUNICATION MODULE
// ====================================================================================
// MPI communication for Y-slice distribution and z-slab assembly
//
// Depends on: config.h, precision.h, types.h, utils/decomposition.h
// ====================================================================================

#include <stdint.h>
#include "../config.h"
#include "../precision.h"
#include "../types.h"
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

// ====================================================================================
// FUNCTIONS
// ====================================================================================

// Exchange metadata between all ranks: num_my_slices and y_global_map
// MULTI-BATCH VERSION: Each rank shares the COMPLETE list of Y-values
// it will process across ALL batches (not just one batch).
//
// Method A (METADATA_EXCHANGE_METHOD=0): Each rank calculates + MPI_Allgatherv
// Method B (METADATA_EXCHANGE_METHOD=1): Rank 0 calculates + MPI_Bcast
//
// Outputs:
// - out_all_num_my_slices: Array of slice counts per rank [num_y_ranks]
// - out_all_y_global_maps: Array of Y-index arrays per rank [num_y_ranks][variable]
void exchange_metadata(int rank, int num_y_ranks, int N,
                      int num_my_slices, int *y_global_map,
                      int **out_all_num_my_slices, int ***out_all_y_global_maps);

// Pack my Y-slices into send buffer for MPI_Ialltoallv
// For each destination rank, extract the (X,Z) region they own from ALL my Y-slices
// Packs all narray arrays for each slice
//
// MULTI-BATCH COMPATIBLE: Works unchanged for multi-batch processing.
// - For single-batch: num_my_slices = 1 or 2 (all slices for this rank)
// - For multi-batch: num_my_slices = 1 or 2 (slices for THIS batch only)
//
// Packing order: [array][slice][x][z] per destination
// Uses PERIODIC_X/Z macros for boundary handling (V14)
void pack_slices_to_send_buffer(
    int rank, int num_ranks, int N, int narray,
    fftw_complex_t *local_y_slices,  // Flat buffer with all arrays
    int num_my_slices, int *y_global_map,
    fftw_complex_t *send_buffer, int64_t *sendcounts, int64_t *sdispls);

// Unpack received data into local pencils
// Each pencil is a Y-column for one (X,Z) point in my owned region
// Data arrives from all Y-slice ranks, organized by source rank
// Unpacks ALL narray arrays for each pencil
//
// Memory layout: local_pencils[Pencil][Array][Y]
// Packing order in recv_buffer: [array][slice][x][z] per source rank
void unpack_recv_buffer_to_pencils(
    int rank, int num_ranks, int N, int narray,
    fftw_complex_t *recv_buffer, int *recvcounts, int *rdispls,
    int **all_y_global_maps, int *all_num_my_slices,
    fftw_complex_t *local_pencils);

#ifdef __cplusplus
}
#endif

#endif

