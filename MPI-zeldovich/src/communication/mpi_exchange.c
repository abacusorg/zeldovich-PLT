#include "mpi_exchange.h"
#include "../mpi_topology.h"
#include "../utils/decomposition.h"
#include "../utils/verification.h"
#include "../types.h" 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

// ====================================================================================
//  Exchange metadata so that each rank knows which Y-slices it is responsible for
//  and which Y-indices it is responsible for during unpacking
// ====================================================================================

void exchange_metadata(int rank, int num_y_ranks, int N,
                      int num_my_slices, int *y_global_map,
                      int **out_all_num_my_slices, int ***out_all_y_global_maps)
{
    // Suppress unused parameter warnings for Method B
    (void)num_my_slices;
    (void)y_global_map;
    
    #if METADATA_EXCHANGE_METHOD == 0
    // ========== METHOD A: Each rank calculates + Allgatherv ==========
    // Unused!! - kept for future reference. See Method B below.
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[METADATA] Using Method A: Each rank calculates + MPI_Allgatherv\n");
    }
    
    // Step 1: Share num_my_slices using MPI_Allgather
    // array containing the number of Y-slices for each rank
    int *all_num_my_slices = (int*)malloc(sizeof(int) * num_y_ranks);
    MPI_Allgather(&num_my_slices, 1, MPI_INT, 
                  all_num_my_slices, 1, MPI_INT, comm_2d);
    
    // Step 2: Share y_global_map using MPI_Allgatherv
    // Calculate displacements and total size
    int *recvcounts = (int*)malloc(sizeof(int) * num_y_ranks);
    int *displs = (int*)malloc(sizeof(int) * num_y_ranks);
    int total_slices = 0;
    for (int i = 0; i < num_y_ranks; i++) {
        recvcounts[i] = all_num_my_slices[i]; // how many Y-slices to receive from each rank
        displs[i] = total_slices;             // cumulative offset of Y-slices for each rank
        total_slices += all_num_my_slices[i]; 
    }
    
    // Allocate flat buffer for all y_global_maps
    int *all_y_maps_flat = (int*)malloc(sizeof(int) * total_slices);
    
    MPI_Allgatherv(y_global_map, num_my_slices, MPI_INT,
                   all_y_maps_flat, recvcounts, displs, MPI_INT, comm_2d);
    
    // Step 3: Rebuild pointer array (allocate separately for each rank to avoid pointer issues)
    int **all_y_global_maps = (int**)malloc(sizeof(int*) * num_y_ranks);
    for (int i = 0; i < num_y_ranks; i++) {
        all_y_global_maps[i] = (int*)malloc(sizeof(int) * all_num_my_slices[i]);
        memcpy(all_y_global_maps[i], &all_y_maps_flat[displs[i]],
               sizeof(int) * all_num_my_slices[i]);
    }
    
    free(recvcounts);
    free(displs);
    free(all_y_maps_flat);  // Now safe to free
    
    *out_all_num_my_slices = all_num_my_slices;
    *out_all_y_global_maps = all_y_global_maps;
    
    #else
    // ========== METHOD B: Rank 0 calculates + Broadcast ==========
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[METADATA] Using Method B: Rank 0 calculates + MPI_Bcast\n");
    }
    
    // All ranks: allocate arrays for # of Y-slices & Y-indices it is responsible for
    int *all_num_my_slices = (int*)malloc(sizeof(int) * num_y_ranks); // number of Y-slices for each rank
    int **all_y_global_maps = (int**)malloc(sizeof(int*) * num_y_ranks); // Y-indices for each rank
    
    // MULTI-BATCH: Each rank can process multiple pairs
    if (rank == 0) {
        // Rank 0: Calculate metadata for all ranks using same distribution as Phase 2
        int total_pairs = N/2 + 1;
        int pairs_per_rank = total_pairs / num_y_ranks;
        int remainder = total_pairs % num_y_ranks;
        
        for (int i = 0; i < num_y_ranks; i++) {
            // Calculate this rank's pair assignment
            int rank_num_pairs;
            int rank_pair_start;
            
            if (i < remainder) {
                // This rank gets one more pair than the average (remainder)
                rank_num_pairs = pairs_per_rank + 1;
                rank_pair_start = i * (pairs_per_rank + 1);
            } else {
                // This rank gets the average number of pairs
                rank_num_pairs = pairs_per_rank;
                rank_pair_start = remainder * (pairs_per_rank + 1) + (i - remainder) * pairs_per_rank;
            }
            
            // Count total Y-values for this rank
            int total_y_values = 0;
            for (int j = 0; j < rank_num_pairs; j++) {
                int pair_idx = rank_pair_start + j;
                int y_primary = pair_idx;
                bool is_self_conjugate = (y_primary == 0 || y_primary == N/2);
                total_y_values += is_self_conjugate ? 1 : 2;
            }
            
            all_num_my_slices[i] = total_y_values;
            all_y_global_maps[i] = (int*)malloc(sizeof(int) * total_y_values);
            
            // Populate all_y_global_maps[i] with y-indices that rank i will own
            int idx = 0;
            for (int j = 0; j < rank_num_pairs; j++) {
                int pair_idx = rank_pair_start + j;
                int y_primary = pair_idx;
                int y_mirror = (y_primary == 0 || y_primary == N/2) ? y_primary : N - y_primary;
                bool is_self_conjugate = (y_primary == y_mirror);
                
                all_y_global_maps[i][idx++] = y_primary;
                if (!is_self_conjugate) {
                    all_y_global_maps[i][idx++] = y_mirror;
                }
            }
        }
    }
    
    // Broadcast all_num_my_slices
    MPI_Bcast(all_num_my_slices, num_y_ranks, MPI_INT, 0, comm_2d);
    
    // Broadcast each y_global_map
    for (int i = 0; i < num_y_ranks; i++) {
        if (rank != 0) {
            all_y_global_maps[i] = (int*)malloc(sizeof(int) * all_num_my_slices[i]);
        }
        MPI_Bcast(all_y_global_maps[i], all_num_my_slices[i], MPI_INT, 0, comm_2d);
    }
    
    // Assign outputs to the function args ("return values")
    *out_all_num_my_slices = all_num_my_slices;
    *out_all_y_global_maps = all_y_global_maps;
    
    #endif
    
    // Debug: Print metadata on rank 0
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[METADATA] Metadata exchange complete:\n");
        for (int i = 0; i < (num_y_ranks < 5 ? num_y_ranks : 5); i++) {
            printf("  Rank %d: %d slices, Y-indices = [", i, (*out_all_num_my_slices)[i]);
            for (int j = 0; j < (*out_all_num_my_slices)[i]; j++) {
                printf("%d%s", (*out_all_y_global_maps)[i][j], 
                       (j < (*out_all_num_my_slices)[i] - 1) ? ", " : "");
            }
            printf("]\n");
        }
        if (num_y_ranks > 5) {
            printf("  ... (%d more ranks)\n", num_y_ranks - 5);
        }
    }
}

// ====================================================================================
//  Pack the Y-slices into send buffer in rank-contiguous manner for MPI comm
//  For each destination rank, extracts the (X,Z) region it owns from all local Y-slices
// ====================================================================================

void pack_slices_to_send_buffer(
    int rank, int num_ranks, int N, int narray,
    fftw_complex_t *local_y_slices,
    int num_my_slices, int *y_global_map,
    fftw_complex_t *send_buffer, int64_t *sendcounts, int64_t *sdispls)
{
    (void)rank;  // Unused, kept for consistency
    (void)y_global_map;  // Unused currently, kept for future use
    
    // Phase 1 (sequential): Compute sdispls, sendcounts, and per-dest bounds.
    // Use int64_t to avoid overflow: region_size * num_my_slices * narray can exceed INT_MAX for large N.
    int64_t total_send_size = 0;
    int64_t offset = 0;
    GridBounds *bounds_arr = (GridBounds *)malloc((size_t)num_ranks * sizeof(GridBounds));
    int *z_count_arr = (int *)malloc((size_t)num_ranks * sizeof(int));
    if (bounds_arr == NULL || z_count_arr == NULL) {
        fprintf(stderr, "[PACK ERROR] Rank %d: Failed to allocate bounds/z_count arrays\n", rank);
        MPI_Abort(comm_2d, 1);
    }
    
    for (int dest = 0; dest < num_ranks; dest++) {
        sdispls[dest] = offset;
        
        GridBounds bounds = get_padded_bounds_simple(dest, N, num_ranks);
        bounds_arr[dest] = bounds;
        int64_t region_size = (int64_t)(bounds.x_end - bounds.x_start) * (int64_t)(bounds.z_end - bounds.z_start);
        
        sendcounts[dest] = region_size * (int64_t)num_my_slices * (int64_t)narray;
        z_count_arr[dest] = bounds.z_end - bounds.z_start;
        
        offset += sendcounts[dest];
    }
    total_send_size = offset;
    
    // Phase 2: One parallel region, one thread team; for each dest, workshare the inner loops.
    // This replaces 32 fork/joins with 1 fork + 32 lightweight barriers (~20x less overhead).
    #pragma omp parallel
    for (int dest = 0; dest < num_ranks; dest++) {
        int64_t dest_offset = sdispls[dest];
        GridBounds bounds = bounds_arr[dest];
        int z_count = z_count_arr[dest];
        int64_t region_size = (int64_t)(bounds.x_end - bounds.x_start) * (int64_t)(bounds.z_end - bounds.z_start);
        
        #pragma omp for collapse(4) nowait
        for (int array_idx = 0; array_idx < narray; array_idx++) {
            for (int slice_idx = 0; slice_idx < num_my_slices; slice_idx++) {
                for (int x = bounds.x_start; x < bounds.x_end; x++) {
                    for (int z = bounds.z_start; z < bounds.z_end; z++) {
                        int x_actual = PERIODIC_X(x, N);
                        int z_actual = PERIODIC_Z(z, N);
                        
                        int local_x = x - bounds.x_start;
                        int local_z = z - bounds.z_start;
                        
                        int64_t pack_idx = (int64_t)array_idx * num_my_slices * region_size
                                         + slice_idx * region_size
                                         + local_x * z_count + local_z;

                        int64_t buffer_idx = dest_offset + pack_idx;
                        if (buffer_idx < 0 || buffer_idx >= total_send_size) {
                            fprintf(stderr, "[PACK ERROR] Rank %d: buffer_idx=%ld out of bounds [0, %lld) for dest=%d\n",
                                   rank, (long)buffer_idx, (long long)total_send_size, dest);
                            fprintf(stderr, "  array_idx=%d, slice_idx=%d, x=%d, z=%d, local_x=%d, local_z=%d\n",
                                   array_idx, slice_idx, x, z, local_x, local_z);
                            fprintf(stderr, "  pack_idx=%ld, offset=%d, region_size=%d, sendcounts[dest]=%d\n",
                                   (long)pack_idx, dest_offset, region_size, sendcounts[dest]);
                            MPI_Abort(comm_2d, 1);
                        }
                        if (pack_idx < 0 || pack_idx >= sendcounts[dest]) {
                            fprintf(stderr, "[PACK ERROR] Rank %d: pack_idx=%ld out of dest region [0, %lld) for dest=%d\n",
                                   rank, (long)pack_idx, (long long)sendcounts[dest], dest);
                            MPI_Abort(comm_2d, 1);
                        }

                        send_buffer[dest_offset + pack_idx][0] = 
                            Y_SLICE(slice_idx, array_idx, x_actual, z_actual, N, narray)[0];
                        send_buffer[dest_offset + pack_idx][1] = 
                            Y_SLICE(slice_idx, array_idx, x_actual, z_actual, N, narray)[1];
                    }
                }
            }
        }
    }
    
    free(bounds_arr);
    free(z_count_arr);
    
    if (DEBUG_PRINTS && rank < 3) {
        printf("[PACK] Rank %d: Packed %lld total elements to send_buffer\n", rank, (long long)offset);
        printf("       Sending to ranks: ");
        for (int dest = 0; dest < (num_ranks < 4 ? num_ranks : 4); dest++) {
            printf("%d elements to rank %d%s", sendcounts[dest], dest, 
                   (dest < num_ranks - 1) ? ", " : "");
        }
        if (num_ranks > 4) printf("...");
        printf("\n");
    }
}

// ====================================================================================
//  Unpack recv buffer (different for each rank) into pencils, 
//  keepin track of global y-indices to maintain correct Hermitian symmetry
// ====================================================================================

void unpack_recv_buffer_to_pencils(
    int rank, int num_ranks, int N, int narray,
    fftw_complex_t *recv_buffer, int *recvcounts, int *rdispls,
    int **all_y_global_maps, int *all_num_my_slices,
    fftw_complex_t *local_pencils)
{
    (void)recvcounts;  // Unused, kept for future extensibility?
    
    // V13: Use PADDED bounds for unpacking (includes overlap regions)
    GridBounds my_bounds = get_padded_bounds_simple(rank, N, num_ranks);
    int my_pencils = (my_bounds.x_end - my_bounds.x_start) * (my_bounds.z_end - my_bounds.z_start);
    int z_count = my_bounds.z_end - my_bounds.z_start; // number of z points in my owned chunk
    
    // Initialize pencil tracking for verification (track for all arrays)
    #if VERIFY_Y_FILLED
    char *y_filled = (char*)calloc(my_pencils * narray * N, sizeof(char));
    #endif
    
    // Unpack data from each source rank
    // Use single parallel region to avoid ~8000 fork/join operations
    // (64 sources x 4 arrays x ~32 slices = ~8000 combinations)
    #pragma omp parallel
    for (int src = 0; src < num_ranks; src++) {
        // Get offset of source's data in the recv buffer
        int src_offset = rdispls[src];

        // Get the number of Y-slices from this source
        int src_num_slices = all_num_my_slices[src];
        
        // STAGE 5: Unpack all narray arrays for each source
        for (int array_idx = 0; array_idx < narray; array_idx++) {
            for (int slice_idx = 0; slice_idx < src_num_slices; slice_idx++) {
                // Get global Y index of slice
                int y_global = all_y_global_maps[src][slice_idx]; // <-- keeps track of global pos of y-slice!!
                
                // OpenMP: Workshare over (X,Z) points in my owned chunk
                // Using existing thread team (no fork/join overhead)
                #pragma omp for collapse(2) nowait
                for (int x = my_bounds.x_start; x < my_bounds.x_end; x++) {
                    for (int z = my_bounds.z_start; z < my_bounds.z_end; z++) {
                        // Calculate unpack and pencil indices for this thread's iteration
                        int local_x = x - my_bounds.x_start;
                        int local_z = z - my_bounds.z_start;
                        
                        int pencil_idx = local_x * z_count + local_z;
                        
                        // Unpack index matches pack index: [array][slice][pencil]
                        // Layout: array_idx * (src_num_slices * my_pencils) <-- each array has # src_num_slices
                        //        + slice_idx * my_pencils <-- each slice contains # my_pencils
                        //        + pencil_idx
                        int64_t unpack_idx = (int64_t)array_idx * src_num_slices * my_pencils
                                           + slice_idx * my_pencils
                                           + pencil_idx;

                        // Use PENCIL macro to unpack into correct position
                        // Unpack the data into the local pencils
                        PENCIL(pencil_idx, array_idx, y_global, N, narray)[0] = 
                            recv_buffer[src_offset + unpack_idx][0];
                        PENCIL(pencil_idx, array_idx, y_global, N, narray)[1] = 
                            recv_buffer[src_offset + unpack_idx][1];
                        
                        #if VERIFY_Y_FILLED
                        y_filled[pencil_idx * narray * N + array_idx * N + y_global] = 1;
                        #endif
                    }
                }
            }
        }
    }
    
    if (DEBUG_PRINTS && rank < 3) {
        printf("[UNPACK] Rank %d: Unpacked %d pencils (Y-columns), each with %d arrays, %d Y-values\n", 
               rank, my_pencils, narray, N);
    }
    
    // Verify completeness (updated for narray)
    #if VERIFY_Y_FILLED
    verify_pencil_completeness_with_flags(y_filled, my_pencils, N, narray, rank);
    free(y_filled);
    #endif
}

