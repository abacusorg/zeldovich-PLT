// ====================================================================================
// Z-STREAMING MODULE
// ====================================================================================

#include "z_streaming.h"
#include "../mpi_topology.h"
#include "../types.h"  // For ZSLAB macro
#include "../config.h"  // For DEBUG_PRINTS, SKIP_VERIFICATION
#include "../precision.h"  // For real_t, fabs_t, fmax_t
#include <stdio.h>
#include <math.h>    // For isinf, isnan
#include <mpi.h>

// Optional safety checks for streaming offsets (lightweight).
// Set to 1 temporarily when debugging packing/unpacking issues.
#ifndef VERIFY_STREAMING_OFFSETS
#define VERIFY_STREAMING_OFFSETS 0
#endif

// ====================================================================================
// Z-slab stream unpacking for Abcus-compatible output (?)
// Unpacks one Z-slab (cuz of mem) from the MPI receive buffer
// Applies 1D FFT along the Y-direction
// ** To-do: Put into format for zeldovich-PLT writing interface 

// ** NOTE: Batch-aware unpacking: The receive buffer is organized by batches, 
// so the function computes cumulative batch offsets to find the correct data location
// ====================================================================================

void z_streaming_unpack(
    int rank, int N, int narray,
    int z_global,                          // Which Z-slab to process
    GridBounds my_bounds,                  // My (X,Z) region bounds
    int my_pencils,                        // Total pencils in my region
    fftw_complex_t *recv_buffer,           // Source data (source-grouped)
    int64_t *recv_displs_src,              // Base offset per source
    int *src_total_slices,                 // Total slices per source (unused, kept for API)
    int *y_owner_src,                      // Y -> src mapping
    int *y_src_local_idx,                  // Y -> local_idx mapping (unused, kept for API)
    int *y_batch_idx,                      // Y -> batch mapping
    int *y_slice_idx_in_batch,             // Y -> slice_idx in batch mapping
    int **src_batch_slice_counts,          // [src][batch] -> slice count
    int global_max_batches,                // Total number of batches
    fftw_complex_t *local_z_slab,          // Destination buffer (one Z-slab)
    fftw_plan_t plan_1d_y)                 // FFT plan for Y-direction
{
    (void)rank;                // Unused in normal builds (used in debug checks)
    (void)global_max_batches;  // Unused but kept for API consistency
    (void)src_total_slices;    // Unused in normal builds (used only in VERIFY_STREAMING_OFFSETS)
    (void)y_src_local_idx;     // Unused in normal builds (used only in VERIFY_STREAMING_OFFSETS)
    
    int x_count = my_bounds.x_end - my_bounds.x_start;
    int z_count = my_bounds.z_end - my_bounds.z_start;
    int z_idx = z_global - my_bounds.z_start;
    
    // Sanity check
    if (z_global < my_bounds.z_start || z_global >= my_bounds.z_end) {
        fprintf(stderr, "[ERROR] Rank %d: z_global=%d out of bounds [%d,%d)\n",
                rank, z_global, my_bounds.z_start, my_bounds.z_end);
        MPI_Abort(comm_2d, 1);
    }
    
    // ========== UNPACKING: Extract this Z-slab from recv_buffer ==========
    // Loop over all arrays, all X in my region, all Y
    // v14: use PERIODIC_X for actual x_idx access
    // Fixed bug: Use batch and slice_idx information to compute correct offset
    // Data arrives as [array][slice_batch_local][pencil] per batch
    // We need to compute cumulative offset of all previous batches + current batch offset
    // Write to local_z_slab in [Array][X][Y] format (Y stride-1 for FFT, better cache locality)
    #pragma omp parallel for collapse(3)
    for (int array_idx = 0; array_idx < narray; array_idx++) {
        for (int x_idx = 0; x_idx < x_count; x_idx++) {
            for (int y = 0; y < N; y++) {

                // Use pre-built lookup arrays to find which rank owns
                // this Y-slice and which batch it came from
                int src = y_owner_src[y];
                int batch = y_batch_idx[y];
                int slice_idx = y_slice_idx_in_batch[y];
                
                // Calculate pencil index in full region: (x_local, z_local)
                int pencil_idx = x_idx * z_count + z_idx;
                
                // Compute cumulative offset of all previous batches from this source
                int64_t batch_offset = 0;
                for (int b = 0; b < batch; b++) {
                    int slices_in_batch = src_batch_slice_counts[src][b];
                    batch_offset += (int64_t)slices_in_batch * my_pencils * narray;
                }
                
                // Calculate position in recv_buffer
                // Actual layout: [src][array][batch_0_slice_0][pencil], [src][array][batch_0_slice_1][pencil], ...
                // Packing order per batch: [array][slice][pencil]
                // So for this batch: offset = batch_offset + array_idx * (slices_in_this_batch * my_pencils) + slice_idx * my_pencils + pencil_idx
                int slices_in_this_batch = src_batch_slice_counts[src][batch];
                int64_t recv_offset = recv_displs_src[src]
                                    + batch_offset
                                    + (int64_t)array_idx * slices_in_this_batch * my_pencils
                                    + (int64_t)slice_idx * my_pencils
                                    + pencil_idx;

#if VERIFY_STREAMING_OFFSETS
                // Lightweight safety checks:
                // 1) src_batch_slice_counts[src][batch] must not exceed src_total_slices[src].
                if (slices_in_this_batch < 0 || slices_in_this_batch > src_total_slices[src]) {
                    fprintf(stderr,
                            "[OFFSET-VERIFY] Rank %d: Invalid slices_in_this_batch=%d for src=%d, batch=%d "
                            "(src_total_slices=%d)\n",
                            rank, slices_in_this_batch, src, batch, src_total_slices[src]);
                    MPI_Abort(comm_2d, 1);
                }

                // 2) recv_offset must lie within this source's allocated region in recv_buffer.
                int64_t rel_src = recv_offset - recv_displs_src[src];
                int64_t max_src = (int64_t)src_total_slices[src] * (int64_t)my_pencils * (int64_t)narray;
                if (rel_src < 0 || rel_src >= max_src) {
                    fprintf(stderr,
                            "[OFFSET-VERIFY] Rank %d: recv_offset out of bounds for src=%d, batch=%d, slice_idx=%d, "
                            "array=%d, pencil=%d (x_idx=%d, z_idx=%d)\n"
                            "  rel_src=%lld, max_src=%lld, src_total_slices=%d, my_pencils=%d, narray=%d\n",
                            rank, src, batch, slice_idx, array_idx, pencil_idx, x_idx, z_idx,
                            (long long)rel_src, (long long)max_src,
                            src_total_slices[src], my_pencils, narray);
                    MPI_Abort(comm_2d, 1);
                }
#endif
                
                // V14: Apply periodic BC when writing to local_z_slab
                // x_idx is offset from my_bounds.x_start, which can be negative
                // No periodic needed here since ZSLAB uses local x_idx (not global X)
                // The periodic wrapping was already handled during packing
                
                // Extract to local_z_slab: [Array][X][Y] format (Y stride-1 for FFT)
                ZSLAB(array_idx, x_idx, y, N, narray, x_count)[0] = recv_buffer[recv_offset][0];
                ZSLAB(array_idx, x_idx, y, N, narray, x_count)[1] = recv_buffer[recv_offset][1];
            }
        }
    }
    
    // // ========== DEBUG: Check for Inf values at specific indices BEFORE 1D FFT ==========
    // // Known problematic locations: (j=133, k_local=126), (j=133, k_local=127), (j=134, k_local=0)
    // // Pattern: even z_global -> array 1 has Inf; odd z_global -> array 3 has Inf
    // #if 1  // Always enable for debugging
    // if (rank == 0) {
    //     // Check specific (y, x_idx) locations where Inf was found
    //     int debug_y_vals[] = {133, 133, 134};
    //     int debug_x_vals[] = {126, 127, 0};
    //     int num_debug_points = 3;
        
    //     int found_inf_before = 0;
    //     for (int d = 0; d < num_debug_points; d++) {
    //         int y = debug_y_vals[d];
    //         int x_idx = debug_x_vals[d];
            
    //         if (x_idx < x_count) {
    //             for (int array_idx = 0; array_idx < narray; array_idx++) {
    //                 real_t re = ZSLAB(array_idx, x_idx, y, N, narray, x_count)[0];
    //                 real_t im = ZSLAB(array_idx, x_idx, y, N, narray, x_count)[1];
                    
    //                 int is_inf = (isinf(re) || isinf(im) || isnan(re) || isnan(im));
    //                 if (is_inf) {
    //                     fprintf(stderr, "[INF-DEBUG Z=%d BEFORE FFT] array=%d x_idx=%d y=%d: re=%.6e im=%.6e\n",
    //                             z_global, array_idx, x_idx, y, (double)re, (double)im);
    //                     found_inf_before = 1;
    //                 }
    //             }
    //         }
    //     }
        
    //     // Also scan entire buffer for any Inf values
    //     int inf_count_before[4] = {0, 0, 0, 0};
    //     for (int array_idx = 0; array_idx < narray; array_idx++) {
    //         for (int x_idx = 0; x_idx < x_count; x_idx++) {
    //             for (int y = 0; y < N; y++) {
    //                 real_t re = ZSLAB(array_idx, x_idx, y, N, narray, x_count)[0];
    //                 real_t im = ZSLAB(array_idx, x_idx, y, N, narray, x_count)[1];
    //                 if (isinf(re) || isinf(im) || isnan(re) || isnan(im)) {
    //                     inf_count_before[array_idx]++;
    //                 }
    //             }
    //         }
    //     }
        
    //     int total_inf_before = 0;
    //     for (int a = 0; a < narray; a++) total_inf_before += inf_count_before[a];
    //     if (total_inf_before > 0) {
    //         fprintf(stderr, "[INF-DEBUG Z=%d BEFORE FFT] Total Inf/NaN count per array: ", z_global);
    //         for (int a = 0; a < narray; a++) {
    //             fprintf(stderr, "A%d=%d ", a, inf_count_before[a]);
    //         }
    //         fprintf(stderr, "\n");
    //     }
    // }
    // #endif
    
    #if DEBUG_PRINTS && !SKIP_VERIFICATION
    if (rank == 0 && z_global == my_bounds.z_start) {
        real_t debug_max_real_before[4] = {0, 0, 0, 0};
        real_t debug_max_imag_before[4] = {0, 0, 0, 0};
        for (int array_idx = 0; array_idx < narray; array_idx++) {
            for (int x_idx = 0; x_idx < x_count; x_idx++) {
                for (int y = 0; y < N; y++) {
                    double re = fabs_t(ZSLAB(array_idx, x_idx, y, N, narray, x_count)[0]);
                    double im = fabs_t(ZSLAB(array_idx, x_idx, y, N, narray, x_count)[1]);
                    debug_max_real_before[array_idx] = fmax_t(debug_max_real_before[array_idx], re);
                    debug_max_imag_before[array_idx] = fmax_t(debug_max_imag_before[array_idx], im);
                }
            }
        }
        printf("[DEBUG Z=%d] Before 1D FFT (after unpack) - Max imag per array: ", z_global);
        for (int a = 0; a < narray; a++) {
            printf("A%d=%.3e ", a, debug_max_imag_before[a]);
        }
        printf("\n");
    }
    #endif
    
    // ========== 1D FFT: Apply along Y-direction for each (Array, X) ==========
    // In [Array][X][Y] format, Y is stride-1 for fixed (array_idx, x_idx)
    #pragma omp parallel for collapse(2) schedule(static)
    for (int array_idx = 0; array_idx < narray; array_idx++) {
        for (int x_idx = 0; x_idx < x_count; x_idx++) {
            fftw_complex_t *y_data = &ZSLAB(array_idx, x_idx, 0, N, narray, x_count);
            FFTW_EXECUTE_DFT(plan_1d_y, y_data, y_data);
        }
    }
    
    // // ========== DEBUG: Extract real(FFT(D + i*F)) or real(FFT(D)) for comparison ==========
    // // After 3D FFT, Array 0 contains either FFT(D + i*F) or FFT(D) depending on just_density
    // // Print real parts for test coordinates to compare between runs
    // // Enable for N <= 16 to match zeldovich code debug output
    // // Match zeldovich.cpp format: [REAL-FFT-DEBUG] N=%d just_density=%d Z=%d (x,y)=(%d,%d): real(Array0)=%.10e imag(Array0)=%.10e
    // // Print for all Z slabs (not just z_global == 0) to match previous behavior
    // // Only print from rank 0 to avoid potential issues with concurrent writes
    // #if DEBUG_RNG_CONSISTENCY
    // if (rank == 0 && N <= 16) {
    //     int just_density_flag = (narray == 1) ? 1 : 0;  // narray == 1 means just_density mode
    //     int max_test_coord = (N <= 16) ? N : 4;
    //     for (int x_idx = 0; x_idx < x_count && x_idx < max_test_coord; x_idx++) {
    //         int x_global = my_bounds.x_start + x_idx;
    //         for (int y = 0; y < max_test_coord && y < N; y++) {
    //             // Array 0: Contains FFT(D + i*F) when just_density=false, or FFT(D) when just_density=true
    //             // Use ZSLAB: (array_idx, x_idx, y, N, narray, x_count)
    //             real_t re = ZSLAB(0, x_idx, y, N, narray, x_count)[0];
    //             real_t im = ZSLAB(0, x_idx, y, N, narray, x_count)[1];
    //             fprintf(stderr, "[REAL-FFT-DEBUG] N=%d just_density=%d Z=%d (x,y)=(%d,%d): "
    //                     "real(Array0)=%.10e imag(Array0)=%.10e\n",
    //                     N, just_density_flag, z_global, x_global, y, (double)re, (double)im);
    //         }
    //     }
    //     fflush(stderr);
    // }
    // #endif
    
    // ========== DEBUG: Check for Inf values at specific indices AFTER 1D FFT ==========
    #if 1  // Always enable for debugging
    if (rank == 0) {
        // Check specific (y, x_idx) locations where Inf was found
        int debug_y_vals[] = {133, 133, 134};
        int debug_x_vals[] = {126, 127, 0};
        int num_debug_points = 3;
        
        for (int d = 0; d < num_debug_points; d++) {
            int y = debug_y_vals[d];
            int x_idx = debug_x_vals[d];
            
            if (x_idx < x_count) {
                for (int array_idx = 0; array_idx < narray; array_idx++) {
                    real_t re = ZSLAB(array_idx, x_idx, y, N, narray, x_count)[0];
                    real_t im = ZSLAB(array_idx, x_idx, y, N, narray, x_count)[1];
                    
                    int is_inf = (isinf(re) || isinf(im) || isnan(re) || isnan(im));
                    if (is_inf) {
                        fprintf(stderr, "[INF-DEBUG Z=%d AFTER FFT] array=%d x_idx=%d y=%d: re=%.6e im=%.6e\n",
                                z_global, array_idx, x_idx, y, (double)re, (double)im);
                    }
                }
            }
        }
        
        // Also scan entire buffer for any Inf values
        int inf_count_after[4] = {0, 0, 0, 0};
        for (int array_idx = 0; array_idx < narray; array_idx++) {
            for (int x_idx = 0; x_idx < x_count; x_idx++) {
                for (int y = 0; y < N; y++) {
                    real_t re = ZSLAB(array_idx, x_idx, y, N, narray, x_count)[0];
                    real_t im = ZSLAB(array_idx, x_idx, y, N, narray, x_count)[1];
                    if (isinf(re) || isinf(im) || isnan(re) || isnan(im)) {
                        inf_count_after[array_idx]++;
                    }
                }
            }
        }
        
        int total_inf_after = 0;
        for (int a = 0; a < narray; a++) total_inf_after += inf_count_after[a];
        if (total_inf_after > 0) {
            fprintf(stderr, "[INF-DEBUG Z=%d AFTER FFT] Total Inf/NaN count per array: ", z_global);
            for (int a = 0; a < narray; a++) {
                fprintf(stderr, "A%d=%d ", a, inf_count_after[a]);
            }
            fprintf(stderr, "\n");
        }
    }
    #endif
    
    // local_z_slab now contains FFT'd data for this Z-slab, ready for writing
}

