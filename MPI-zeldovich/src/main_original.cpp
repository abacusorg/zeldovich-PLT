// ====================================================================================
// VERSION 14: PERIODIC BOUNDARY CONDITIONS FOR PADDING
// ====================================================================================
// - Real MPI implementation (not OpenMP mock)
// - Each MPI rank generates multiple Y-slice pairs across batches
// - OpenMP used WITHIN each rank to parallelize X-Z loops during generation
// - 2D grid decomposition for redistribution to pencil layout
// - **NEW in v14**: PERIODIC BOUNDARY CONDITIONS FOR PADDING
//   * Padding uses periodic boundaries (wrap-around) instead of clamping
//   * Example: rank 0: [0,100] core → [-10,110] padded (wraps from right side!)
//             rank 1: [100,200] core → [90,210] padded (normal overlap)
//   * Physics-correct for cosmological simulations (periodic box)
//   * X=-10 wraps to X=N-10, X=N+5 wraps to X=5 (modulo arithmetic)
//   * Efficient indexing: PERIODIC_X macro with fast-path for common case
// - **RETAINED from v13**: OVERLAPPING X-REGIONS FOR ABACUS GRID READING
//   * Each rank stores extra X-values beyond its core region (X_PADDING on each side)
//   * Handles Abacus's imperfect grid reading (may access X outside strict bounds)
//   * Creates redundant storage at boundaries for seamless Abacus integration
//   * File output includes padded regions (some X-values written by multiple ranks)
// - **RETAINED from v12**: Z-SLAB STREAMING (ZELDOVICH-COMPATIBLE OUTPUT)
//   * Process one Z-slab at a time instead of X-rows
//   * Output format: [Z][Array][Y][X] matches Zeldovich's AZYX layout
//   * Each file contains all (X,Y) for a given Z-slab
//   * Direct compatibility with WriteParticlesSlab interface
// - **RETAINED from v11**: PERSISTENT RECV_BUFFER WITH SOURCE-GROUPED LAYOUT
//   * recv_buffer layout: [src_rank][y_idx_for_src][pencil_idx][array_idx]
//   * Direct receive into persistent buffer (NO per-batch recv_buffer_batch)
//   * Per-source write cursors track progress across batches
//   * Eliminates ~64 GB transient allocation per batch (single precision, N=32K)
// - Hermitian symmetry verification at multiple stages
// - Final result: Real-space data (imaginary parts ≈ 0)
// ====================================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>  // For INT_MAX
#include <mpi.h>
#include <omp.h>  // For hybrid MPI+OpenMP within each rank
#include <fftw3.h>
#include <sys/stat.h>
#include <errno.h>

// Include PCG RNG and STimer
#include "pcg-rng/pcg_random.hpp"
#ifdef __cplusplus
extern "C" {
#endif
#include "STimer.h"
#ifdef __cplusplus
}
#endif

// ====================================================================================
// CONFIGURATION AND TYPES
// ====================================================================================
// All compile-time flags and constants are now in config.h
// All type definitions and precision handling are in precision.h and types.h
#include "config.h"
#include "precision.h"
#include "types.h"

// ====================================================================================
// UTILITY MODULES
// ====================================================================================
// Utility functions extracted to separate modules
#include "utils/printing.h"
#include "utils/verification.h"
#include "utils/decomposition.h"
#include "utils/batch_helpers.h"

// ====================================================================================
// CORE MODULES
// ====================================================================================
// Core algorithm modules extracted in Phase 5
#include "fft/fft_setup.h"
#include "generation/hermitian_generation.h"
#include "communication/mpi_exchange.h"
#include "streaming/z_streaming.h"

// ====================================================================================
// PCG RNG MODULE
// ====================================================================================

// NOTE: Using local RNG module (src/utils/rng.c)
// The RNG module uses pcg64 with advance() for sequential initialization
// This matches the original approach used in other versions

#include "utils/rng.h"
#include "utils/power_spectrum.h"
#include "utils/zeldovich_wrapper.h"  // For zeldovich-PLT integration

// Note: Precision types (real_t, fftw_complex_t, etc.) now defined in precision.h
// Note: Data structures (GridBounds, ExtendedGridBounds) now defined in types.h
// Note: Indexing macros (Y_SLICE, PENCIL, ZSLAB) now defined in types.h
// Note: Periodic boundary macros (PERIODIC_X, PERIODIC_Z) now defined in types.h

// ====================================================================================
// HELPER FUNCTION PROTOTYPES
// ====================================================================================
// Note: All function prototypes are now in module headers:
//   - FFT functions: fft/fft_setup.h
//   - Generation functions: generation/hermitian_generation.h
//   - Communication functions: communication/mpi_exchange.h
//   - Streaming functions: streaming/z_streaming.h
//   - Utility functions: utils/*.h
//
// Legacy streaming functions (x_streaming_unpack, unpack_and_FFT_x_row) are
// kept in main for backward compatibility but may be extracted in future phases.
// - Real MPI implementation (not OpenMP mock)
// - Each MPI rank generates multiple Y-slice pairs across batches
// - OpenMP used WITHIN each rank to parallelize X-Z loops during generation
// - 2D grid decomposition for redistribution to pencil layout
// - **NEW in v14**: PERIODIC BOUNDARY CONDITIONS FOR PADDING
//   * Padding uses periodic boundaries (wrap-around) instead of clamping
//   * Example: rank 0: [0,100] core → [-10,110] padded (wraps from right side!)
//             rank 1: [100,200] core → [90,210] padded (normal overlap)
//   * Physics-correct for cosmological simulations (periodic box)
//   * X=-10 wraps to X=N-10, X=N+5 wraps to X=5 (modulo arithmetic)
//   * Efficient indexing: PERIODIC_X macro with fast-path for common case
// - **RETAINED from v13**: OVERLAPPING X-REGIONS FOR ABACUS GRID READING
//   * Each rank stores extra X-values beyond its core region (X_PADDING on each side)
//   * Handles Abacus's imperfect grid reading (may access X outside strict bounds)
//   * Creates redundant storage at boundaries for seamless Abacus integration
//   * File output includes padded regions (some X-values written by multiple ranks)
// - **RETAINED from v12**: Z-SLAB STREAMING (ZELDOVICH-COMPATIBLE OUTPUT)
//   * Process one Z-slab at a time instead of X-rows
//   * Output format: [Z][Array][Y][X] matches Zeldovich's AZYX layout
//   * Each file contains all (X,Y) for a given Z-slab
//   * Direct compatibility with WriteParticlesSlab interface
// - **RETAINED from v11**: PERSISTENT RECV_BUFFER WITH SOURCE-GROUPED LAYOUT
//   * recv_buffer layout: [src_rank][y_idx_for_src][pencil_idx][array_idx]
//   * Direct receive into persistent buffer (NO per-batch recv_buffer_batch)
//   * Per-source write cursors track progress across batches
//   * Eliminates ~64 GB transient allocation per batch (single precision, N=32K)
// - Hermitian symmetry verification at multiple stages
// - Final result: Real-space data (imaginary parts ≈ 0)
// ====================================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>  // For INT_MAX
#include <mpi.h>
#include <omp.h>  // For hybrid MPI+OpenMP within each rank
#include <fftw3.h>
#include <sys/stat.h>
#include <errno.h>

// Include PCG RNG and STimer
#include "pcg-rng/pcg_random.hpp"
#ifdef __cplusplus
extern "C" {
#endif
#include "STimer.h"
#ifdef __cplusplus
}
#endif

// ====================================================================================
// CONFIGURATION AND TYPES
// ====================================================================================
// All compile-time flags and constants are now in config.h
// All type definitions and precision handling are in precision.h and types.h
#include "config.h"
#include "precision.h"
#include "types.h"

// ====================================================================================
// UTILITY MODULES
// ====================================================================================
// Utility functions extracted to separate modules
#include "utils/printing.h"
#include "utils/verification.h"
#include "utils/decomposition.h"
#include "utils/batch_helpers.h"

// ====================================================================================
// CORE MODULES
// ====================================================================================
// Core algorithm modules extracted in Phase 5
#include "fft/fft_setup.h"
#include "generation/hermitian_generation.h"
#include "communication/mpi_exchange.h"
#include "streaming/z_streaming.h"

// ====================================================================================
// PCG RNG WRAPPER FUNCTIONS
// ====================================================================================

// NOTE: Using matrix.c implementation (from fftw_openmp/matrix.c)
// The matrix.c version uses pcg64 with advance() for sequential initialization
// This matches the original approach used in other versions
//
// To switch back to the inline version below, comment out the matrix.c link
// in the Makefile and uncomment the code below.

// ========== INLINE VERSION (COMMENTED OUT) ==========
// static pcg32* global_pcg_array = NULL;
// static int global_pcg_L = 0, global_pcg_M = 0, global_pcg_N = 0;
// 
// void initialize_global_pcg(int L, int M, int N, uint64_t seed) {
//     global_pcg_L = L;
//     global_pcg_M = M;
//     global_pcg_N = N;
//     int total = L * M * N;
//     global_pcg_array = new pcg32[total];
//     for (int i = 0; i < total; i++) {
//         global_pcg_array[i].seed(seed, i);
//     }
// }
// 
// void cleanup_global_pcg() {
//     if (global_pcg_array) {
//         delete[] global_pcg_array;
//         global_pcg_array = NULL;
//     }
// }
// 
// double random_real_pcg_global(int slice_index) {
//     if (slice_index < 0 || slice_index >= global_pcg_L * global_pcg_M * global_pcg_N) {
//         return 0.0;
//     }
//     uint32_t val = global_pcg_array[slice_index]();
//     return (double)val / (double)UINT32_MAX;
// }

// ========== RNG MODULE ==========
// RNG functions are now in src/utils/rng.c (included via utils/rng.h above)

// Note: Precision types (real_t, fftw_complex_t, etc.) now defined in precision.h
// Note: Data structures (GridBounds, ExtendedGridBounds) now defined in types.h
// Note: Indexing macros (Y_SLICE, PENCIL, ZSLAB) now defined in types.h
// Note: Periodic boundary macros (PERIODIC_X, PERIODIC_Z) now defined in types.h

// ====================================================================================
// HELPER FUNCTION PROTOTYPES
// ====================================================================================
// Note: All function prototypes are now in module headers:
//   - FFT functions: fft/fft_setup.h
//   - Generation functions: generation/hermitian_generation.h
//   - Communication functions: communication/mpi_exchange.h
//   - Streaming functions: streaming/z_streaming.h
//   - Utility functions: utils/*.h
//
// Legacy streaming functions (x_streaming_unpack, unpack_and_FFT_x_row) are
// kept in main for backward compatibility but may be extracted in future phases.

// ====================================================================================
// HELPER FUNCTION IMPLEMENTATIONS
// ====================================================================================
// UTILITY FUNCTIONS
// ====================================================================================
// Note: The following utility functions have been extracted to separate modules:
//   - Printing functions: utils/printing.c
//   - Verification functions: utils/verification.c
//   - Decomposition functions: utils/decomposition.c
//   - Batch helper functions: utils/batch_helpers.c
// All function prototypes are available via the headers included above.
// ====================================================================================

// Removed: get_rank_batch_slice_count - now in utils/batch_helpers.c
// Removed: get_rank_batch_y_values - now in utils/batch_helpers.c
// Removed: calculate_batch_send_recv_counts - now in utils/batch_helpers.c
// Removed: print_y_slice_fourier - now in utils/printing.c
// Removed: verify_hermitian_pair - now in utils/verification.c
// Removed: get_grid_bounds - now in utils/decomposition.c
// Removed: get_extended_grid_bounds - now in utils/decomposition.c
// Removed: get_padded_bounds_simple - now in utils/decomposition.c
// Removed: verify_initial_fourier_hermitian_symmetry - now in utils/verification.c
// Removed: verify_self_conjugate_constraints - now in utils/verification.c
// Removed: count_missing_y_values - now in utils/verification.c (static helper)
// Removed: verify_pencil_completeness_with_flags - now in utils/verification.c
// Removed: verify_real_space_symmetry - now in utils/verification.c
// Removed: print_3d_matrix_visual - now in utils/printing.c

// All utility function implementations have been moved to utils/ modules

// ====================================================================================
// EXTRACTED FUNCTION IMPLEMENTATIONS
// ====================================================================================
// The following functions have been extracted to separate modules:
//   - setup_fftw_plans_full() → fft/fft_setup.c
//   - generate_hermitian_slice_pair_local() → generation/hermitian_generation.c
//   - exchange_metadata() → communication/mpi_exchange.c
//   - pack_slices_to_send_buffer() → communication/mpi_exchange.c
//   - unpack_recv_buffer_to_pencils() → communication/mpi_exchange.c
//   - z_streaming_unpack() → streaming/z_streaming.c
//
// All implementations are available via the module headers included above.
// ====================================================================================

// NOTE: generate_hermitian_slice_pair_local() is now implemented in 
// generation/hermitian_generation.c and declared in generation/hermitian_generation.h
// The old implementation has been removed to avoid duplicate definitions.

// All verification function implementations have been moved to utils/verification.c

// ----------------------------- unpack_and_FFT_x_row -----------------------------
// NEW in v10: Unpack ONE X-row from recv_buffer, apply 1D FFT, ready for writing
// This function processes a single X value across all Z values in the rank's region
// Enables streaming processing to reduce memory (one X-row at a time)
// ------------------------------------------------------------------------------------
void unpack_and_FFT_x_row(
    int rank, int num_ranks, int N, int narray,
    int x_global,                          // Which X-row to process (global coordinate)
    GridBounds my_bounds,                  // My (X,Z) region bounds
    int my_pencils,                        // Total pencils in my region (for indexing)
    fftw_complex_t *recv_buffer,            // Source data (from MPI_Ialltoallv)
    int *rdispls,                         // Offsets per source rank
    int **all_y_global_maps,              // Y-indices per source rank
    int *all_num_my_slices,               // Number of slices per source rank
    fftw_complex_t *local_pencils,          // Destination buffer (one X-row: z_count pencils)
    fftw_plan_t plan_1d_y)                  // FFT plan for 1D Y-direction
{
    /*
    Unpack ONE X-row from recv_buffer into local_pencils and apply 1D FFT.
    
    Key concept:
    - recv_buffer contains data organized by source rank: [src0][src1][src2]...
    - Each source's data: [array][slice][pencil] where pencil = (X,Z) in MY region
    - We extract data for x_global across all Z in my region
    - Apply 1D FFT along Y for each pencil in this X-row
    
    Memory savings: local_pencils holds only ONE X-row (z_count pencils),
    not all pencils (x_count × z_count pencils).
    */
    
    // Calculate dimensions
    int z_count = my_bounds.z_end - my_bounds.z_start;
    int local_x = x_global - my_bounds.x_start;
    
    // Sanity check
    if (x_global < my_bounds.x_start || x_global >= my_bounds.x_end) {
        fprintf(stderr, "[ERROR] Rank %d: x_global=%d out of bounds [%d,%d)\n",
                rank, x_global, my_bounds.x_start, my_bounds.x_end);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // ========== UNPACKING: Extract this X-row from recv_buffer ==========
    // Loop over all source ranks to collect Y-values for this X-row
    for (int src = 0; src < num_ranks; src++) {
        int src_offset = rdispls[src];
        int src_num_slices = all_num_my_slices[src];
        
        // Loop over arrays and slices from this source
        for (int array_idx = 0; array_idx < narray; array_idx++) {
            for (int slice_idx = 0; slice_idx < src_num_slices; slice_idx++) {
                int y_global = all_y_global_maps[src][slice_idx];
                
                // Loop over Z values in my region for this X
                #pragma omp parallel for
                for (int z = my_bounds.z_start; z < my_bounds.z_end; z++) {
                    int local_z = z - my_bounds.z_start;
                    
                    // Calculate pencil_idx in full region (for recv_buffer indexing)
                    // This matches the indexing used during packing
                    int pencil_idx = local_x * z_count + local_z;
                    
                    // Calculate position in recv_buffer
                    int64_t unpack_idx = (int64_t)array_idx * src_num_slices * my_pencils
                                       + slice_idx * my_pencils
                                       + pencil_idx;
                    
                    // Store in local_pencils using ONLY z_idx (since we only have one X-row)
                    // local_pencils layout: [Z][Array][Y] where Z is 0 to z_count-1
                    int z_idx = local_z;
                    PENCIL(z_idx, array_idx, y_global, N, narray)[0] = 
                        recv_buffer[src_offset + unpack_idx][0];
                    PENCIL(z_idx, array_idx, y_global, N, narray)[1] = 
                        recv_buffer[src_offset + unpack_idx][1];
                }
            }
        }
    }
    
    // ========== 1D FFT: Apply along Y-direction for each pencil ==========
    #pragma omp parallel for collapse(2)
    for (int z_idx = 0; z_idx < z_count; z_idx++) {
        for (int array_idx = 0; array_idx < narray; array_idx++) {
            // Get pointer to start of this pencil's Y-data
            fftw_complex_t *pencil_start = &PENCIL(z_idx, array_idx, 0, N, narray);
            FFTW_EXECUTE_DFT(plan_1d_y, pencil_start, pencil_start);
        }
    }
    
    // local_pencils now contains FFT'd data for this X-row, ready for writing
}

// ----------------------------- x_streaming_unpack -----------------------------
// V11: Unpack ONE X-row from recv_buffer (source-grouped) and apply 1D FFT
// recv_buffer layout: [src_rank][y_idx_for_src][pencil_idx][array_idx]
// local_pencils layout: [Z][Array][Y] (one X-row only)
// ---------------------------------------------------------------------------------
void x_streaming_unpack(
    int rank, int N, int narray,
    int x_global,                          // Which X-row to process
    GridBounds my_bounds,                  // My (X,Z) region bounds
    int my_pencils,                        // Total pencils in my region
    fftw_complex_t *recv_buffer,            // Source data (source-grouped)
    int64_t *recv_displs_src,              // Base offset per source
    int *src_total_slices,                 // Total slices per source
    int *y_owner_src,                      // Y → src mapping
    int *y_src_local_idx,                  // Y → local_idx mapping
    fftw_complex_t *local_pencils,          // Destination buffer (one X-row)
    fftw_plan_t plan_1d_y)                  // FFT plan for Y-direction
{
    (void)rank;  // Unused
    (void)src_total_slices;  // Unused but kept for API consistency
    
    int z_count = my_bounds.z_end - my_bounds.z_start;
    int local_x = x_global - my_bounds.x_start;
    
    // Sanity check
    if (x_global < my_bounds.x_start || x_global >= my_bounds.x_end) {
        fprintf(stderr, "[ERROR] Rank %d: x_global=%d out of bounds [%d,%d)\n",
                rank, x_global, my_bounds.x_start, my_bounds.x_end);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // ========== UNPACKING: Extract this X-row from recv_buffer ==========
    #pragma omp parallel for collapse(3)
    for (int z_idx = 0; z_idx < z_count; z_idx++) {
        for (int y = 0; y < N; y++) {
            for (int array_idx = 0; array_idx < narray; array_idx++) {
                // Find source and local index for this Y
                int src = y_owner_src[y];
                int y_local = y_src_local_idx[y];
                
                // How many total Y-slices from this source
                int src_num_slices = src_total_slices[src];
                
                // Calculate pencil index in full region
                int pencil_idx = local_x * z_count + z_idx;
                
                // Calculate position in recv_buffer
                // Layout from packing: [src][array][y_local][pencil]
                // This matches pack_slices_to_send_buffer() order
                int64_t recv_offset = recv_displs_src[src]
                                    + (int64_t)array_idx * src_num_slices * my_pencils
                                    + (int64_t)y_local * my_pencils
                                    + pencil_idx;
                
                // Copy to local_pencils: [Z][Array][Y]
                PENCIL(z_idx, array_idx, y, N, narray)[0] = recv_buffer[recv_offset][0];
                PENCIL(z_idx, array_idx, y, N, narray)[1] = recv_buffer[recv_offset][1];
            }
        }
    }
    
    // ========== 1D FFT: Apply along Y-direction for each pencil ==========
    #pragma omp parallel for collapse(2)
    for (int z_idx = 0; z_idx < z_count; z_idx++) {
        for (int array_idx = 0; array_idx < narray; array_idx++) {
            fftw_complex_t *pencil_start = &PENCIL(z_idx, array_idx, 0, N, narray);
            FFTW_EXECUTE_DFT(plan_1d_y, pencil_start, pencil_start);
        }
    }
    
    // local_pencils now contains FFT'd data for this X-row, ready for writing
}

// ====================================================================================
// MAIN FUNCTION
// ====================================================================================

int main(int argc, char **argv)
{
    // ========================================================================
    // MPI INITIALIZATION
    // ========================================================================
    
    MPI_Init(&argc, &argv);
    
    int rank, num_ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    
    // ========================================================================
    // PARSE ARGUMENTS
    // ========================================================================
    
    int N = 64;
    const char* param_file = NULL;
    
    if (argc < 2) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s N [param_file.par]\n", argv[0]);
            fprintf(stderr, "  N: Grid size (e.g., 256)\n");
            fprintf(stderr, "  param_file.par: Optional zeldovich-PLT parameter file\n");
            fprintf(stderr, "                 If provided, enables power spectrum mode\n");
        }
        MPI_Finalize();
        return 1;
    }
    
    N = atoi(argv[1]);
    if (N <= 0 || N % 2 != 0) { 
        if (rank == 0) {
            fprintf(stderr, "Error: N must be positive even integer\n");
        }
        MPI_Finalize();
        return 1; 
    }
    
    // Optional parameter file (enables power spectrum mode)
    if (argc >= 3) {
        param_file = argv[2];
    }
    
    // ========================================================================
    // VERIFY RANK COUNT
    // ========================================================================
    
    int total_pairs = N/2 + 1;  // Number of Y-slice pairs to distribute
    
    // Allow any number of ranks (within reason)
    if (num_ranks < 1 || num_ranks > total_pairs) {
        if (rank == 0) {
            fprintf(stderr, "Error: Invalid num_ranks=%d for N=%d (total_pairs=%d)\n",
                    num_ranks, N, total_pairs);
        }
        MPI_Finalize();
        return 1;
    }
    
    // Calculate pair distribution for informational purposes
    int pairs_per_rank_base = total_pairs / num_ranks;
    int remainder_pairs = total_pairs % num_ranks;
    
    // Warn about potential load imbalance
    if (num_ranks < total_pairs && rank == 0) {
        printf("[INFO] num_ranks=%d < total_pairs=%d: Each rank will process %d-%d pairs\n",
               num_ranks, total_pairs, pairs_per_rank_base, 
               pairs_per_rank_base + (remainder_pairs > 0 ? 1 : 0));
    }
    
    // NOTE: In multi-batch mode, all ranks participate in the computation
    // Each rank processes multiple Y-slice pairs as needed
    
    // ========================================================================
    // PRINT CONFIGURATION
    // ========================================================================
    
    if (rank == 0) {
        printf("====================================================================================\n");
        printf("VERSION 14: PERIODIC BOUNDARY CONDITIONS FOR PADDING\n");
        printf("====================================================================================\n");
        printf("Precision: %s (%d bytes per complex number)\n", PRECISION_NAME, BYTES_PER_COMPLEX);
        printf("Matrix size: N = %d³, Total elements: %zu\n", N, (size_t)N * (size_t)N * (size_t)N);
        printf("MPI ranks: %d (total_pairs: %d for N=%d)\n", num_ranks, total_pairs, N);
        printf("Multi-batch processing: Each rank processes multiple Y-slice pairs\n");
#if USE_X_PADDING
        printf("V14 feature: Periodic boundary conditions (X_PADDING=%d per side)\n", X_PADDING);
        printf("            Physics-correct wrap-around at boundaries (not clamped!)\n");
        printf("            Example: rank 0 core [0,100] → padded [-10,110] (wraps from right!)\n");
        printf("                     rank 1 core [100,200] → padded [90,210] (normal overlap)\n");
        printf("                     X=-10 wraps to X=N-10, X=N+5 wraps to X=5\n");
        printf("            Efficient PERIODIC_X macro: fast-path for common case x∈[0,N)\n");
#else
        printf("Grid decomposition: Core grid only (no padding, X_PADDING=0)\n");
        printf("            Each rank owns non-overlapping X-region\n");
        printf("            No periodic boundary wrapping - standard decomposition\n");
#endif
        printf("V12 feature: Z-slab streaming (Zeldovich-compatible output)\n");
        printf("            Output format: [Z][Array][Y][X] matches Zeldovich AZYX\n");
        printf("            Each file contains all (X,Y) for one Z-slab\n");
        printf("V11 feature: Direct receive into persistent source-grouped buffer\n");
        printf("Communication method: %s\n",
               OVERLAP_COMMUNICATION == 0 ? "MPI_Ialltoallv (Phase A)" : "Isend/Irecv with overlap (Phase B)");
        printf("====================================================================================\n\n");
    }
    
    // ========================================================================
    // MULTI-BATCH SLICE DISTRIBUTION
    // ========================================================================
    
    // Structure to hold a pair of Y-slices (primary + conjugate)
    typedef struct {
        int y_primary;    // Primary Y-index
        int y_mirror;     // Conjugate Y-index (may equal y_primary if self-conjugate)
        bool is_self_conjugate;  // True if y_primary == y_mirror
    } YSlicePair;
    
    // STEP 1: Calculate total number of pairs to distribute
    // Pairs: (0,0), (1,N-1), (2,N-2), ..., up to (N/2, N/2)
    // Note: total_pairs already calculated above in rank verification
    
    // STEP 2: Distribute pairs across ranks
    int pairs_per_rank = total_pairs / num_ranks;
    int remainder = total_pairs % num_ranks;
    
    // Each rank gets pairs_per_rank, and first 'remainder' ranks get +1
    int my_num_pairs;
    int my_pair_start;
    
    if (rank < remainder) {
        my_num_pairs = pairs_per_rank + 1;
        my_pair_start = rank * (pairs_per_rank + 1);
    } else {
        my_num_pairs = pairs_per_rank;
        my_pair_start = remainder * (pairs_per_rank + 1) + (rank - remainder) * pairs_per_rank;
    }
    
    // STEP 3: Create list of pairs assigned to this rank
    YSlicePair *my_pair_list = NULL;
    bool is_idle_rank = (my_num_pairs == 0);
    
    if (!is_idle_rank) {
        my_pair_list = (YSlicePair*)malloc(sizeof(YSlicePair) * my_num_pairs);
        
        for (int i = 0; i < my_num_pairs; i++) {
            int pair_idx = my_pair_start + i;
            int y_primary = pair_idx;
            int y_mirror = (y_primary == 0 || y_primary == N/2) ? y_primary : N - y_primary;
            
            my_pair_list[i].y_primary = y_primary;
            my_pair_list[i].y_mirror = y_mirror;
            my_pair_list[i].is_self_conjugate = (y_primary == y_mirror);
        }
        
        // Debug output for first few ranks
        if (DEBUG_PRINTS && rank < 5) {
            printf("[INIT] Rank %d: Assigned %d pairs:\n", rank, my_num_pairs);
            for (int i = 0; i < my_num_pairs && i < 3; i++) {  // Show first 3 pairs
                printf("  Pair %d: Y=%d", i, my_pair_list[i].y_primary);
                if (!my_pair_list[i].is_self_conjugate) {
                    printf(" and Y=%d (conjugate)", my_pair_list[i].y_mirror);
                } else {
                    printf(" (self-conjugate)");
                }
                printf("\n");
            }
            if (my_num_pairs > 3) {
                printf("  ... (%d more pairs)\n", my_num_pairs - 3);
            }
        }
    } else {
        // IDLE RANK: No pairs assigned
        if (rank == num_ranks - 1 || (rank < num_ranks && rank + 1 >= num_ranks)) {
            printf("[WARNING] Rank %d is IDLE (total_pairs=%d, num_ranks=%d)\n",
                   rank, total_pairs, num_ranks);
        }
    }
    
    // STEP 4: Set up data structures (compatible with existing code)
    // Note: These will be used during batch processing
    int num_my_slices = 0;  // Will be set per-batch (1 or 2)
    // V11: No y_global_map needed (passing NULL to pack_slices_to_send_buffer)
    fftw_complex_t *local_y_slices = NULL;  // Allocated once, reused per-batch
    fftw_complex_t *local_z_slab = NULL;    // V12: Allocated once for Z-slab streaming
    ExtendedGridBounds my_extended_bounds;  // V13: Extended bounds with padding (if enabled)
    int my_pencils = 0;                     // Will be set to padded or core pencil count
    
    // ========================================================================
    // V11: NEW DATA STRUCTURES FOR PERSISTENT RECV_BUFFER
    // ========================================================================
    int *src_total_slices = NULL;       // [num_ranks] - total Y-slices from each source
    int64_t *recv_displs_src = NULL;    // [num_ranks] - base offset per source (elements)
    int *src_write_cursor = NULL;       // [num_ranks] - current write position per source
    int *y_owner_src = NULL;            // [N] - which rank generated Y=i
    int *y_src_local_idx = NULL;        // [N] - index within that source's chunk
    int *y_batch_idx = NULL;            // [N] - which batch Y=i came from (for unpacking)
    int *y_slice_idx_in_batch = NULL;     // [N] - slice_idx (0 or 1) in that batch (for unpacking)
    int **src_batch_slice_counts = NULL; // [num_ranks][global_max_batches] - slices per batch per source
    fftw_complex_t *recv_buffer = NULL; // Persistent recv buffer (source-grouped layout)
    int64_t recv_total_elems = 0;       // Total elements in recv_buffer
    
    if (is_idle_rank) {
        // Set up dummy structures for idle rank
        num_my_slices = 0;
        // V11: No y_global_map needed
        
        // V13: Initialize extended bounds to zero
        my_extended_bounds.core.x_start = my_extended_bounds.core.x_end = 0;
        my_extended_bounds.core.z_start = my_extended_bounds.core.z_end = 0;
        my_extended_bounds.padded.x_start = my_extended_bounds.padded.x_end = 0;
        my_extended_bounds.padded.z_start = my_extended_bounds.padded.z_end = 0;
        my_extended_bounds.num_pencils_core = 0;
        my_extended_bounds.num_pencils_padded = 0;
        my_pencils = 0;
    }
    
    // STEP 5: Print distribution summary
    if (rank == 0) {
        printf("\n[DISTRIBUTION] Multi-batch slice distribution:\n");
        printf("  Total Y-slice pairs: %d\n", total_pairs);
        printf("  Active ranks: %d\n", num_ranks);
        printf("  Pairs per rank: %d base + %d ranks with +1\n", 
               pairs_per_rank, remainder);
        printf("  Rank 0 processes: %d pairs (Y indices %d-%d)\n",
               my_num_pairs, my_pair_start, my_pair_start + my_num_pairs - 1);
        int total_slices_estimate = total_pairs * 2 - 2;  // Account for 2 self-conjugate
        printf("  Estimated total Y-slices: ~%d\n\n", total_slices_estimate);
    }
    
    // ========================================================================
    // INITIALIZE RNG (All ranks have identical RNG arrays)
    // ========================================================================
    
    uint64_t seed = 4;
    initialize_global_pcg(N, N, N, seed);
    // NOTE: PCG array = (N/2) × 64 bytes = O(N) memory per rank
    //       N=1024: 32 KB,  N=32,000: 1 MB  ← Manageable!
    // Each rank allocates identical array (same seed → reproducible RNG)
    // All OpenMP threads within rank share this array (thread-safe: different indices)
    // Can free after Stage 1 (after slice generation)
    
    // ========================================================================
    // ZELDOVICH-PLT PARAMETERS AND POWER SPECTRUM (v15.2)
    // ========================================================================
    // Load parameters from file if provided, create PowerSpectrum object
    // Each rank creates its own objects (same file → identical objects)
    
    ParametersHandle params = NULL;
    PowerSpectrumHandle ps = NULL;
    
    // Legacy power spectrum parameters (for backward compatibility)
    power_spectrum_params_t ps_params;
    power_spectrum_params_t *ps_params_ptr = NULL;  // NULL = use uniform RNG (backward compatible)
    
    // TODO: Read these from config file 
    // For now, use simple power law defaults (can be enabled by uncommenting):
    init_power_spectrum_params(&ps_params, 
        -2.0,      // powerlaw_index: P(k) ~ k^-2 (scale-invariant)
        1.0,       // normalization: adjust based on desired amplitude
        0.0,       // Pk_smooth: no smoothing
        0          // fixed_power: random amplitude (not fixed)
    );
    ps_params_ptr = &ps_params;  // Enable power spectrum mode
    
    if (param_file != NULL) {
        // Load parameters from file
        if (rank == 0) {
            printf("[INIT] Loading zeldovich-PLT parameters from: %s\n", param_file);
        }
        
        params = zeldovich_params_create(param_file);
        if (!params) {
            if (rank == 0) {
                fprintf(stderr, "ERROR: Failed to load parameter file: %s\n", param_file);
                fprintf(stderr, "       Check that file exists and has correct format\n");
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        // Verify N matches ppd from parameter file
        int64_t ppd = zeldovich_params_get_ppd(params);
        if (ppd != N) {
            if (rank == 0) {
                fprintf(stderr, "WARNING: Parameter file ppd=%ld != N=%d\n", ppd, N);
                fprintf(stderr, "         Using N=%d from command line\n", N);
            }
        }
        
        // Create PowerSpectrum object (each rank creates its own)
        // Spline resolution: typically 128 for good accuracy
        ps = zeldovich_ps_create(128, params);
        if (!ps) {
            if (rank == 0) {
                fprintf(stderr, "ERROR: Failed to create PowerSpectrum object\n");
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        // Initialize power spectrum (power law mode)
        // Get power law index from parameters
        double powerlaw_index = zeldovich_params_get_Pk_powerlaw_index(params);
        
        // Check if power law index is valid (not the default invalid value)
        if (powerlaw_index == 1000.0) {
            if (rank == 0) {
                fprintf(stderr, "ERROR: ZD_Pk_powerlaw_index not specified in parameter file\n");
                fprintf(stderr, "       Add 'ZD_Pk_powerlaw_index = -2.0' (or desired value) to parameter file\n");
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (zeldovich_ps_init_powerlaw(ps, powerlaw_index, params) != 0) {
            if (rank == 0) {
                fprintf(stderr, "ERROR: Failed to initialize power spectrum (power law index: %.2f)\n", powerlaw_index);
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (rank == 0) {
            double boxsize = zeldovich_params_get_boxsize(params);
            double fundamental = zeldovich_params_get_fundamental(params);
            int seed = zeldovich_params_get_seed(params);
            printf("[INIT] PowerSpectrum initialized:\n");
            printf("       Power law index: %.2f\n", powerlaw_index);
            printf("       Box size: %.2f\n", boxsize);
            printf("       Fundamental wavenumber: %.6e\n", fundamental);
            printf("       Seed: %d\n", seed);
            printf("       Power spectrum mode: ENABLED\n");
        }
    } else {
        if (rank == 0) {
            printf("[INIT] No parameter file provided, using uniform random mode\n");
        }
    }

    // =========================================================================
    // TEMPORARY DEBUG OPTION: Disable zeldovich-PLT power spectrum (ps/params)
    // =========================================================================
    //
    // For heap corruption / segmentation fault isolation, it is useful to
    // confirm whether the issue originates in the zeldovich-PLT integration
    // (PowerSpectrum / Parameters / v2rng) or in the MPI / packing code.
    //
    // When this block is enabled, the code runs in legacy fallback mode:
    //   - ps_handle     = NULL
    //   - params_handle = NULL
    //   - ps_params_ptr = NULL (uniform RNG mode)
    //
    // This allows testing whether crashes disappear when zeldovich-PLT is
    // completely bypassed.
    //
    // To re-enable power spectrum mode, change the #if from 1 to 0.
    //
    #if 0  // DEBUG: Force-disable power spectrum / zeldovich-PLT
    // Set to 1 to disable power spectrum for debugging
    if (rank == 0) {
        printf("[DEBUG] Power spectrum mode DISABLED for debugging "
               "(ps_handle=NULL, params_handle=NULL, ps_params_ptr=NULL)\n");
    }
    ps_params_ptr = NULL;
    ps = NULL;
    params = NULL;
    #endif
    
    // ========================================================================
    // SETUP NARRAY (number of arrays per slice)
    // ========================================================================
    
    int narray = NARRAY;  // Use the constant defined above (available throughout main)
    
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[SETUP] narray = %d (arrays per Y-slice)\n", narray);
    }
    
    // ========================================================================
    // DECLARE FFT PLANS (available throughout main)
    // ========================================================================
    fftw_plan_t plan_2d, plan_1d_y;
    
    // ========================================================================
    // STAGE 3: SETUP FFT PLANS (BEFORE DATA ALLOCATION)
    // ========================================================================
    // CRITICAL: Create FFT plans using dummy memory before allocating actual data
    // This prevents data destruction during planning (FFTW_MEASURE/PATIENT modes)
    
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[SETUP] Creating FFT plans with dummy memory...\n");
    }
    
    // Create plans using dummy memory (before allocating actual data)
    setup_fftw_plans_full(N, narray, &plan_2d, &plan_1d_y);
    
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[SETUP] FFT plans created successfully (2D and 1D)\n");
    }
    
    // ========================================================================
    // ALLOCATE MY Y-SLICES (active ranks only)
    // ========================================================================
    
    if (!is_idle_rank) {
        // MEMORY: num_my_slices × narray × N² × 16 bytes (e.g., 131.2 GB for N=32K, 2 slices, narray=4)
        // Can free after Stage 3 packing complete (after send_buffer filled)
        // NEW: Single flat allocation for all arrays
        int64_t total_size = (int64_t)num_my_slices * narray * N * N;
        
        if (posix_memalign((void**)&local_y_slices, ALIGN_BYTES, 
                           sizeof(fftw_complex_t) * total_size) != 0) {
            fprintf(stderr, "Rank %d: posix_memalign failed for local_y_slices\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        memset(local_y_slices, 0, sizeof(fftw_complex_t) * total_size);
        
        if (rank == 0 && DEBUG_PRINTS) {
            size_t total_bytes = (size_t)total_size * sizeof(fftw_complex_t);
            printf("[MEMORY] Allocated local_y_slices: %zu bytes (%.2f GB) for %d slices, %d arrays\n",
                   total_bytes, total_bytes / (1024.0 * 1024.0 * 1024.0), num_my_slices, narray);
        }
        
        // ========================================================================
        // STEP 1d: VERIFY MACRO INDEXING
        // ========================================================================
        #if DEBUG_PRINTS
        // Test macro indexing (add after allocation, before generation)
        if (rank < 3) {  // Only test on first few ranks to reduce output
            // Fill test pattern
            for (int s = 0; s < num_my_slices; s++) {
                for (int a = 0; a < narray; a++) {
                    Y_SLICE(s, a, 5, 3, N, narray)[0] = s * 100 + a * 10;
                    Y_SLICE(s, a, 5, 3, N, narray)[1] = s * 100 + a * 10 + 1;
                }
            }
            
            // Verify
            bool indexing_ok = true;
            for (int s = 0; s < num_my_slices; s++) {
                for (int a = 0; a < narray; a++) {
                    double expected_re = s * 100 + a * 10;
                    double expected_im = s * 100 + a * 10 + 1;
                    double actual_re = Y_SLICE(s, a, 5, 3, N, narray)[0];
                    double actual_im = Y_SLICE(s, a, 5, 3, N, narray)[1];
                    
                    if (fabs_t(expected_re - actual_re) > 1e-10 || 
                        fabs_t(expected_im - actual_im) > 1e-10) {
                        printf("[ERROR] Rank %d: Macro indexing failed at (slice=%d, array=%d)\n", 
                               rank, s, a);
                        printf("        Expected: %.2f + %.2fi, Got: %.2f + %.2fi\n",
                               expected_re, expected_im, actual_re, actual_im);
                        indexing_ok = false;
                    }
                }
            }
            if (indexing_ok) {
                printf("[OK] Rank %d: Macro indexing verified (slices=%d, arrays=%d)\n", 
                       rank, num_my_slices, narray);
            } else {
                fprintf(stderr, "[ERROR] Rank %d: Macro indexing test FAILED!\n", rank);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            // Clear test pattern
            memset(local_y_slices, 0, sizeof(fftw_complex_t) * total_size);
        }
        #endif
    }
    
    // ========================================================================
    // STAGE 1: GENERATE MY PAIR (with 2D FFT)
    // ========================================================================
    
    // ========================================================================
    // NOTE: OLD STAGE 1 (SINGLE-PAIR GENERATION) REMOVED
    // ========================================================================
    // Generation and communication are now integrated in the MULTI-BATCH LOOP below
    // (after metadata exchange and before streaming processing)
    
    // V11: No y_global_map needed - using y_owner_src + y_src_local_idx instead
    
    // ========================================================================
    // STAGE 2: METADATA EXCHANGE
    // ========================================================================
    
    // NOTE: Stage 2 (Metadata exchange) is INTEGRATED into multi-batch processing
    // Each rank already knows its own Y-values from y_global_map (built from my_pair_list)
    // No need for separate exchange_metadata() call
    
    // V14: Verify grid decomposition (show both core and padded with periodic BC)
    // Calculate grid factors first (needed for get_extended_grid_bounds)
    int grid_x_verify, grid_z_verify;
    calculate_grid_factors(num_ranks, &grid_x_verify, &grid_z_verify);
    
    // Validate Abacus compatibility (if enabled)
#if ENABLE_ABACUS_VALIDATION
    if (rank == 0) {
        if (!validate_abacus_compatibility(N, num_ranks, grid_x_verify, grid_z_verify)) {
            fprintf(stderr, "[WARNING] Domain decomposition may not be Abacus-compatible\n");
            fprintf(stderr, "          Abacus requires: N divisible by grid_x and grid_z (exact division)\n");
        }
    }
#endif
    
    if (rank == 0 && DEBUG_PRINTS) {
        printf("\n[V14-DEBUG] Grid decomposition verification (PERIODIC BOUNDARIES):\n");
        int num_to_print = (num_ranks < 4) ? num_ranks : 4;
        for (int dest = 0; dest < num_to_print; dest++) {
            ExtendedGridBounds ext_b = get_extended_grid_bounds(dest, N, num_ranks, grid_x_verify, grid_z_verify);
            printf("  Rank %d:\n", dest);
            printf("    Core:   X=[%d,%d), Z=[%d,%d), Pencils=%d\n",
                   ext_b.core.x_start, ext_b.core.x_end, 
                   ext_b.core.z_start, ext_b.core.z_end, ext_b.num_pencils_core);
            printf("    Padded: X=[%d,%d), Z=[%d,%d), Pencils=%d",
                   ext_b.padded.x_start, ext_b.padded.x_end, 
                   ext_b.padded.z_start, ext_b.padded.z_end, ext_b.num_pencils_padded);
            
            // Show if this rank has periodic wrapping
            if (ext_b.padded.x_start < 0 || ext_b.padded.x_end > N) {
                printf(" [PERIODIC WRAP]");
            }
            printf("\n");
        }
        if (num_ranks > 4) {
            printf("  ... (showing first 4 ranks only)\n");
        }
        printf("  Note: Negative X-start or X-end > N indicates periodic wrap-around\n");
    }
    
    if (rank == 0) {
        printf("[Stage 2] Complete.\n");
    }
    
    // ========================================================================
    // MULTI-BATCH PROCESSING LOOP (NEW IN MULTI-BATCH V10)
    // ========================================================================
    // Process each Y-slice pair assigned to this rank in batches
    // Each batch: Generate → 2D FFT → Pack → Communicate → Accumulate
    // After ALL batches: Streaming unpack → 1D FFT → Write
    
    if (rank == 0) {
        printf("\n[MULTI-BATCH] Starting batch processing...\n");
        printf("              Ranks will process %d pairs each (approx)\n", 
               is_idle_rank ? 0 : my_num_pairs);
    }
    
    // ========================================================================
    // STEP 0: CALCULATE PENCIL REGION (WITH PADDING IF ENABLED, NEEDED FOR BUFFER SIZES)
    // ========================================================================
    
    // V13: Need to calculate grid factors first to call get_extended_grid_bounds
    int grid_x, grid_z;
    calculate_grid_factors(num_ranks, &grid_x, &grid_z);
    
    // Validate Abacus compatibility (if enabled)
#if ENABLE_ABACUS_VALIDATION
    if (rank == 0 && !validate_abacus_compatibility(N, num_ranks, grid_x, grid_z)) {
        fprintf(stderr, "[WARNING] Domain decomposition may not be Abacus-compatible\n");
        fprintf(stderr, "          Abacus requires: N divisible by grid_x and grid_z (exact division)\n");
    }
#endif
    
    if (!is_idle_rank) {
        // V13: Get extended bounds with X-padding for Abacus compatibility
        my_extended_bounds = get_extended_grid_bounds(rank, N, num_ranks, grid_x, grid_z);
        
        // Use appropriate pencil count (padded if enabled, core otherwise)
#if USE_X_PADDING
        my_pencils = my_extended_bounds.num_pencils_padded;
#else
        my_pencils = my_extended_bounds.num_pencils_core;
#endif
        
        if (DEBUG_PRINTS && rank < 2) {
#if USE_X_PADDING
            printf("[V14-PENCILS] Rank %d (PERIODIC BOUNDARIES):\n", rank);
            printf("  Core region:   X=[%d,%d), Z=[%d,%d), pencils=%d\n",
                   my_extended_bounds.core.x_start, my_extended_bounds.core.x_end,
                   my_extended_bounds.core.z_start, my_extended_bounds.core.z_end,
                   my_extended_bounds.num_pencils_core);
            printf("  Padded region: X=[%d,%d), Z=[%d,%d), pencils=%d (STORAGE SIZE)\n",
                   my_extended_bounds.padded.x_start, my_extended_bounds.padded.x_end,
                   my_extended_bounds.padded.z_start, my_extended_bounds.padded.z_end,
                   my_extended_bounds.num_pencils_padded);
            
            // Show periodic wrapping info
            if (my_extended_bounds.padded.x_start < 0) {
                printf("  ⚠ Left wrap: X=[%d,0) wraps to X=[%d,%d) (periodic BC)\n",
                       my_extended_bounds.padded.x_start, 
                       N + my_extended_bounds.padded.x_start, N);
            }
            if (my_extended_bounds.padded.x_end > N) {
                printf("  ⚠ Right wrap: X=[%d,%d) wraps to X=[0,%d) (periodic BC)\n",
                       N, my_extended_bounds.padded.x_end, 
                       my_extended_bounds.padded.x_end - N);
            }
            
            int x_overlap_left = my_extended_bounds.core.x_start - my_extended_bounds.padded.x_start;
            int x_overlap_right = my_extended_bounds.padded.x_end - my_extended_bounds.core.x_end;
            printf("  X-padding: %d values on left, %d values on right\n", 
                   x_overlap_left, x_overlap_right);
#else
            printf("[GRID] Rank %d (NO PADDING - CORE GRID ONLY):\n", rank);
            printf("  Core region:   X=[%d,%d), Z=[%d,%d), pencils=%d\n",
                   my_extended_bounds.core.x_start, my_extended_bounds.core.x_end,
                   my_extended_bounds.core.z_start, my_extended_bounds.core.z_end,
                   my_extended_bounds.num_pencils_core);
#endif
        }
    }
    
    // V11: No recv_buffer_accumulator, y_value_received, or y_to_local_idx needed
    // Receiving directly into persistent recv_buffer (source-grouped layout)
    
    // ========================================================================
    // STEP 1: CALCULATE MAXIMUM BATCH COUNT
    // ========================================================================
    // All ranks must loop over the same number of batches (some may have no data in later batches)
    
    int max_batches = 0;
    if (!is_idle_rank && my_num_pairs > max_batches) {
        max_batches = my_num_pairs;
    }
    
    // Share maximum across all ranks
    int global_max_batches = 0;
    MPI_Allreduce(&max_batches, &global_max_batches, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[BATCHES] Total batches to process: %d\n", global_max_batches);
    }
    
    // ========================================================================
    // V11: PHASE 1 - PRE-COMPUTATION (BEFORE BATCH LOOP)
    // ========================================================================
    // Calculate per-source totals and allocate persistent recv_buffer
    
    // ========================================================================
    // V11: STEP 1 - COMPUTE PER-SOURCE TOTALS ACROSS ALL BATCHES
    // ========================================================================
    
    if (rank == 0 && DEBUG_PRINTS) {
        printf("\n[V11-SETUP] Computing per-source Y-slice totals...\n");
    }
    
    src_total_slices = (int*)calloc(num_ranks, sizeof(int));
    
    // Count how many Y-slices each source will send across ALL batches
    for (int src = 0; src < num_ranks; src++) {
        for (int batch = 0; batch < global_max_batches; batch++) {
            int count = get_rank_batch_slice_count(src, batch, N, num_ranks);
            src_total_slices[src] += count;
        }
    }
    
    if (DEBUG_PRINTS && rank < 2) {
        printf("[V11-SETUP] Rank %d: Source totals: ", rank);
        for (int src = 0; src < (num_ranks < 4 ? num_ranks : 4); src++) {
            printf("src%d=%d ", src, src_total_slices[src]);
        }
        if (num_ranks > 4) printf("...");
        printf("\n");
    }
    
    // ========================================================================
    // V11: STEP 2 - COMPUTE RECEIVE DISPLACEMENTS (PREFIX SUM)
    // ========================================================================
    
    recv_displs_src = (int64_t*)malloc(sizeof(int64_t) * num_ranks);
    recv_total_elems = 0;
    
    for (int src = 0; src < num_ranks; src++) {
        recv_displs_src[src] = recv_total_elems;
        recv_total_elems += (int64_t)src_total_slices[src] * my_pencils * narray;
    }
    
    if (rank == 0 && DEBUG_PRINTS) {
        size_t recv_bytes = (size_t)recv_total_elems * sizeof(fftw_complex_t);
        printf("[V11-SETUP] Persistent recv_buffer size: %zu bytes (%.2f GB)\n",
               recv_bytes, recv_bytes / (1024.0 * 1024.0 * 1024.0));
        printf("            Total elements: %lld = sum of (src_total_slices[src] * my_pencils * narray)\n",
               (long long)recv_total_elems);
    }
    
    // ========================================================================
    // V11: STEP 3 - ALLOCATE PERSISTENT RECV_BUFFER (SOURCE-GROUPED LAYOUT)
    // ========================================================================
    // Layout: [src_rank][y_idx_for_src][pencil_idx][array_idx]
    // This replaces V10's recv_buffer_accumulator
    
    if (!is_idle_rank && recv_total_elems > 0) {
        if (posix_memalign((void**)&recv_buffer, ALIGN_BYTES,
                           sizeof(fftw_complex_t) * recv_total_elems) != 0) {
            fprintf(stderr, "Rank %d: posix_memalign failed for recv_buffer\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        memset(recv_buffer, 0, sizeof(fftw_complex_t) * recv_total_elems);
        
        if (rank == 0 && DEBUG_PRINTS) {
            printf("[V11-MEMORY] Allocated recv_buffer: %.2f GB (source-grouped layout)\n",
                   (recv_total_elems * sizeof(fftw_complex_t)) / (1024.0 * 1024.0 * 1024.0));
        }
    }
    
    // ========================================================================
    // V11: STEP 4 - INITIALIZE WRITE CURSORS
    // ========================================================================
    
    src_write_cursor = (int*)calloc(num_ranks, sizeof(int));
    
    if (DEBUG_PRINTS && rank < 2) {
        printf("[V11-SETUP] Rank %d: Write cursors initialized to 0\n", rank);
    }
    
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[V11-SETUP] Phase 1 pre-computation complete.\n\n");
    }
    
    // ========================================================================
    // V11: PHASE 2 - BUILD Y → (SRC, LOCAL_IDX) MAPPING
    // ========================================================================
    
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[V11-SETUP] Building Y-value ownership mapping...\n");
    }
    
    y_owner_src = (int*)malloc(sizeof(int) * N);
    y_src_local_idx = (int*)malloc(sizeof(int) * N);
    y_batch_idx = (int*)malloc(sizeof(int) * N);
    y_slice_idx_in_batch = (int*)malloc(sizeof(int) * N);
    
    // Initialize to invalid
    for (int y = 0; y < N; y++) {
        y_owner_src[y] = -1;
        y_src_local_idx[y] = -1;
        y_batch_idx[y] = -1;
        y_slice_idx_in_batch[y] = -1;
    }
    
    // Allocate src_batch_slice_counts to track slices per batch per source
    src_batch_slice_counts = (int**)malloc(sizeof(int*) * num_ranks);
    for (int src = 0; src < num_ranks; src++) {
        src_batch_slice_counts[src] = (int*)calloc(global_max_batches, sizeof(int));
    }
    
    // Build mapping by simulating what each source generates
    int *src_y_counter = (int*)calloc(num_ranks, sizeof(int));
    
    for (int src = 0; src < num_ranks; src++) {
        for (int batch = 0; batch < global_max_batches; batch++) {
            int y_values[2];
            int count;
            get_rank_batch_y_values(src, batch, N, num_ranks, y_values, &count);
            
            // Store slice count for this batch from this source
            src_batch_slice_counts[src][batch] = count;
            
            for (int i = 0; i < count; i++) {
                int y_global = y_values[i];
                
                // Ensure not already assigned (debug check)
                if (y_owner_src[y_global] != -1) {
                    fprintf(stderr, "[ERROR] Y=%d already assigned to src %d, now src %d wants it!\n",
                           y_global, y_owner_src[y_global], src);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                
                y_owner_src[y_global] = src;
                y_src_local_idx[y_global] = src_y_counter[src];
                y_batch_idx[y_global] = batch;
                y_slice_idx_in_batch[y_global] = i;  // i is the slice_idx (0 or 1) in this batch
                src_y_counter[src]++;
            }
        }
    }
    
    free(src_y_counter);
    
    // Verify all Y-values mapped
    int unmapped_count = 0;
    for (int y = 0; y < N; y++) {
        if (y_owner_src[y] == -1) {
            if (unmapped_count < 5) {
                fprintf(stderr, "[ERROR] Rank %d: Y=%d not mapped to any source!\n", rank, y);
            }
            unmapped_count++;
        }
    }
    
    if (unmapped_count > 0) {
        fprintf(stderr, "[ERROR] Rank %d: %d Y-values unmapped!\n", rank, unmapped_count);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (DEBUG_PRINTS && rank < 2) {
        printf("[V11-SETUP] Rank %d: Y-mapping complete. Sample: Y=0→src%d[%d], Y=%d→src%d[%d]\n",
               rank, y_owner_src[0], y_src_local_idx[0], N/2, y_owner_src[N/2], y_src_local_idx[N/2]);
    }
    
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[V11-SETUP] Phase 2 Y-mapping complete.\n\n");
    }
    
    // V11: No y_global_map needed - using y_owner_src + y_src_local_idx instead
    
    // ========================================================================
    // STEP 4: ALLOCATE local_y_slices (REUSED PER BATCH)
    // ========================================================================
    
    if (!is_idle_rank) {
        // Allocate buffer for maximum 2 slices (conjugate pair)
        int max_slices_per_batch = 2;
        int64_t slice_buffer_size = (int64_t)max_slices_per_batch * narray * N * N;
        
        if (posix_memalign((void**)&local_y_slices, ALIGN_BYTES,
                           sizeof(fftw_complex_t) * slice_buffer_size) != 0) {
            fprintf(stderr, "Rank %d: posix_memalign failed for local_y_slices\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (rank == 0 && DEBUG_PRINTS) {
            size_t slice_bytes = (size_t)slice_buffer_size * sizeof(fftw_complex_t);
            printf("[MEMORY] Allocated local_y_slices (reusable): %zu bytes (%.2f GB)\n",
                   slice_bytes, slice_bytes / (1024.0 * 1024.0 * 1024.0));
        }
    }
    
    // ========================================================================
    // MULTI-BATCH LOOP: GENERATION + COMMUNICATION + ACCUMULATION
    // ========================================================================
    
    STimer t_gen, t_comm;
    t_gen.Start();
    
    if (rank == 0) {
        printf("\n[MULTI-BATCH] Processing %d batches...\n", global_max_batches);
    }
    
    for (int batch_idx = 0; batch_idx < global_max_batches; batch_idx++) {
        // ===== BATCH STEP 1: Get this batch's data =====
        int my_batch_slice_count = 0;
        int y_batch_primary = -1, y_batch_mirror = -1;
        bool batch_is_self_conjugate = false;
        
        if (!is_idle_rank && batch_idx < my_num_pairs) {
            y_batch_primary = my_pair_list[batch_idx].y_primary;
            y_batch_mirror = my_pair_list[batch_idx].y_mirror;
            batch_is_self_conjugate = my_pair_list[batch_idx].is_self_conjugate;
            my_batch_slice_count = batch_is_self_conjugate ? 1 : 2;
        }
        
        if (DEBUG_PRINTS && rank < 2 && batch_idx < 2) {
            printf("[BATCH %d] Rank %d: Processing %d Y-slices", 
                   batch_idx, rank, my_batch_slice_count);
            if (my_batch_slice_count > 0) {
                printf(" (Y=%d", y_batch_primary);
                if (!batch_is_self_conjugate) printf(", Y=%d", y_batch_mirror);
                printf(")");
            }
            printf("\n");
        }
        
        // ===== BATCH STEP 2: Generate + 2D FFT =====
        if (!is_idle_rank && my_batch_slice_count > 0) {
            // Clear buffer for reuse
            memset(local_y_slices, 0, sizeof(fftw_complex_t) * 2 * narray * N * N);
            
            // Generate this batch's pair
            fftw_complex_t *primary_ptr = &local_y_slices[0 * narray * N * N];
            fftw_complex_t *conjugate_ptr = batch_is_self_conjugate
                ? &local_y_slices[0 * narray * N * N]
                : &local_y_slices[1 * narray * N * N];
            
            generate_hermitian_slice_pair_local(
                N, y_batch_primary, y_batch_mirror,
                primary_ptr, conjugate_ptr,
                narray, plan_2d, rank,
                ps_params_ptr,  // Legacy power spectrum parameters (NULL = use uniform RNG)
                ps,             // v15.2: zeldovich-PLT PowerSpectrum handle (NULL = use legacy or uniform RNG)
                params          // v15.2: zeldovich-PLT Parameters handle (NULL = use legacy)
            );
            
            // DEBUG: Memory guard after generation + 2D FFT
            if (rank < 4 || DEBUG_PRINTS) {
                fprintf(stderr, "[Rank %d] After generate_hermitian_slice_pair_local (Y=%d): Verifying local_y_slices buffer...\n",
                       rank, y_batch_primary);
                fflush(stderr);
                
                // Verify local_y_slices buffer bounds
                size_t expected_size = (size_t)2 * narray * N * N * sizeof(fftw_complex_t);
                if (local_y_slices != NULL) {
                    // Touch memory to trigger segfault if corrupted
                    volatile char checksum = 0;
                    for (size_t i = 0; i < expected_size && i < 1024*1024; i += 4096) {
                        checksum ^= ((char*)local_y_slices)[i];
                    }
                    fprintf(stderr, "[Rank %d] local_y_slices buffer OK (size=%zu bytes, checksum=%d)\n",
                           rank, expected_size, (int)checksum);
                } else {
                    fprintf(stderr, "[Rank %d] WARNING: local_y_slices is NULL!\n", rank);
                }
                fflush(stderr);
            }
        }
        
        // ===== BATCH STEP 3: Calculate send/recv counts (V11 MODIFIED) =====
        int *sendcounts_batch, *sdispls_batch, *recvcounts_batch, *rdispls_batch;
        int total_send_batch, total_recv_batch;
        
        calculate_batch_send_recv_counts(
            rank, num_ranks, N, narray, batch_idx,
            my_batch_slice_count, my_pencils,
            &sendcounts_batch, &sdispls_batch,
            &recvcounts_batch, &rdispls_batch,
            &total_send_batch, &total_recv_batch
        );
        
        // V11: Adjust rdispls to point into persistent recv_buffer
        // CRITICAL: Check for integer overflow! rdispls_batch is int*, but recv_displs_src is int64_t*
        for (int src = 0; src < num_ranks; src++) {
            int64_t new_displ = recv_displs_src[src] + src_write_cursor[src];
            
            // CHECK FOR OVERFLOW!
            if (new_displ > INT_MAX) {
                fprintf(stderr, "[Rank %d] ERROR: rdispls overflow! src=%d, displ=%lld > INT_MAX (%d)\n",
                       rank, src, (long long)new_displ, INT_MAX);
                fprintf(stderr, "  recv_displs_src[%d] = %lld\n", src, (long long)recv_displs_src[src]);
                fprintf(stderr, "  src_write_cursor[%d] = %d\n", src, src_write_cursor[src]);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            if (new_displ < 0) {
                fprintf(stderr, "[Rank %d] ERROR: Negative rdispls! src=%d, displ=%lld\n",
                       rank, src, (long long)new_displ);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            rdispls_batch[src] = (int)new_displ;
        }
        
        // ===== BATCH STEP 4: Allocate per-batch send buffer (V11: No recv_buffer_batch!) =====
        fftw_complex_t *send_buffer_batch = NULL;
        // V11: recv_buffer_batch REMOVED - receiving directly into persistent recv_buffer
        
        if (total_send_batch > 0) {
            if (posix_memalign((void**)&send_buffer_batch, ALIGN_BYTES,
                               sizeof(fftw_complex_t) * total_send_batch) != 0) {
                fprintf(stderr, "Rank %d: posix_memalign failed for send_buffer_batch\n", rank);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        
        // ===== BATCH STEP 5: Pack =====
        if (!is_idle_rank && my_batch_slice_count > 0 && total_send_batch > 0) {
            pack_slices_to_send_buffer(
                rank, num_ranks, N, narray,
                local_y_slices, my_batch_slice_count, NULL,
                send_buffer_batch, sendcounts_batch, sdispls_batch
            );
        }
        
        // ===== BATCH STEP 6: MPI_Ialltoallv (V11: Direct to recv_buffer) =====
        // DEBUG: Memory guards to detect heap corruption before MPI call
        // Verbose output controlled by VERBOSE_MPI_BUFFER_CHECKS flag - error checks always active
        #if VERBOSE_MPI_BUFFER_CHECKS
        if (rank < 4 || DEBUG_PRINTS) {  // Only print for first few ranks to avoid spam
            fprintf(stderr, "[Rank %d] Verifying buffers before MPI_Ialltoallv...\n", rank);
            fflush(stderr);
            
            // Verify send buffer
            if (send_buffer_batch != NULL && total_send_batch > 0) {
                size_t total_send_bytes = (size_t)total_send_batch * sizeof(fftw_complex_t);
                fprintf(stderr, "[Rank %d] Send buffer: %zu bytes (%d elements)\n", 
                       rank, total_send_bytes, total_send_batch);
                
                // Touch memory to trigger segfault if corrupted
                volatile char checksum = 0;
                for (size_t i = 0; i < total_send_bytes && i < 1024*1024; i += 4096) {  // Limit to 1MB check
                    checksum ^= ((char*)send_buffer_batch)[i];
                }
                fprintf(stderr, "[Rank %d] Send buffer OK (checksum=%d)\n", rank, (int)checksum);
            }
            
            // Verify recv buffer
            if (recv_buffer != NULL) {
                // Calculate total recv size for this batch
                int total_recv_batch = 0;
                for (int i = 0; i < num_ranks; i++) {
                    total_recv_batch += recvcounts_batch[i];
                }
                size_t total_recv_bytes = (size_t)total_recv_batch * sizeof(fftw_complex_t);
                fprintf(stderr, "[Rank %d] Recv buffer: %zu bytes (%d elements) for this batch\n",
                       rank, total_recv_bytes, total_recv_batch);
                fprintf(stderr, "[Rank %d] Recv buffer bounds OK\n", rank);
            }
            
            // CRITICAL: Validate ALL displacements fit in int and don't exceed buffers
            int64_t max_send_displ = 0, max_recv_displ = 0;
            for (int i = 0; i < num_ranks; i++) {
                int64_t send_end = (int64_t)sdispls_batch[i] + sendcounts_batch[i];
                int64_t recv_end = (int64_t)rdispls_batch[i] + recvcounts_batch[i];
                
                if (send_end > max_send_displ) max_send_displ = send_end;
                if (recv_end > max_recv_displ) max_recv_displ = recv_end;
            }
            
            fprintf(stderr, "[Rank %d] Max displacements: send=%lld, recv=%lld\n",
                   rank, (long long)max_send_displ, (long long)max_recv_displ);
            
            // CRITICAL CHECK: Print displacement details
            fprintf(stderr, "[Rank %d] CRITICAL CHECK:\n", rank);
            fprintf(stderr, "  recv_total_elems = %lld\n", (long long)recv_total_elems);
            fprintf(stderr, "  total_send_batch = %d\n", total_send_batch);
            for (int i = 0; i < num_ranks && i < 4; i++) {
                fprintf(stderr, "  src %d: rdispls=%d + recvcounts=%d = %d (end=%lld)\n",
                       i, rdispls_batch[i], recvcounts_batch[i], 
                       rdispls_batch[i] + recvcounts_batch[i],
                       (long long)rdispls_batch[i] + recvcounts_batch[i]);
            }
            
            fprintf(stderr, "[Rank %d] All buffer checks passed, calling MPI_Ialltoallv...\n", rank);
            fflush(stderr);
        }
        #endif
        
        // Silent error checks (always active, no verbose output)
        // Verify sendcounts sum matches total_send_batch
        if (send_buffer_batch != NULL && total_send_batch > 0) {
            int sum_sendcounts = 0;
            for (int i = 0; i < num_ranks; i++) {
                sum_sendcounts += sendcounts_batch[i];
            }
            if (sum_sendcounts != total_send_batch) {
                fprintf(stderr, "[Rank %d] ERROR: sendcounts sum (%d) != total_send_batch (%d)!\n",
                       rank, sum_sendcounts, total_send_batch);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        
        // Verify recvcounts sum matches
        if (recv_buffer != NULL) {
            int total_recv_batch = 0;
            for (int i = 0; i < num_ranks; i++) {
                total_recv_batch += recvcounts_batch[i];
            }
            int sum_recvcounts = 0;
            for (int i = 0; i < num_ranks; i++) {
                sum_recvcounts += recvcounts_batch[i];
            }
            if (sum_recvcounts != total_recv_batch) {
                fprintf(stderr, "[Rank %d] ERROR: recvcounts sum (%d) != total_recv_batch (%d)!\n",
                       rank, sum_recvcounts, total_recv_batch);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            // Verify rdispls don't exceed recv_buffer bounds
            int64_t recv_buffer_size = recv_total_elems;
            for (int i = 0; i < num_ranks; i++) {
                if (recvcounts_batch[i] > 0) {
                    int64_t end_offset = rdispls_batch[i] + recvcounts_batch[i];
                    if (end_offset > recv_buffer_size) {
                        fprintf(stderr, "[Rank %d] ERROR: rdispls[%d] + recvcounts[%d] = %ld exceeds buffer size %ld!\n",
                               rank, i, i, (long)end_offset, (long)recv_buffer_size);
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                }
            }
        }
        
        // CRITICAL: Validate ALL displacements fit in int and don't exceed buffers
        int64_t max_send_displ = 0, max_recv_displ = 0;
        for (int i = 0; i < num_ranks; i++) {
            int64_t send_end = (int64_t)sdispls_batch[i] + sendcounts_batch[i];
            int64_t recv_end = (int64_t)rdispls_batch[i] + recvcounts_batch[i];
            
            if (send_end > max_send_displ) max_send_displ = send_end;
            if (recv_end > max_recv_displ) max_recv_displ = recv_end;
            
            // Check individual displacements
            if (sdispls_batch[i] < 0 || rdispls_batch[i] < 0) {
                fprintf(stderr, "[Rank %d] ERROR: Negative displacement! sdispls[%d]=%d, rdispls[%d]=%d\n",
                       rank, i, sdispls_batch[i], i, rdispls_batch[i]);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            // Check for overflow in individual displacements
            if (sdispls_batch[i] > INT_MAX || rdispls_batch[i] > INT_MAX) {
                fprintf(stderr, "[Rank %d] ERROR: Displacement exceeds INT_MAX! sdispls[%d]=%d, rdispls[%d]=%d\n",
                       rank, i, sdispls_batch[i], i, rdispls_batch[i]);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        
        // Check that displacements don't exceed buffer bounds
        if (max_send_displ > total_send_batch) {
            fprintf(stderr, "[Rank %d] ERROR: Send displacement exceeds buffer! max=%lld > total=%d\n",
                   rank, (long long)max_send_displ, total_send_batch);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (max_recv_displ > recv_total_elems) {
            fprintf(stderr, "[Rank %d] ERROR: Recv displacement exceeds buffer! max=%lld > total=%lld\n",
                   rank, (long long)max_recv_displ, (long long)recv_total_elems);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        t_comm.Start();
        MPI_Request comm_request_batch;
        MPI_Ialltoallv(
            send_buffer_batch, sendcounts_batch, sdispls_batch, MPI_COMPLEX_TYPE,
            recv_buffer, recvcounts_batch, rdispls_batch, MPI_COMPLEX_TYPE,  // ← V11: Changed!
            MPI_COMM_WORLD, &comm_request_batch
        );
        
        // ===== BATCH STEP 7: MPI_Wait (SYNCHRONIZATION POINT - KEEP THIS!) =====
        MPI_Wait(&comm_request_batch, MPI_STATUS_IGNORE);
        t_comm.Stop();
        
        // ===== BATCH STEP 8: Update write cursors (V11: No accumulation) =====
        for (int src = 0; src < num_ranks; src++) {
            src_write_cursor[src] += recvcounts_batch[src];
        }
        
        if (DEBUG_PRINTS && rank == 0 && batch_idx < 2) {
            printf("[V11-BATCH %d] Cursors advanced. Sample: src0=%d, src1=%d\n",
                   batch_idx, src_write_cursor[0], num_ranks > 1 ? src_write_cursor[1] : 0);
        }
        
        // ===== BATCH STEP 9: Free per-batch buffers (V11: No recv_buffer_batch!) =====
        if (send_buffer_batch) free(send_buffer_batch);
        // V11: recv_buffer_batch removed - no longer allocated
        free(sendcounts_batch);
        free(sdispls_batch);
        free(recvcounts_batch);
        free(rdispls_batch);
        
        // Progress indicator
        if (rank == 0 && global_max_batches > 5) {
            if ((batch_idx + 1) % (global_max_batches / 10 + 1) == 0 || batch_idx == global_max_batches - 1) {
                printf("         ... Completed %d/%d batches (%.1f%%)\n",
                       batch_idx + 1, global_max_batches, 
                       100.0 * (batch_idx + 1) / global_max_batches);
            }
        }
    }
    
    t_gen.Stop();
    
    if (rank == 0) {
        printf("[MULTI-BATCH] All batches complete.\n");
        printf("              Generation time: %.6f s\n", t_gen.Elapsed());
        printf("              Communication time: %.6f s\n", t_comm.Elapsed());
    }
    
    // Free local_y_slices (no longer needed)
    if (local_y_slices != NULL) {
        free(local_y_slices);
        local_y_slices = NULL;
    }
    
    // ========================================================================
    // V11: VERIFY CURSORS AFTER ALL BATCHES
    // ========================================================================
    
    if (!is_idle_rank) {
        bool cursor_error = false;
        for (int src = 0; src < num_ranks; src++) {
            int64_t expected = (int64_t)src_total_slices[src] * my_pencils * narray;
            if (src_write_cursor[src] != expected) {
                fprintf(stderr, "[ERROR] Rank %d: src %d cursor mismatch! Got %d, expected %lld\n",
                       rank, src, src_write_cursor[src], (long long)expected);
                cursor_error = true;
            }
        }
        
        if (cursor_error) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        } else if (DEBUG_PRINTS && rank < 2) {
            printf("[V11-VERIFY] Rank %d: All cursors match expected totals\n", rank);
        }
    }
    
    // ========================================================================
    // FREE PCG RNG (NO LONGER NEEDED)
    // ========================================================================
    
    #if !FREE_PCG_AFTER_STAGE1
    cleanup_global_pcg();
    if (rank == 0 && DEBUG_PRINTS) {
        printf("[MEMORY] Freed PCG RNG after batch processing\n");
    }
    #endif
    
    // ========================================================================
    // V11: STAGE 3 - STREAMING X-ROW PROCESSING (FROM SOURCE-GROUPED RECV_BUFFER)
    // ========================================================================
    
    if (rank == 0) {
        printf("\n[V11-Stage 3] Streaming X-row processing from recv_buffer...\n");
    }
    
    // Note: my_bounds and my_pencils already calculated above (before accumulator allocation)
    
    // ========================================================================
    // NEW IN V10: STREAMING X-ROW PROCESSING (FROM ACCUMULATOR)
    // ========================================================================
    // MEMORY: z_count × narray × N × 16 bytes (ONE X-ROW ONLY, e.g., ~4 GB for N=32K, narray=4)
    // Memory reduced by ~80× compared to v9 (storing ONE X-row instead of ALL pencils)
    // Memory layout: [Z][Array][Y] where Z is 0 to z_count-1 (one X-row)
    
    STimer t_streaming;
    t_streaming.Start();
    
    if (rank == 0) {
        printf("\n[Stage 3] Z-slab streaming: Unpack → FFT → Write (Zeldovich format)...\n");
    }
    
    // V12: Allocate local_z_slab for ONE Z-SLAB ONLY (active ranks only)
    // V13: Use appropriate bounds for allocation (padded if enabled, core otherwise)
    int64_t elements_per_z_slab = 0;
    if (!is_idle_rank) {
#if USE_X_PADDING
        int x_count = my_extended_bounds.padded.x_end - my_extended_bounds.padded.x_start;
        int z_count = my_extended_bounds.padded.z_end - my_extended_bounds.padded.z_start;
#else
        int x_count = my_extended_bounds.core.x_end - my_extended_bounds.core.x_start;
        int z_count = my_extended_bounds.core.z_end - my_extended_bounds.core.z_start;
#endif
        elements_per_z_slab = (int64_t)x_count * narray * N;
        
        if (posix_memalign((void**)&local_z_slab, ALIGN_BYTES, 
                           sizeof(fftw_complex_t) * elements_per_z_slab) != 0) {
            fprintf(stderr, "Rank %d: posix_memalign failed for local_z_slab (one Z-slab)\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (rank == 0 && DEBUG_PRINTS) {
            size_t total_bytes = (size_t)elements_per_z_slab * sizeof(fftw_complex_t);
#if USE_X_PADDING
            printf("[V13-MEMORY] Allocated local_z_slab (one Z-slab with padding): %zu bytes (%.2f GB)\n",
                   total_bytes, total_bytes / (1024.0 * 1024.0 * 1024.0));
            printf("         Processing %d Z-slabs sequentially (Z=[%d,%d))\n", 
                   z_count, my_extended_bounds.padded.z_start, my_extended_bounds.padded.z_end);
            printf("         Each Z-slab: %d X-values (PADDED) × %d arrays × %d Y-values\n",
                   x_count, narray, N);
            printf("         Output format: [Z][Array][Y][X] (Zeldovich-compatible)\n");
            printf("         NOTE: X-count includes %d padding for Abacus compatibility\n", X_PADDING);
#else
            printf("[MEMORY] Allocated local_z_slab (one Z-slab, core grid): %zu bytes (%.2f GB)\n",
                   total_bytes, total_bytes / (1024.0 * 1024.0 * 1024.0));
            printf("         Processing %d Z-slabs sequentially (Z=[%d,%d))\n", 
                   z_count, my_extended_bounds.core.z_start, my_extended_bounds.core.z_end);
            printf("         Each Z-slab: %d X-values (CORE) × %d arrays × %d Y-values\n",
                   x_count, narray, N);
            printf("         Output format: [Z][Array][Y][X] (Zeldovich-compatible)\n");
#endif
        }
    } else {
        // Idle ranks: local_z_slab stays NULL
        local_z_slab = NULL;
    }
    
    // Create directory for this rank (before Z-loop)
    char dirname[64];
    snprintf(dirname, sizeof(dirname), "rank_%d", rank);
    int mkdir_result = mkdir(dirname, 0755);
    if (mkdir_result != 0 && errno != EEXIST) {
        fprintf(stderr, "Rank %d: ERROR creating directory %s (errno=%d)\n", rank, dirname, errno);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // V12: MAIN Z-LOOP: Process one Z-slab at a time (Zeldovich-compatible)
    int files_written = 0;
    size_t total_bytes_written = 0;
    
    // Verification: Track max real/imag values per array across all Z-slabs
    real_t *local_max_real = (real_t*)calloc(narray, sizeof(real_t));
    real_t *local_max_imag = (real_t*)calloc(narray, sizeof(real_t));
    
    if (!is_idle_rank && my_pencils > 0) {
        // V13: Use appropriate bounds for Z-loop (padded if enabled, core otherwise)
#if USE_X_PADDING
        int x_count = my_extended_bounds.padded.x_end - my_extended_bounds.padded.x_start;
        int z_count = my_extended_bounds.padded.z_end - my_extended_bounds.padded.z_start;
        
        for (int z = my_extended_bounds.padded.z_start; z < my_extended_bounds.padded.z_end; z++) {
#else
        int x_count = my_extended_bounds.core.x_end - my_extended_bounds.core.x_start;
        int z_count = my_extended_bounds.core.z_end - my_extended_bounds.core.z_start;
        
        for (int z = my_extended_bounds.core.z_start; z < my_extended_bounds.core.z_end; z++) {
#endif
            // Clear buffer for this Z-slab
            memset(local_z_slab, 0, sizeof(fftw_complex_t) * elements_per_z_slab);
            
            // V12: Unpack + FFT this Z-slab from recv_buffer (source-grouped)
            // V14+: Fixed batch-aware unpacking formula
            z_streaming_unpack(
                rank, N, narray,
                z,                               // Current Z to process
#if USE_X_PADDING
                my_extended_bounds.padded,       // V13: MY PADDED grid bounds
#else
                my_extended_bounds.core,        // Core grid bounds (no padding)
#endif
                my_pencils,                      // Total pencils in region
                recv_buffer,                     // Source: persistent recv_buffer
                recv_displs_src,                 // Base offsets per source
                src_total_slices,                // Totals per source (for API)
                y_owner_src,                     // Y → src mapping
                y_src_local_idx,                 // Y → local_idx mapping
                y_batch_idx,                     // Y → batch mapping
                y_slice_idx_in_batch,            // Y → slice_idx in batch mapping
                src_batch_slice_counts,          // [src][batch] → slice count
                global_max_batches,              // Total number of batches
                local_z_slab,                    // Destination buffer
                plan_1d_y                        // FFT plan
            );
            
            // ========== DEBUG: Check imaginary parts BEFORE final verification ==========
            #if DEBUG_PRINTS && !SKIP_VERIFICATION
            if (rank == 0 && z == my_extended_bounds.core.z_start) {
                real_t debug_max_real[4] = {0, 0, 0, 0};
                real_t debug_max_imag[4] = {0, 0, 0, 0};
                for (int x_idx = 0; x_idx < x_count; x_idx++) {
                    for (int array_idx = 0; array_idx < narray; array_idx++) {
                        for (int y = 0; y < N; y++) {
                            double re = fabs_t(ZSLAB(x_idx, array_idx, y, N, narray)[0]);
                            double im = fabs_t(ZSLAB(x_idx, array_idx, y, N, narray)[1]);
                            debug_max_real[array_idx] = fmax_t(debug_max_real[array_idx], re);
                            debug_max_imag[array_idx] = fmax_t(debug_max_imag[array_idx], im);
                        }
                    }
                }
                printf("[DEBUG Z=%d] After 1D FFT - Max imag per array: ", z);
                for (int a = 0; a < narray; a++) {
                    printf("A%d=%.3e ", a, debug_max_imag[a]);
                }
                printf("\n");
            }
            #endif
            
            // ========== VERIFICATION: Check this Z-slab after FFT ==========
            #if !SKIP_VERIFICATION
            for (int x_idx = 0; x_idx < x_count; x_idx++) {
                for (int array_idx = 0; array_idx < narray; array_idx++) {
                    for (int y = 0; y < N; y++) {
                        double re = fabs_t(ZSLAB(x_idx, array_idx, y, N, narray)[0]);
                        double im = fabs_t(ZSLAB(x_idx, array_idx, y, N, narray)[1]);
                        local_max_real[array_idx] = fmax_t(local_max_real[array_idx], re);
                        local_max_imag[array_idx] = fmax_t(local_max_imag[array_idx], im);
                    }
                }
            }
            #endif
            
            // V12: Write this Z-slab in Zeldovich order: [Array][Y][X]
            // Each file contains all (X,Y) for one Z-slab (matches ZeldovichXY output)
            #ifndef SKIP_FILE_WRITE
            char filename[256];
            snprintf(filename, sizeof(filename), "rank_%d/z%d_slab_N%d.bin", 
                    rank, z, N);
            
            FILE *fp = fopen(filename, "wb");
            if (fp) {
                // Write in Zeldovich AZYX order: [Array][Y][X]
                for (int array_idx = 0; array_idx < narray; array_idx++) {
                    for (int y = 0; y < N; y++) {
                        for (int x_idx = 0; x_idx < x_count; x_idx++) {
                            size_t written = fwrite(&ZSLAB(x_idx, array_idx, y, N, narray), 
                                                   sizeof(fftw_complex_t), 1, fp);
                            if (written != 1) {
                                fprintf(stderr, "Rank %d: Write error in %s at (array=%d,y=%d,x_idx=%d)\n",
                                       rank, filename, array_idx, y, x_idx);
                            }
                        }
                    }
                }
                fclose(fp);
                
                files_written++;
                size_t slab_bytes = (size_t)x_count * narray * N * sizeof(fftw_complex_t);
                total_bytes_written += slab_bytes;
                
                if (DEBUG_PRINTS && rank == 0 && files_written <= 3) {
                    printf("[DEBUG] Rank %d: Wrote %s (%d X × %d arrays × %d Y, %.2f MB)\n", 
                           rank, filename, x_count, narray, N, 
                           slab_bytes / (1024.0 * 1024.0));
                    printf("        Layout: [Array][Y][X] (Zeldovich AZYX compatible)\n");
                }
            } else {
                fprintf(stderr, "Rank %d: ERROR opening %s for writing (errno=%d)\n", 
                       rank, filename, errno);
            }
            #else
            // Skip file writing for large N to avoid disk space issues
            files_written++;
            size_t slab_bytes = (size_t)x_count * narray * N * sizeof(fftw_complex_t);
            total_bytes_written += slab_bytes;
            #endif
            
            // Progress indicator for large grids
            if (DEBUG_PRINTS && rank == 0 && z_count > 10) {
#if USE_X_PADDING
                int z_processed = z - my_extended_bounds.padded.z_start + 1;
#else
                int z_processed = z - my_extended_bounds.core.z_start + 1;
#endif
                if (z_processed % (z_count / 10) == 0 || z_processed == z_count) {
                    printf("         ... Processed %d/%d Z-slabs (%.1f%%)\n", 
                           z_processed, z_count, 100.0 * z_processed / z_count);
                }
            }
        }
    }
    // Idle ranks: files_written = 0, total_bytes_written = 0 (already initialized)
    
    // NOTE: In V12 multi-batch mode:
    // - recv_buffer (persistent) is freed in cleanup section (not here)
    // - local_z_slab is freed in cleanup section (not here)
    // - Per-batch send buffers are freed within batch loop
    
    t_streaming.Stop();
    
    // ========================================================================
    // VERIFICATION: Report real-space statistics (after all Z-slabs processed)
    // ========================================================================
    #if !SKIP_VERIFICATION
    {
        // Gather global max for each array (all ranks participate)
        real_t *global_max_real = (real_t*)malloc(sizeof(real_t) * narray);
        real_t *global_max_imag = (real_t*)malloc(sizeof(real_t) * narray);
        MPI_Reduce(local_max_real, global_max_real, narray, MPI_REAL_TYPE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(local_max_imag, global_max_imag, narray, MPI_REAL_TYPE, MPI_MAX, 0, MPI_COMM_WORLD);
        
        if (rank == 0) {
            printf("\n========== VERIFICATION: REAL-SPACE AFTER 3D FFT ==========\n");
            bool all_real = true;
            for (int array_idx = 0; array_idx < narray; array_idx++) {
                real_t ratio = global_max_imag[array_idx] / (global_max_real[array_idx] + 1e-16);
                printf("Array %d: Max real=%.6e, Max imag=%.6e, Ratio=%.6e", 
                       array_idx, global_max_real[array_idx], global_max_imag[array_idx], ratio);
                if (ratio < 1e-7) {
                    printf(" [REAL]\n");
                } else {
                    printf(" [NOT REAL, %.2f%%]\n", 100.0 * ratio);
                    all_real = false;
                }
            }
            if (all_real) {
                printf("RESULT: All %d arrays are PURELY REAL after forward 3D FFT!\n", narray);
            } else {
                printf("RESULT: Some arrays have significant imaginary components\n");
            }
            printf("============================================================\n\n");
        }
        
        free(global_max_real);
        free(global_max_imag);
    }
    #endif
    
    // Free verification arrays
    free(local_max_real);
    free(local_max_imag);
    
    // Gather statistics across ranks (all ranks participate)
    int total_files_written;
    size_t total_bytes_all_ranks;
    MPI_Reduce(&files_written, &total_files_written, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_bytes_written, &total_bytes_all_ranks, 1, MPI_UNSIGNED_LONG, 
               MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("[Stage 3] Streaming complete. Time: %.6f s\n", t_streaming.Elapsed());
        printf("          (Includes unpacking, FFT, and I/O for all Z-slabs)\n");
        printf("          Total files written: %d (across all ranks)\n", total_files_written);
        printf("          Total data written: %.3f GB\n", 
               total_bytes_all_ranks / (1024.0 * 1024.0 * 1024.0));
        
        // Calculate typical file size (varies by rank due to grid decomposition)
        // Each rank writes z_count files (one per Z-slab it owns)
        // Each file contains: x_count X × narray arrays × N Y-values
        int typical_x_count = N / (int)sqrt((double)num_ranks);
        size_t typical_file_size = (size_t)typical_x_count * narray * N * sizeof(fftw_complex_t);
        int typical_files_per_rank = total_files_written / num_ranks;
        printf("          Files per rank: ~%d (z_count per rank, for %d ranks)\n", 
               typical_files_per_rank, num_ranks);
        printf("          Each file contains one Z-slab: x_count X × %d arrays × %d Y-values\n", 
               narray, N);
        printf("          Typical file size: ~%.2f MB (varies by rank's X-extent)\n",
               typical_file_size / (1024.0 * 1024.0));
        printf("          File naming: rank_<rank>/z<z>_slab_N%d.bin\n", N);
        printf("          File format: Binary [Array][Y][X] (Zeldovich AZYX order, fftw_complex = %d bytes)\n", BYTES_PER_COMPLEX);
        printf("          3D FFT COMPLETE (2D XZ + 1D Y)\n\n");
    }
    
    // NOTE: In multi-batch mode, communication buffers are freed per-batch within the batch loop
    // (sendcounts_batch, sdispls_batch, recvcounts_batch, rdispls_batch)
    
    // ========================================================================
    // VERIFICATION: Enabled for Z-slab streaming
    // ========================================================================
    // Verification checks real-space properties as each Z-slab is processed.
    // Statistics are accumulated across all Z-slabs and reported after completion.
    // Verifies that final result is purely real (imaginary parts ≈ 0).
    // ========================================================================
    
    // ========================================================================
    // STAGE 7: INVERSE FFT VERIFICATION (DISABLED)
    // ========================================================================
    // NOTE: Inverse FFT code has been moved to:
    //       hermitian_3d_matrix_mpi_real_inverse_fft.c.archived
    // The forward FFT already produces a purely real result, so inverse
    // FFT verification is not needed for production.
    
    // ========================================================================
    // TIMING SUMMARY
    // ========================================================================
    
    #if DETAILED_TIMING
    if (rank == 0) {
        printf("\n====================================================================================\n");
        printf("TIMING SUMMARY (v11 - Persistent Recv Buffer + Streaming)\n");
        printf("====================================================================================\n");
        printf("Stage 1 (Y-slice generation + 2D FFT): %.6f s\n", t_gen.Elapsed());
        printf("Stage 2 (Metadata exchange):            %.6f s\n", 0.0);  // Minimal time
        printf("Stage 3 (Communication: Alltoallv):    %.6f s\n", t_comm.Elapsed());
        printf("Stage 3 (Streaming: Unpack+FFT+Write):  %.6f s\n", t_streaming.Elapsed());
        printf("------------------------------------------------------------------------------------\n");
        printf("Total 3D FFT time (Gen + Comm + FFT):   %.6f s\n", 
               t_gen.Elapsed() + t_comm.Elapsed());
        printf("Total time (including all stages):      %.6f s\n", 
               t_gen.Elapsed() + t_comm.Elapsed());
        printf("\nNote: In v11, persistent recv_buffer eliminates per-batch allocation\n");
        printf("      Unpacking, 1D FFT, and I/O are combined in streaming stage\n");
        printf("====================================================================================\n");
    }
    #endif
    
    // ========================================================================
    // CLEANUP
    // ========================================================================
    
    // FREE: FFT plans (destroy after all FFTs complete)
    FFTW_DESTROY_PLAN(plan_2d);
    FFTW_DESTROY_PLAN(plan_1d_y);
    
    // FREE: Y-slices (already freed after packing, but check for safety)
    // Note: local_y_slices is freed immediately after packing (Stage 3) to reduce peak memory
    if (!is_idle_rank && local_y_slices != NULL) {
        free(local_y_slices);  // Safety check (should already be NULL)
        local_y_slices = NULL;
    }
    // V11: No y_global_map, recv_buffer_accumulator, or y_to_local_idx to free
    
    // V12: FREE: Z-slab buffer (final data - freed after streaming complete)
    if (!is_idle_rank && local_z_slab != NULL) {
        free(local_z_slab);
        local_z_slab = NULL;
    }
    
    // FREE: Y-mapping arrays (used for unpacking)
    if (y_owner_src != NULL) {
        free(y_owner_src);
        y_owner_src = NULL;
    }
    if (y_src_local_idx != NULL) {
        free(y_src_local_idx);
        y_src_local_idx = NULL;
    }
    if (y_batch_idx != NULL) {
        free(y_batch_idx);
        y_batch_idx = NULL;
    }
    if (y_slice_idx_in_batch != NULL) {
        free(y_slice_idx_in_batch);
        y_slice_idx_in_batch = NULL;
    }
    if (src_batch_slice_counts != NULL) {
        for (int src = 0; src < num_ranks; src++) {
            if (src_batch_slice_counts[src] != NULL) {
                free(src_batch_slice_counts[src]);
            }
        }
        free(src_batch_slice_counts);
        src_batch_slice_counts = NULL;
    }
    if (recv_buffer != NULL) {
        free(recv_buffer);
        recv_buffer = NULL;
    }
    
    // FREE: PCG RNG array (no longer needed after slice generation)
    // FREE: PCG RNG (may already be freed after Stage 1 if FREE_PCG_AFTER_STAGE1=1)
    // If kept for future generation (FREE_PCG_AFTER_STAGE1=0), free it here
    // cleanup_global_pcg() is safe to call even if already freed
    cleanup_global_pcg();
    
    // Cleanup zeldovich-PLT objects
    // PowerSpectrum destructor now has NULL checks to prevent double-free
    // The wrapper tracks destroyed objects to prevent double-free errors
    if (ps) {
        zeldovich_ps_destroy(ps);
        ps = NULL;
    }
    if (params) {
        zeldovich_params_destroy(params);
        params = NULL;
    }
    
    MPI_Finalize();
    
    if (rank == 0) {
        printf("\n====================================================================================\n");
        printf("ALL STAGES COMPLETE!\n");
        printf("====================================================================================\n");
    }
    
    return 0;
}
