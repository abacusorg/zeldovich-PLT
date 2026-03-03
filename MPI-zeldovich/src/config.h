#ifndef HERMITIAN_CONFIG_H
#define HERMITIAN_CONFIG_H

// Modify flags or override them at compile time using -D flags.
// Ex. make CFLAGS="-DUSE_DOUBLE_PRECISION -DDEBUG_PRINTS=0"

// Use X-direction padding with periodic boundary conditions
// 1 = Enable X-padding with periodic wrap-around
#ifndef USE_X_PADDING
#define USE_X_PADDING 0
#endif

// Enable Abacus domain decomposition validation
// 1 = Validate that N is divisible by grid_x and grid_z (exact division required)
#ifndef ENABLE_ABACUS_VALIDATION
#define ENABLE_ABACUS_VALIDATION 0
#endif

// ====================================================================================
// PRECISION
// ====================================================================================
// NOTE:
//       - Single precision (default): Don't define USE_DOUBLE_PRECISION (undefined = single)
//       - Double precision: Define via CFLAGS: make CFLAGS="-DUSE_DOUBLE_PRECISION"
// The Makefile checks CFLAGS to automatically link the correct FFTW libraries.
// DO NOT define USE_DOUBLE_PRECISION here (even as 0) - precision.h uses #ifdef which
// checks for definition, not value. Leave it undefined for single precision!!
// #define USE_DOUBLE_PRECISION  // Uncomment this to use double precision (but prefer CFLAGS)

// ====================================================================================
// Debugging
// ====================================================================================

// Print matrix (small N only)
// 1 = Print matrices for debugging
#ifndef PRINT_MATRICES
#define PRINT_MATRICES 1
#endif

// Print Y-slices before/after communication
// 1 = Print detailed Y-slice information (small N only)
#ifndef PRINT_DETAILED_SLICES
#define PRINT_DETAILED_SLICES 0
#endif

// Print all final output Z-slabs after 3D FFT
// 1 = Print all Z-slabs for debugging (small N only, N <= 16)
#ifndef PRINT_Z_SLABS
#define PRINT_Z_SLABS 1
#endif

// Show detailed debug information
// 1 = Verbose
#ifndef DEBUG_PRINTS
#define DEBUG_PRINTS 1
#endif

// Maximum coordinate value to test for RNG skip debugging (DEBUG_RNG_SKIP)
#ifndef MAX_DEBUG_COORD
#define MAX_DEBUG_COORD 10
#endif

// Maximum coord for RNG consistency debugging
#ifndef MAX_DEBUG_BOUNDARY_COORD
#define MAX_DEBUG_BOUNDARY_COORD 128
#endif

// Sampling stride for RNG debugging when N > 16
// For N > 16, only print coords where x, y, z are multiples of DEBUG_SAMPLE_STRIDE
#ifndef DEBUG_SAMPLE_STRIDE
#define DEBUG_SAMPLE_STRIDE 10
#endif

// Maximum N for full coordinate printing (no sampling)
// For N <= DEBUG_FULL_PRINT_MAX_N, print all coordinates in test range
// For N > DEBUG_FULL_PRINT_MAX_N, use sampling with DEBUG_SAMPLE_STRIDE
#ifndef DEBUG_FULL_PRINT_MAX_N
#define DEBUG_FULL_PRINT_MAX_N 16
#endif

// Debug RNG skip logic: Print detailed skip info to trace RNG state
#ifndef DEBUG_RNG_SKIP
#define DEBUG_RNG_SKIP 1
#endif

// Verbose for MPI buffer verification (send/recv buffer bounds, displacements, etc.)
#ifndef VERBOSE_MPI_BUFFER_CHECKS
#define VERBOSE_MPI_BUFFER_CHECKS 0  // default
#endif

// ====================================================================================
// TIMING CONFIGURATION
// ====================================================================================

// Show timing for each stage
// 0 = No timing output
// 1 = Show detailed timing breakdown
#ifndef DETAILED_TIMING
#define DETAILED_TIMING 1
#endif

// Include I/O timing in measurements
// 0 = Skip I/O timing
// 1 = Include I/O timing
#ifndef TIMING_IO
#define TIMING_IO 1
#endif

// ====================================================================================
// VERIFICATION AND TESTING
// ====================================================================================

// Check Hermitian symmetry during generation
// 0 = Skip Hermitian symmetry checks (production mode)
// 1 = Verify Hermitian symmetry (adds overhead) AND set F=0, H=0 for test
//     When enabled: Sets F=0 and H=0 so conjugate slices are true conjugates of primary slices
//     After 3D FFT, the result should be purely real (imaginary parts ≈ 0)
// 2 = Verify Hermitian symmetry (adds overhead) AND set D=0, G=0 for test
//     When enabled: Sets D=0 and G=0 so the matrix is purely imaginary
//     After 3D FFT, the result should be purely imaginary (real parts ≈ 0)
#ifndef VERIFY_HERMITIAN_SYMMETRY
#define VERIFY_HERMITIAN_SYMMETRY 0
#endif

// Skip verification checks entirely
// 0 = Include verification (recommended for testing)
// 1 = Skip all verification (production mode)
#ifndef SKIP_VERIFICATION
#define SKIP_VERIFICATION 0
#endif

// Verify all Y values present after unpacking
// 0 = Skip Y-filled verification
// 1 = Verify all Y values received (recommended for testing)
#ifndef VERIFY_Y_FILLED
#define VERIFY_Y_FILLED 1
#endif

// Reconstruct global matrix for verification
// 0 = Realistic MPI mode (distributed data, no global reconstruction)
// 1 = Verification mode (reconstruct global matrix, memory intensive)
#ifndef RECONSTRUCT_GLOBAL_FOR_VERIFICATION
#define RECONSTRUCT_GLOBAL_FOR_VERIFICATION 1
#endif

// Verify RNG call counts match expected (similar to zeldovich.cpp assertion)
// 0 = Skip RNG verification (production mode)
// 1 = Verify total RNG calls = 2 * MAX_PPD * MAX_PPD (testing mode)
#ifndef VERIFY_RNG_CALLS
#define VERIFY_RNG_CALLS 0
#endif

// ====================================================================================
// I/O CONFIGURATION
// ====================================================================================

// Write individual pencil files vs global matrix
// 0 = Write global matrix (if RECONSTRUCT_GLOBAL_FOR_VERIFICATION=1)
// 1 = Write individual pencil files per rank
#ifndef WRITE_PENCILS
#define WRITE_PENCILS 0
#endif

// Particle output mode (unified flag for all output paths)
// 0 = Write particle ICs directly (no transpose needed)
//     - Uses WriteParticlesSlab_range with [array][x][y] layout (ZSLAB format)
//     - No transpose overhead, no extra allocation
//     - Works directly with main.cpp's data layout
//     - Writes particle ICs directly
// 1 = Write .bin files for later re-assembly
//     - Writes complex .bin files (rank_*/i*_slab_N*.bin)
//     - Files can be reassembled by write_particles_from_reassembled_mpi.cpp
//     - No particle IC writing in main.cpp
// 2 = Write .bin files then immediately read back and write particle ICs
//     - Writes .bin files first (same as Mode 1)
//     - Then reads them back and calls WriteParticlesSlab_range on local data
//     - Writes particle ICs from .bin file data (useful for verifying .bin format)
#ifndef PARTICLE_OUTPUT_MODE
#define PARTICLE_OUTPUT_MODE 1  // Default: Write .bin files for later re-assembly
#endif

// ====================================================================================
// COMMUNICATION METHOD
// ====================================================================================

// Communication overlap strategy
// 0 = MPI_Ialltoallv
// 1 = Isend/Irecv with overlap
#ifndef OVERLAP_COMMUNICATION
#define OVERLAP_COMMUNICATION 0
#endif

// Metadata exchange method
// 0 = Method A (Each rank calculates + MPI_Allgatherv)
// 1 = Method B (Rank 0 calculates + MPI_Bcast)
#ifndef METADATA_EXCHANGE_METHOD
#define METADATA_EXCHANGE_METHOD 1
#endif

// ====================================================================================
// CONSTANTS
// ====================================================================================

// Maximum grid size for RNG consistency across different N values
// When N < MAX_PPD, skip RNG calls for missing grid points to maintain
// consistency with what a full MAX_PPD * MAX_PPD grid would generate
#ifndef MAX_PPD
#define MAX_PPD 65536
#endif

// Power spectrum spline interpolation resolution
#ifndef SPLINE_RESOLUTION
#define SPLINE_RESOLUTION 128
#endif

// Memory alignment for FFTW
#define ALIGN_BYTES 4096

// FFT direction
// FFTW_FORWARD = -1 (real -> Fourier)
// FFTW_BACKWARD = +1 (Fourier -> real)
#define FFT_SIGN FFTW_BACKWARD

// # arrays per Y-slice
#define NARRAY 4

// X-direction padding
// Only used if USE_X_PADDING=1
#if USE_X_PADDING
#define X_PADDING 10
#else
#define X_PADDING 0  // No padding
#endif

// Free PCG RNG after Stage 1 to save memory?
// 1 = Free PCG RNG after Stage 1 (saves ~1 MB for N=32K)
#define FREE_PCG_AFTER_STAGE1 1

// Phase 4: Parallelize z-loop in Y-slice generation
// 0 = Serial z-loop (original, reference implementation)
// 1 = Parallel z-loop with thread-local RNG copies
#ifndef PARALLELIZE_Z_LOOP
#define PARALLELIZE_Z_LOOP 0
#endif

// When SKIP_FILE_WRITE is defined and non-zero, Stage 3 will skip writing rank_*/z*_slab_N*.bin
// Use -DSKIP_FILE_WRITE=0 to enable file writes, or -USKIP_FILE_WRITE to undefine it
#ifndef SKIP_FILE_WRITE
#define SKIP_FILE_WRITE 1
#elif SKIP_FILE_WRITE == 0
#undef SKIP_FILE_WRITE
#endif

#ifdef PRODUCTION_MODE
    #undef PRINT_MATRICES
    #define PRINT_MATRICES 0
    #undef PRINT_DETAILED_SLICES
    #define PRINT_DETAILED_SLICES 0
    #undef PRINT_Z_SLABS
    #define PRINT_Z_SLABS 0
    #undef DEBUG_PRINTS
    #define DEBUG_PRINTS 0
    #undef DEBUG_RNG_SKIP
    #define DEBUG_RNG_SKIP 0
    #define DEBUG_EIGENVECTOR 0
    #undef VERBOSE_MPI_BUFFER_CHECKS
    #define VERBOSE_MPI_BUFFER_CHECKS 0
    #undef VERIFY_HERMITIAN_SYMMETRY
    #define VERIFY_HERMITIAN_SYMMETRY 0
    #undef VERIFY_Y_FILLED
    #define VERIFY_Y_FILLED 0
    #undef SKIP_VERIFICATION
    #define SKIP_VERIFICATION 1
    #undef DETAILED_TIMING
    #define DETAILED_TIMING 1  // Keep ?
#endif

#ifdef DEBUG_MODE
    #undef PRINT_MATRICES
    #define PRINT_MATRICES 1
    #undef PRINT_DETAILED_SLICES
    #define PRINT_DETAILED_SLICES 1
    #undef DEBUG_PRINTS
    #define DEBUG_PRINTS 1
    #undef VERIFY_HERMITIAN_SYMMETRY
    #define VERIFY_HERMITIAN_SYMMETRY 1
    #undef VERIFY_Y_FILLED
    #define VERIFY_Y_FILLED 1
    #undef SKIP_VERIFICATION
    #define SKIP_VERIFICATION 0
    #undef RECONSTRUCT_GLOBAL_FOR_VERIFICATION
    #define RECONSTRUCT_GLOBAL_FOR_VERIFICATION 1
    #undef VERBOSE_MPI_BUFFER_CHECKS
    #define VERBOSE_MPI_BUFFER_CHECKS 1
#endif

// ====================================================================================
// CONFIGURATION SUMMARY
// ====================================================================================

#define CONFIG_SUMMARY() do { \
    if (rank == 0) { \
        printf("Configuration:\n"); \
        printf("  USE_DOUBLE_PRECISION: %d\n", USE_DOUBLE_PRECISION); \
        printf("  USE_X_PADDING: %d (X_PADDING=%d)\n", USE_X_PADDING, X_PADDING); \
        printf("  DEBUG_PRINTS: %d\n", DEBUG_PRINTS); \
        printf("  VERIFY_HERMITIAN_SYMMETRY: %d\n", VERIFY_HERMITIAN_SYMMETRY); \
        printf("  SKIP_VERIFICATION: %d\n", SKIP_VERIFICATION); \
        printf("  NARRAY: %d\n", NARRAY); \
    } \
} while(0)

#endif 

