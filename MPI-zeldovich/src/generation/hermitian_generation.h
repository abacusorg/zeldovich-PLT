#ifndef HERMITIAN_GENERATION_H
#define HERMITIAN_GENERATION_H

// ====================================================================================
// HERMITIAN Y-SLICES GENERATION+FFT MODULE
// ====================================================================================
// Depends on: config.h, precision.h, types.h, fft/fft_setup.h, utils/verification.h
// External: PCG RNG, OpenMP, FFTW3
// ====================================================================================

#include "../config.h"
#include "../precision.h"
#include "../types.h"

// External RNG functions (from utils/rng.h)
#include "../utils/rng.h"
// zeldovich-PLT wrapper (for v15.2)
#include "../utils/zeldovich_wrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

// ====================================================================================
// Generate one pair of Hermitian Y-slices (primary + conjugate) 
// with 2D FFT and returns slices as flat arrays
// Arguments:
//      - N: Grid size
//      - global_y: Primary Y-index
//      - y_mirror: Mirror Y-index (N-y for conjugate pairs, y for self-conjugate)
//      - primary_slices: Output buffer for primary slice [narray][N][N]
//      - conjugate_slices: Output buffer for conjugate slice [narray][N][N] (same as primary if self-conjugate)
//      - narray: Number of arrays per slice
//      - plan_2d: Precomputed 2D FFT plan
//      - rank: MPI rank (for debug output only)
//      - ps_handle: zeldovich-PLT PowerSpectrum handle
//      - params_handle: zeldovich-PLT Parameters handle
//
// Uses OpenMP to parallelize X-Z loops within the rank (not recommended)
// Uses global_y for deterministic RNG seeding (thread-safe)
// ====================================================================================

void generate_hermitian_slice_pair_local(
    int N,
    int global_y,
    int y_mirror,
    fftw_complex_t *primary_slices,   // Flat array (all narray arrays for primary slice)
    fftw_complex_t *conjugate_slices, // Flat array (all narray arrays for conjugate slice)
    int narray,                       // Number of arrays per slice
    fftw_plan_t plan_2d,              // 2D FFT plan
    int rank,                         // MPI rank (for debug output)
    PowerSpectrumHandle ps_handle,   // zeldovich-PLT PowerSpectrum handle
    ParametersHandle params_handle,  // zeldovich-PLT Parameters handle
    void** thread_rng_buffers);       // Pre-allocated RNG buffers [nthreads] (NULL = use malloc)

#ifdef __cplusplus
}
#endif

#endif
