#ifndef HERMITIAN_PRINTING_H
#define HERMITIAN_PRINTING_H

// ====================================================================================
// PRINTING UTILITIES
// ====================================================================================
// All functions respect PRINT_DETAILED_SLICES and PRINT_MATRICES configuration flags.
// Depends on: config.h, precision.h, types.h
// ====================================================================================

#include "config.h"
#include "precision.h"
#include "types.h"

// ====================================================================================
// Print a Y-slice in Fourier space (after 2D FFT)
// Only prints for small N (N <= 16) and if PRINT_DETAILED_SLICES is enabled
// label for output (e.g., "Primary slice", "Conjugate slice")

void print_y_slice_fourier(int rank, int y_global, fftw_complex_t *slice, int N, int array_idx, const char *label);
// ====================================================================================

// Print 3D matrix for visual inspection (small N only)
// Only prints if PRINT_MATRICES is enabled

void print_3d_matrix_visual(int N, fftw_complex_t *global_matrix, const char* title);
// ====================================================================================

// Print a Z-slab (final output after 3D FFT)
// Only prints for small N (N <= 16) and if PRINT_Z_SLABS is enabled
// rank: MPI rank number
// z: Z-index of this slab
// local_z_slab: Pointer to Z-slab data (layout: [X][Array][Y])
// x_count: Number of X-values in this slab
// N: Grid size in Y-direction
// narray: Number of arrays
// x_start: Starting X-index for this rank's region

void print_z_slab(int rank, int z, fftw_complex_t *local_z_slab, int x_count, int N, int narray, int x_start);
// ====================================================================================

#endif 

