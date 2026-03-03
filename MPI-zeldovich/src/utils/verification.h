#ifndef HERMITIAN_VERIFICATION_H
#define HERMITIAN_VERIFICATION_H

// ====================================================================================
// VERIFICATION UTILITIES
// ====================================================================================
// Depends on: config.h, precision.h, types.h, mpi.h
// ====================================================================================

#include "config.h"
#include "precision.h"
#include "types.h"
#include <mpi.h>
#include <stdbool.h>

// ====================================================================================
// Verifying Hermitian symmetry, completeness, and correctness at various stages
// ====================================================================================

// Verify Hermitian symmetry between a pair of Y-slices
// Checks: F(x,y,z) == conj(F(N-x, N-y, N-z))
// is_after_2d_fft: true if checking after 2D FFT, false if before
void verify_hermitian_pair(int rank, int y1, int y2, fftw_complex_t *slice1, fftw_complex_t *slice2, 
                           int N, int array_idx, bool is_after_2d_fft);

// Verify Hermitian symmetry in initial Fourier space (before 2D FFT)
// Checks symmetry within a slice (self-conjugate) or between primary and conjugate slices
void verify_initial_fourier_hermitian_symmetry(int N, fftw_complex_t *slice_y, int global_y, 
                                               int y_mirror, fftw_complex_t *conjugate_slice, int array_idx);

// Verify self-conjugate planes (Y=0, Y=N/2) have zero imaginary parts at special points
// Checks points: (0,0), (N/2,0), (0,N/2), (N/2,N/2)
void verify_self_conjugate_constraints(int N, fftw_complex_t *slice, int global_y, int rank, int array_idx);

// Verify all Y-values [0, N) are present after unpacking (BEFORE 1D FFT)
// Uses filled flags array to track which Y-values were received
// Includes both local and global verification via MPI_Reduce
void verify_pencil_completeness_with_flags(char *y_filled, int pencils_per_rank, int N, int narray, int rank);

// Verify final real-space matrix has imaginary parts ≈ 0 (AFTER full 3D FFT)
// Checks that the final output is numerically real (imaginary parts < 1e-10)
void verify_real_space_symmetry(int N, fftw_complex_t *global_matrix);

#endif

