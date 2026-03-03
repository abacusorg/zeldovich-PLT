// ====================================================================================
// VERIFICATION UTILITIES
// ====================================================================================

#include "utils/verification.h"
#include "../mpi_topology.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ====================================================================================
// INTERNAL HELPER FUNCTIONS
// ====================================================================================

// Count missing Y-values in filled flags array
// Returns: Number of missing Y-values across all pencils and arrays on this rank
static int count_missing_y_values(char *y_filled, int pencils_per_rank, int N, int narray)
{
    int missing_count = 0;
    
    for (int pencil = 0; pencil < pencils_per_rank; pencil++) {
        for (int array_idx = 0; array_idx < narray; array_idx++) {
            for (int y = 0; y < N; y++) {
                // Index: pencil * narray * N + array_idx * N + y
                if (!y_filled[pencil * narray * N + array_idx * N + y]) {
                    missing_count++;
                }
            }
        }
    }
    
    return missing_count;
}

// ====================================================================================
// FUNCTION IMPLEMENTATIONS
// ====================================================================================

void verify_hermitian_pair(int rank, int y1, int y2, fftw_complex_t *slice1, fftw_complex_t *slice2, 
                           int N, int array_idx, bool is_after_2d_fft) {
    #if PRINT_DETAILED_SLICES
    if (N > 16) return;
    
    const char *stage = is_after_2d_fft ? "after 2D FFT" : "before 2D FFT";
    printf("\n[RANK %d] Checking Hermitian pair Y=%d and Y=%d, Array=%d (%s):\n", rank, y1, y2, array_idx, stage);
    
    double max_error = 0.0;
    int error_count = 0;
    
    for (int x = 0; x < N; x++) {
        for (int z = 0; z < N; z++) {
            int x_conj = (x == 0) ? 0 : N - x;
            int z_conj = (z == 0) ? 0 : N - z;
            
            double re1 = slice1[x * N + z][0];
            double im1 = slice1[x * N + z][1];
            double re2 = slice2[x_conj * N + z_conj][0];
            double im2 = slice2[x_conj * N + z_conj][1];
            
            // Check: F(x,y,z) == conj(F(N-x, N-y, N-z))
            // So slice1(x,z) at Y=y1 should equal conj(slice2(N-x, N-z)) at Y=y2
            double re_error = fabs_t(re1 - re2);
            double im_error = fabs_t(im1 + im2);  // im1 should equal -im2
            double total_error = sqrt_t(re_error*re_error + im_error*im_error);
            
            if (total_error > 1e-10) {
                error_count++;
                if (error_count <= 10) {  // Show first 10 errors
                    printf("  ERROR at (X=%d,Z=%d,Array=%d): F(%d,%d,%d)=%+.6e%+.6ei, F(%d,%d,%d)*=%+.6e%+.6ei, err=%.3e\n",
                           x, z, array_idx, x, y1, z, re1, im1, x_conj, y2, z_conj, re2, -im2, total_error);
                }
                max_error = (total_error > max_error) ? total_error : max_error;
            }
        }
    }
    
    if (error_count == 0) {
        printf("  Array %d: Hermitian symmetry PERFECT (max error < 1e-10)\n", array_idx);
    } else {
        printf("  Array %d: Hermitian symmetry VIOLATED: %d errors, max error = %.3e\n", array_idx, error_count, max_error);
    }
    #else
    (void)rank; (void)y1; (void)y2; (void)slice1; (void)slice2; (void)N; (void)array_idx; (void)is_after_2d_fft;
    #endif
}

void verify_initial_fourier_hermitian_symmetry(int N, fftw_complex_t *slice_y, int global_y, 
                                               int y_mirror, fftw_complex_t *conjugate_slice, int array_idx) 
{
    int errors = 0;
    const double tol = 1e-10;
    
    if (global_y == y_mirror) {
        // Self-conjugate: Check within the same slice
        for (int x = 0; x < N; x++) {
            int xm = (x == 0) ? 0 : N - x;
            for (int z = 0; z < N; z++) {
                int zm = (z == 0) ? 0 : N - z;
                if (x == xm && z == zm) continue;  // Skip self-symmetric points
                int idx = x * N + z;
                int midx = xm * N + zm;
                double rd = fabs_t(slice_y[idx][0] - slice_y[midx][0]);
                double id = fabs_t(slice_y[idx][1] + slice_y[midx][1]);
                if (rd > tol || id > tol) errors++;
            }
        }
    } else {
        // Conjugate pair: Check between primary and conjugate slices
        for (int x = 0; x < N; x++) {
            int xm = (x == 0) ? 0 : N - x;
            for (int z = 0; z < N; z++) {
                int zm = (z == 0) ? 0 : N - z;
                int idx = x * N + z;
                int midx = xm * N + zm;
                double rd = fabs_t(slice_y[idx][0] - conjugate_slice[midx][0]);
                double id = fabs_t(slice_y[idx][1] + conjugate_slice[midx][1]);
                if (rd > tol || id > tol) errors++;
            }
        }
    }
    
    if (errors > 0 && global_y == 0) {
        printf("[WARN] Y=%d, Array=%d: %d Hermitian symmetry errors (before 2D FFT)\n", global_y, array_idx, errors);
    }
}

void verify_self_conjugate_constraints(int N, fftw_complex_t *slice, int global_y, int rank, int array_idx)
{
    if (global_y != 0 && global_y != N/2) {
        return;  // Not a self-conjugate plane
    }
    
    int violations = 0;
    const double tol = 1e-12;
    
    // Check special points: (0,0), (N/2,0), (0,N/2), (N/2,N/2)
    int special_points[4][2] = {{0, 0}, {N/2, 0}, {0, N/2}, {N/2, N/2}};
    
    for (int i = 0; i < 4; i++) {
        int x = special_points[i][0];
        int z = special_points[i][1];
        int idx = x * N + z;
        
        if (fabs_t(slice[idx][1]) > tol) {
            violations++;
            if (violations <= 3) {  // Show first 3 violations
                printf("[WARN] Rank %d, Y=%d, Array=%d: Point (%d,%d) has non-zero imaginary: %.2e\n", 
                       rank, global_y, array_idx, x, z, slice[idx][1]);
            }
        }
    }
    
    if (violations == 0) {
        if (DEBUG_PRINTS && rank == 0) {
            printf("[VERIFY] Rank %d, Y=%d, Array=%d: Self-conjugate constraints satisfied\n", 
                   rank, global_y, array_idx);
        }
    } else {
        fprintf(stderr, "[ERROR] Rank %d, Y=%d, Array=%d: %d self-conjugate violations!\n", 
                rank, global_y, array_idx, violations);
    }
}

void verify_pencil_completeness_with_flags(char *y_filled, int pencils_per_rank, int N, int narray, int rank)
{
    int missing_count = count_missing_y_values(y_filled, pencils_per_rank, N, narray);
    
    // Show first 5 missing Y-values locally (across all arrays)
    if (missing_count > 0) {
        int shown = 0;
        for (int pencil = 0; pencil < pencils_per_rank && shown < 5; pencil++) {
            for (int array_idx = 0; array_idx < narray && shown < 5; array_idx++) {
                for (int y = 0; y < N && shown < 5; y++) {
                    if (!y_filled[pencil * narray * N + array_idx * N + y]) {
                        printf("[ERROR] Rank %d: Pencil %d, Array %d, Y=%d MISSING\n", 
                               rank, pencil, array_idx, y);
                        shown++;
                    }
                }
            }
        }
    }
    
    // Global verification: Count missing across all ranks
    int global_missing;
    MPI_Reduce(&missing_count, &global_missing, 1, MPI_INT, MPI_SUM, 0, comm_2d);
    
    if (rank == 0) {
        if (global_missing > 0) {
            fprintf(stderr, "[ERROR] Global verification FAILED: %d total (pencil,array,Y) values missing across all ranks!\n", 
                    global_missing);
        } else {
            if (DEBUG_PRINTS) {
                printf("[VERIFY] Global completeness check PASSED :) (all %d arrays × %d Y-values present across all ranks)\n",
                       narray, N);
            }
        }
    }
    
    // Local verification reporting
    if (missing_count == 0) {
        if (DEBUG_PRINTS || rank == 0) {
            printf("[VERIFY] Rank %d: Local completeness check PASSED :) (all %d arrays × %d Y-values present)\n", 
                   rank, narray, N);
        }
    } else {
        fprintf(stderr, "[ERROR] Rank %d: Local completeness check FAILED :( %d (pencil,array,Y) values missing!\n", 
                rank, missing_count);
    }
}

void verify_real_space_symmetry(int N, fftw_complex_t *global_matrix) 
{
    printf("Verifying final real-space imaginary parts...\n");
    int errors = 0;
    double max_im = 0.0, rms = 0.0;
    size_t tot = (size_t)N * (size_t)N * (size_t)N;
    
    for (size_t i = 0; i < tot; i++) {
        double im = fabs_t(global_matrix[i][1]);
        rms += im * im;
        if (im > max_im) max_im = im;
        if (im > 1e-10) { 
            errors++; 
            if (errors < 10) {
                size_t idx = i;
                size_t N2 = (size_t)N * (size_t)N;
                int y = (int)(idx / N2);
                int x = (int)((idx % N2) / (size_t)N);
                int z = (int)(idx % (size_t)N);
                printf("  imag at (%d,%d,%d) = %e\n", x, y, z, global_matrix[i][1]);
            }
        }
    }
    
    rms = sqrt_t(rms / (double)tot);
    printf("  max_imag = %e, rms_imag = %e, errors = %d\n", max_im, rms, errors);
    
    if (errors == 0) {
        printf("  :) OK: Final matrix is (numerically) real\n");
    } else {
        printf("  :( NOT OK: Matrix has non-zero imaginary parts\n");
    }
}
