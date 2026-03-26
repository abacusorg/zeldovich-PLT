#include "fft_setup.h"
#ifdef USE_FFTW_WISDOM
#include "fft_wisdom.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

void setup_fftw_plans_full(int N, int narray, fftw_complex_t *plan_buffer,
                           fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out)
{
    fftw_complex_t *buf_2d = NULL;
    fftw_complex_t *dummy_1d = NULL;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef USE_FFTW_WISDOM
    // Import FFTW wisdom on rank 0 and broadcast to all ranks so every
    // process starts from an identical planner state
    fft_wisdom_import_broadcast(rank, MPI_COMM_WORLD);
#endif

    // Ensure N^2 * sizeof(fftw_complex_t) is a multiple of 64 for AVX-512 alignment
    // of consecutive arrays in plan_many_dft (N divisible by 4 is sufficient)
    if (N <= 0 || (N % 4) != 0) {
        fprintf(stderr, "[ERROR] N (PPD) must be positive and divisible by 4 for FFT alignment (got N=%d)\n", N);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // ====================================================================================
    // INITIALIZE FFTW THREADING 
    // ====================================================================================
    // FFTW threading must be initialized before creating plans (to use OMP)

    static int fftw_threads_initialized = 0;
    if (!fftw_threads_initialized) {
        int nthreads = omp_get_max_threads();
        
        // All threads for FFTW (no outer OpenMP over narray).
        int fft_threads = nthreads;
        
        if (FFTW_INIT_THREADS() == 0) {
            fprintf(stderr, "[ERROR] Rank %d: Failed to initialize FFTW threads (%s precision)\n", rank, PRECISION_NAME);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        FFTW_PLAN_WITH_NTHREADS(fft_threads);
        if (rank == 0) {
            printf("[FFTW-THREADING] %s precision: %d FFTW threads (no outer OMP over narray)\n", 
                   PRECISION_NAME, fft_threads);
        }
        
        fftw_threads_initialized = 1;
    }
    // ====================================================================================
    
    // Create 2D batched plan (plan_many_dft): howmany=narray transforms of size NxN each.
    // plan_buffer must be provided (caller allocates local_y_slices before setup).
    if (plan_buffer == NULL) {
        fprintf(stderr, "[ERROR] Failed to provide plan_buffer for 2D batched FFT plan\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    buf_2d = plan_buffer;
    
    {// This creates a single FFTW plan that encodes a batch of narray identical 2D transforms, each of size N×N, over contiguous data in buf_2d
        int n[2] = { N, N };
        *plan_2d_out = FFTW_PLAN_MANY_DFT(
            2,      // 2D transform
            n,      // each transform is size N×N
            narray, // number of transforms in the batch
            buf_2d, NULL, 1, N * N, // input data
            buf_2d, NULL, 1, N * N, // output data
            FFT_SIGN, // FFT direction
            FFTW_MEASURE);
    }

    // Verify thread count that FFTW will use when executing this plan (must match plan_with_nthreads)
    if (rank == 0) {
        int planner_n = FFTW_PLANNER_NTHREADS();
        printf("[FFTW-THREADING] Planner nthreads (used at execute): %d\n", planner_n);
        fflush(stdout);
    }
    
    // Create 1D FFT plan with FFTW_MEASURE
    if (posix_memalign((void**)&dummy_1d, ALIGN_BYTES, 
                       sizeof(fftw_complex_t) * N) != 0) {
        fprintf(stderr, "[ERROR] Failed to allocate dummy_1d for 1D FFT plan creation\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    *plan_1d_out = FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, 
                                     FFT_SIGN, FFTW_MEASURE);
    
    free(dummy_1d);
    
    // Check plan creation
    if (*plan_2d_out == NULL || *plan_1d_out == NULL) {
        fprintf(stderr, "[ERROR] Failed to create FFT plans (one or both plans are NULL)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

#ifdef USE_FFTW_WISDOM
    // After successful plan creation (and any warm-up executes outside this
    // function), rank 0 exports the accumulated wisdom back to disk so
    // subsequent runs can reuse it.
    fft_wisdom_export_rank0(rank);
#endif
}

