/*
 * Rank-0 FFTW planning and wisdom export
 */

#include "wisdom_rank0.h"
#include "fft_wisdom.h"
#include "../config.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

// IMPORT WISDOM FROM FILE
static void wisdom_import_file(void)
{
    const int imported = FFTW_IMPORT_WISDOM_FROM_FILENAME(FFTW_WISDOM_FILENAME);
    if (!imported) {
        FFTW_FORGET_WISDOM();
        printf("[wisdom_rank0] No usable wisdom at '%s' (%s); measuring cold.\n",
               FFTW_WISDOM_FILENAME, PRECISION_NAME);
        fflush(stdout);
    } else {
        printf("[wisdom_rank0] Imported wisdom from '%s' (%s precision).\n", FFTW_WISDOM_FILENAME,
               PRECISION_NAME);
        fflush(stdout);
    }
}

// EXPORT WISDOM TO FILE
static int wisdom_export_file(void)
{
    const int ok = FFTW_EXPORT_WISDOM_TO_FILENAME(FFTW_WISDOM_FILENAME);
    if (!ok) {
        fprintf(stderr, "[wisdom_rank0] Failed to export wisdom to '%s' (%s precision)\n",
                FFTW_WISDOM_FILENAME, PRECISION_NAME);
        fflush(stderr);
        return -1;
    }
    printf("[wisdom_rank0] Exported wisdom to '%s' (%s precision)\n", FFTW_WISDOM_FILENAME,
           PRECISION_NAME);
    fflush(stdout);
    return 0;
}

// PLAN AND EXPORT WISDOM
int wisdom_rank0_plans_and_export(int N, int narray, fftw_complex_t *plan_buffer,
                                  fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out)
{
    fftw_complex_t *dummy_1d = NULL; // dummy 1D buffer for 1D FFT plan

    if (plan_2d_out == NULL || plan_1d_out == NULL) {
        fprintf(stderr, "[wisdom_rank0] output plan pointers must be non-NULL\n");
        return -1;
    }
    *plan_2d_out = NULL;
    *plan_1d_out = NULL;

    if (N <= 0 || (N % 4) != 0) {
        fprintf(stderr,
                "[wisdom_rank0] N (PPD) must be divisible by 4 (got N=%d)\n",
                N);
        return -1;
    }
    if (plan_buffer == NULL) {
        fprintf(stderr, "[wisdom_rank0] plan_buffer must be non-NULL\n");
        return -1;
    }

    /* Keep same order as fft_setup.c: thread init before wisdom import. */
    static int fftw_threads_initialized = 0;
    if (!fftw_threads_initialized) {
        if (FFTW_INIT_THREADS() == 0) {
            fprintf(stderr, "[wisdom_rank0] FFTW_INIT_THREADS failed (%s precision)\n",
                    PRECISION_NAME);
            return -1;
        }
        printf("[wisdom_rank0] %s precision: FFTW threads init (2D=OMP_MAX, 1D=1; match fft_setup.c)\n",
               PRECISION_NAME);
        fflush(stdout);
        fftw_threads_initialized = 1;
    }

    wisdom_import_file();

    // 2D plan: use OMP_MAX threads
    FFTW_PLAN_WITH_NTHREADS(omp_get_max_threads());
    printf("[wisdom_rank0] 2D batched plan: FFTW_PLAN_WITH_NTHREADS(%d)\n", omp_get_max_threads());
    fflush(stdout);

    {
        int n[2] = { N, N };
        *plan_2d_out = FFTW_PLAN_MANY_DFT(
            2,
            n,
            narray,
            plan_buffer,
            NULL,
            1,
            N * N,
            plan_buffer,
            NULL,
            1,
            N * N,
            FFT_SIGN,
            FFTW_MEASURE);
    }

    // 1D Y FFT: single FFTW thread (OpenMP parallelizes across pencils; staged buffers in z_streaming)
    FFTW_PLAN_WITH_NTHREADS(1);
    printf("[wisdom_rank0] 1D Y plan: FFTW_PLAN_WITH_NTHREADS(1)\n");
    fflush(stdout);

    if (posix_memalign((void **)&dummy_1d, ALIGN_BYTES, sizeof(fftw_complex_t) * (size_t)N) != 0) {
        fprintf(stderr, "[wisdom_rank0] posix_memalign failed for 1D dummy\n");
        if (*plan_2d_out) {
            FFTW_DESTROY_PLAN(*plan_2d_out);
            *plan_2d_out = NULL;
        }
        return -1;
    }

    *plan_1d_out =
        FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, FFT_SIGN, FFTW_MEASURE);
    free(dummy_1d);

    if (*plan_2d_out == NULL || *plan_1d_out == NULL) {
        fprintf(stderr, "[wisdom_rank0] plan creation failed (NULL plan)\n");
        if (*plan_2d_out) {
            FFTW_DESTROY_PLAN(*plan_2d_out);
            *plan_2d_out = NULL;
        }
        if (*plan_1d_out) {
            FFTW_DESTROY_PLAN(*plan_1d_out);
            *plan_1d_out = NULL;
        }
        return -1;
    }

    if (wisdom_export_file() != 0) {
        FFTW_DESTROY_PLAN(*plan_2d_out);
        FFTW_DESTROY_PLAN(*plan_1d_out);
        *plan_2d_out = NULL;
        *plan_1d_out = NULL;
        return -1;
    }

    return 0;
}
