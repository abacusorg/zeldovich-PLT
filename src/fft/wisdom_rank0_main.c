/*
 * Standalone binary: ./mpi_zeldovich_wisdom_rank0 <N> <narray>
 * Run before hermitian_3d_matrix with the same OMP_NUM_THREADS; wisdom path is FFTW_WISDOM_FILENAME (cwd).
 */

#include "../config.h"
#include "fft_wisdom.h"
#include "wisdom_rank0.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s <N> <narray>\n",
                argc > 0 ? argv[0] : "mpi_zeldovich_wisdom_rank0");
        fprintf(stderr, "  Use same N, narray, OMP_NUM_THREADS as Zeldovich-MPI; "
                        "run from the directory containing %s.\n",
                FFTW_WISDOM_FILENAME);
        return 1;
    }

    const int N = atoi(argv[1]);
    const int narray = atoi(argv[2]);

    if (N <= 0 || narray <= 0) {
        fprintf(stderr, "%s: N and narray must be positive\n", argv[0]);
        return 1;
    }

    const size_t nbytes = (size_t)narray * (size_t)N * (size_t)N * sizeof(fftw_complex_t);
    fftw_complex_t *plan_buffer = NULL;
    if (posix_memalign((void **)&plan_buffer, ALIGN_BYTES, nbytes) != 0) {
        fprintf(stderr, "%s: posix_memalign failed (%zu bytes)\n", argv[0], nbytes);
        return 1;
    }

    fftw_plan_t plan_2d = NULL;
    fftw_plan_t plan_1d = NULL;

    const int rc = wisdom_rank0_plans_and_export(N, narray, plan_buffer, &plan_2d, &plan_1d);
    if (plan_2d) {
        FFTW_DESTROY_PLAN(plan_2d);
    }
    if (plan_1d) {
        FFTW_DESTROY_PLAN(plan_1d);
    }
    free(plan_buffer);

    return rc == 0 ? 0 : 1;
}
