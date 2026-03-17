#ifndef FFT_SETUP_H
#define FFT_SETUP_H

#include "../config.h"
#include "../precision.h"
#include "../types.h"

#ifdef __cplusplus
extern "C" {
#endif

// ====================================================================================
// Create FFTW plans for 2D and 1D transforms
// - plan_buffer: Buffer for 2D batched plan (narray * N * N elements). NULL = allocate dummy.
// - plan_2d_out: Batched plan_many_dft (howmany=narray, NxN each) for Y-slice X-Z plane
// - plan_1d_out: For Y-direction FFT on pencils (size N)
// ====================================================================================

void setup_fftw_plans_full(int N, int narray, fftw_complex_t *plan_buffer,
                           fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out);

#ifdef __cplusplus
}
#endif

#endif 

