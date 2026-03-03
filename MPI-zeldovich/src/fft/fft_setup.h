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
// - plan_2d_out: For Y-slice processing (X-Z plane, size N×N)
// - plan_1d_out: For Y-direction FFT on pencils (size N)
// Both plans use FFTW_ESTIMATE for fast planning (change this later)
// ====================================================================================

void setup_fftw_plans_full(int N, int narray, fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out);

#ifdef __cplusplus
}
#endif

#endif 

