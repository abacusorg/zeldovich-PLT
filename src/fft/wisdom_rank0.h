#ifndef WISDOM_RANK0_H
#define WISDOM_RANK0_H

#include "../precision.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Import/export wisdom at FFTW_WISDOM_FILENAME (see fft_wisdom.h),
 * initialize FFTW threads, create 2D many-DFT + 1D DFT plans on plan_buffer (in-place 2D),
 * export accumulated wisdom to the same path.
 *
 * plan_buffer must be posix_memalign(ALIGN_BYTES, narray * N * N * sizeof(fftw_complex_t)).
 * Caller destroys returned plans after this returns.
 */
int wisdom_rank0_plans_and_export(int N, int narray, fftw_complex_t *plan_buffer,
                                  fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out);

#ifdef __cplusplus
}
#endif

#endif
