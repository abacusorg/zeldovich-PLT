#ifndef FFT_WISDOM_H
#define FFT_WISDOM_H

#include "../precision.h"

// Wisdom file path (relative to cwd), same for wisdom_rank0 and MPI run.
#define FFTW_WISDOM_FILENAME "fftw_wisdom_float.wisdom"

#ifdef __cplusplus
extern "C" {
#endif

// Each rank reads the same wisdom file from disk (no MPI broadcast).
void fft_wisdom_import_from_file(int rank);

// Rank 0 only: writes wisdom to FFTW_WISDOM_FILENAME.
void fft_wisdom_export_rank0(int rank);

#ifdef __cplusplus
}
#endif

#endif
