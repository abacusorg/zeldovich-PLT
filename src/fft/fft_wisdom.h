#ifndef FFT_WISDOM_H
#define FFT_WISDOM_H

#include <mpi.h>
#include "../precision.h"

// Shared filename for FFTW wisdom (precision-dependent via FFTW_* macros).
// For single-precision (default), this will store fftwf_* wisdom.
#define FFTW_WISDOM_FILENAME "fftw_wisdom_float"

#ifdef __cplusplus
extern "C" {
#endif

// Import FFTW wisdom on rank 0 from FFTW_WISDOM_FILENAME, then broadcast
// the wisdom string to all ranks in 'comm'. If the file does not exist,
// all ranks start from a clean "no wisdom" state and proceed to plan
// from scratch.
void fft_wisdom_import_broadcast(int rank, MPI_Comm comm);

// Export FFTW wisdom from rank 0 to FFTW_WISDOM_FILENAME. Other ranks
// do nothing. This creates the file if it did not exist, and overwrites
// it otherwise.
void fft_wisdom_export_rank0(int rank);

#ifdef __cplusplus
}
#endif

#endif

