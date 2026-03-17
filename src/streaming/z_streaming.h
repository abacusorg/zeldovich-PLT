#ifndef Z_STREAMING_H
#define Z_STREAMING_H

// ====================================================================================
// Z-STREAMING MODULE
// ====================================================================================
// Depends on: config.h, precision.h, types.h
// External: FFTW3 (plan passed as parameter)
// ====================================================================================

#include "../config.h"
#include "../precision.h"
#include "../types.h"

#ifdef __cplusplus
extern "C" {
#endif

// ====================================================================================
// Z-slab streaming: Unpack one Z-slab from recv_buffer, apply 1D FFT, ready to write
// V12+: Processes one Z-slab at a time cuz of memory, format for zeldovich-PLT writing 
// V14+: Fixed batch-aware unpacking formula to preserve Hermitian symmetry

// Memory layout: local_z_slab[X][Array][Y] matches Zeldovich AZYX format
// After FFT: Data is ready for WriteParticlesSlab() or file output
//
// NOTE: The batch-aware parameters (y_batch_idx, y_slice_idx_in_batch, src_batch_slice_counts)
// are required to correctly compute offsets in recv_buffer, which is organized by batches
// See docs/RECV_BUFFER_REINDEXING_EXPLANATION.md for details
// ====================================================================================

void z_streaming_unpack(
    int rank, int N, int narray,
    int z_global,                          // Which Z-slab to process
    GridBounds my_bounds,                  // My (X,Z) region bounds
    int my_pencils,                        // Total pencils in my region
    fftw_complex_t *recv_buffer,            // Source data (source-grouped)
    int64_t *recv_displs_src,              // Base offset per source
    int *src_total_slices,                 // Total slices per source (unused, kept for API)
    int *y_owner_src,                      // Y -> src mapping
    int *y_src_local_idx,                  // Y -> local_idx mapping (unused, kept for API)
    int *y_batch_idx,                      // Y -> batch mapping
    int *y_slice_idx_in_batch,             // Y -> slice_idx in batch mapping
    int **src_batch_slice_counts,           // [src][batch] -> slice count
    int global_max_batches,                // Total number of batches
    fftw_complex_t *local_z_slab,          // Destination buffer (one Z-slab)
    fftw_plan_t plan_1d_y);                 // FFT plan for Y-direction

#ifdef __cplusplus
}
#endif

#endif

