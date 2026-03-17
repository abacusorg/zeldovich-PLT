#ifndef RNG_H
#define RNG_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Initialize global PCG generators array (NxNxN)
void initialize_global_pcg(int L, int M, int N, uint64_t seed);

// Clean up global PCG generators
void cleanup_global_pcg();

// Generate random # using global v2rng array with advancing
double random_real_pcg_global(int slice_index);

// Advance RNG generator for a specific Y-slice
void advance_pcg_global(int slice_index, uint64_t n);

#ifdef __cplusplus
}
#endif

#endif

