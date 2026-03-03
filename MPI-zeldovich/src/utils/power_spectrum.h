#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H

// ====================================================================================
// POWER SPECTRUM AND GAUSSIAN RNG
// ** obsolete! Use zeldovich_wrapper.cpp instead **
// ====================================================================================

#include <stdint.h>
#include <fftw3.h> 

#ifdef __cplusplus
extern "C" {
#endif

// ====================================================================================
// POWER SPECTRUM FUNCTIONS
// ====================================================================================

// Power spectrum evaluation modes
typedef enum {
    POWER_SPECTRUM_POWERLAW,  // P(k) = k^n * normalization
    POWER_SPECTRUM_FROM_FILE   // P(k) from file (not yet implemented)
} power_spectrum_mode_t;

// Power spectrum parameters structure
typedef struct {
    power_spectrum_mode_t mode;
    double powerlaw_index;      // For power law: P(k) ~ k^n
    double normalization;       // Overall normalization factor
    double Pk_smooth2;          // Smoothing scale squared: exp(-k^2 * Pk_smooth2)
    int fixed_power;            // If 1, use fixed amplitude sqrt(P(k)) instead of random
} power_spectrum_params_t;

// Initialize power spectrum parameters
// powerlaw_index: Power law index n (P(k) ~ k^n)
// normalization: Overall normalization (typically includes box volume factor)
// Pk_smooth: Smoothing scale (Gaussian smoothing: exp(-k^2 * Pk_smooth^2 / 2))
// fixed_power: If 1, use fixed amplitude; if 0, use random Gaussian amplitude
void init_power_spectrum_params(
    power_spectrum_params_t *params,
    double powerlaw_index,
    double normalization,
    double Pk_smooth,
    int fixed_power
);

// Evaluate power spectrum P(k) at given wavenumber
// params: Power spectrum parameters
// wavenumber: Magnitude of wavenumber vector k = sqrt(kx^2 + ky^2 + kz^2)
// Returns: P(k) value
double power_spectrum_eval(const power_spectrum_params_t *params, double wavenumber);

// ====================================================================================
// COMPLEX GAUSSIAN RANDOM FIELD GENERATION
// ====================================================================================

// Generate complex Gaussian random field mode with variance P(k)
// See zeldovich-PLT's PowerSpectrum::cgauss<2>(kmag, y) for details
//
// params: Power spectrum parameters
// wavenumber: Magnitude of wavenumber vector k = sqrt(kx^2 + ky^2 + kz^2)
// slice_index: Y-slice index (used to select which PCG RNG generator to use)
//
// Returns: Complex number D = g1 + i*g2 where:
//   - g1 and g2 are independent Gaussian random variables
//   - Variance of each component: P(k)/2
//   - Total variance: P(k)
//
// Algorithm: Box-Muller transform (deterministic, exactly 2 RNG calls)
//  1. Generate two uniform random numbers R, theta
//  2. Convert to polar coordinates: R = sqrt(-P(k) * log(R)), theta = 2π * theta
//  3. Return: R * (cos(theta) + i*sin(theta))
//
// Note: This function uses the global PCG RNG array from rng.h
// Returns: fftw_complex array [real, imag]
void cgauss(
    const power_spectrum_params_t *params,
    double wavenumber,
    int slice_index,
    fftw_complex *result
);

#ifdef __cplusplus
}
#endif

#endif

