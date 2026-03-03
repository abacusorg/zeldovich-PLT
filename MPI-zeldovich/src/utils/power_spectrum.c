// ====================================================================================
// POWER SPECTRUM AND GAUSSIAN RANDOM FIELD GENERATION 
// ** obsolete! Use zeldovich_wrapper.cpp instead **
// ====================================================================================

#include "power_spectrum.h"
#include "rng.h"
#include <math.h>
#include <float.h>
#include <stdio.h>

void init_power_spectrum_params(
    power_spectrum_params_t *params,
    double powerlaw_index,
    double normalization,
    double Pk_smooth,
    int fixed_power
) {
    params->mode = POWER_SPECTRUM_POWERLAW;
    params->powerlaw_index = powerlaw_index;
    params->normalization = normalization;
    params->Pk_smooth2 = Pk_smooth * Pk_smooth;  // Store squared for efficiency
    params->fixed_power = fixed_power;
}

double power_spectrum_eval(const power_spectrum_params_t *params, double wavenumber) {
    if (wavenumber <= 0.0) {
        return 0.0;
    }
    
    if (params->mode == POWER_SPECTRUM_POWERLAW) {
        // Power law: P(k) = k^n * normalization * exp(-k^2 * Pk_smooth2)
        double Pk = pow(wavenumber, params->powerlaw_index);
        Pk *= exp(-wavenumber * wavenumber * params->Pk_smooth2);
        Pk *= params->normalization;
        return Pk;
    } else {
        // File-based power spectrum (not yet implemented)
        fprintf(stderr, "ERROR: File-based power spectrum not yet implemented\n");
        return 0.0;
    }
}

// ====================================================================================
// COMPLEX GAUSSIAN RANDOM FIELD GENERATION
// ====================================================================================

// Helper function: Generate random number in (0, 1] from PCG
// Copied from zeldovich-PLT's one_rand<2>()
static double one_rand_pcg(int slice_index) {
    // Get uniform random number from PCG (returns [0, 1))
    double r = random_real_pcg_global(slice_index);
    
    // Shift from [0, 1) to (0, 1] to avoid log(0) in Box-Muller
    // Copied from zeldovich-PLT's one_rand<2>() (see power_spectrum.cpp:one_rand<2>)
    if (r == 0.0) {
        r = 1.0;  // Avoid exactly 0.0 (would cause log(0))
    } else {
        // Add small epsilon to ensure we're in (0, 1] range
        // In practice, PCG won't return exactly 1.0, so this is safe
        r = r + DBL_EPSILON;
        if (r > 1.0) r = 1.0;
    }
    
    return r;
}

void cgauss(
    const power_spectrum_params_t *params,
    double wavenumber,
    int slice_index,
    fftw_complex *result
) {
    // Get power spectrum value at this wavenumber
    double Pk = power_spectrum_eval(params, wavenumber);
    
    // Generate two uniform random numbers (0, 1]
    double R = one_rand_pcg(slice_index);      // Random radius
    double theta = one_rand_pcg(slice_index);  // Random angle
    
    // Box-Muller transform (deterministic version, exactly 2 RNG calls)
    // See zeldovich-PLT's cgauss<2>() implementation
    
    if (!params->fixed_power) {
        // Random amplitude: R = sqrt(-P(k) * log(R))
        // This gives Gaussian distribution with variance P(k)/2 per component
        R = sqrt(-Pk * log(R));
    } else {
        // Fixed amplitude: R = sqrt(P(k))
        // All modes have same amplitude, only phase is random
        R = sqrt(Pk);
    }
    
    // Convert angle to [0, 2π]
    theta = 2.0 * M_PI * theta;
    
    // Get z1, z2 (normal) and store in fftw_complex format
    // fftw_complex is double[2], result is pointer to it
    // Access as (*result)[0] for real, (*result)[1] for imag
    (*result)[0] = R * cos(theta);  // Real part
    (*result)[1] = R * sin(theta);  // Imaginary part
}

