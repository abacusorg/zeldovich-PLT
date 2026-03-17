// ====================================================================================
// EX: Using cgauss() to generate power-spectrum-weighted random fields
// with variance matching a p(k), similar to zeldovich-PLT's PowerSpectrum::cgauss<2>()
// ====================================================================================

#include "power_spectrum.h"
#include "rng.h"
#include "precision.h" 
#include <math.h>
#include <stdio.h>

void example_usage() {
    int N = 256;
    int slice_index = 5;  // Y-slice index
    
    // Initialize PCG RNG (required before using cgauss)
    initialize_global_pcg(N, N, N, 12345);  // seed = 12345
    
    // Initialize power spectrum parameters
    power_spectrum_params_t ps_params;
    
    // Example: Power law P(k) = k^n with normalization
    double powerlaw_index = -2.0;  // P(k) ~ k^-2
    double normalization = 1.0;    // Adjust amplitude
    double Pk_smooth = 0.0;        // No smoothing
    int fixed_power = 0;           // Random amplitude (not fixed)
    
    init_power_spectrum_params(&ps_params, powerlaw_index, normalization, 
                               Pk_smooth, fixed_power);
    
    // Example: Generate a complex Gaussian mode for a specific wavenumber
    double kx = 10.0, ky = 5.0, kz = 3.0;
    double kmag = sqrt(kx*kx + ky*ky + kz*kz); 
    
    fftw_complex D;
    cgauss(&ps_params, kmag, slice_index, &D);
    
    // D[0] = real part, D[1] = imaginary part
    printf("Generated mode: D = %g + i*%g\n", D[0], D[1]);
    printf("Amplitude: |D| = %g\n", sqrt(D[0]*D[0] + D[1]*D[1]));
    
    // Example: Use in a loop (like zeldovich.cpp)
    for (int z = 0; z < N; z++) {
        for (int x = 0; x < N; x++) {
            // Calculate wavenumber components
            int kx_idx = (x > N/2) ? x - N : x;
            int kz_idx = (z > N/2) ? z - N : z;
            double k2 = kx_idx*kx_idx + ky*ky + kz_idx*kz_idx;
            double kmag = sqrt(k2);
            
            // Generate complex Gaussian mode with variance P(k)
            fftw_complex D;
            cgauss(&ps_params, kmag, slice_index, &D);
            
            // Use D for further processing (e.g., compute displacements)
            // ...
        }
    }
    
    cleanup_global_pcg();
}

