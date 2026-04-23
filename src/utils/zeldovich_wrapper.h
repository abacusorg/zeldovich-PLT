#ifndef ZELDOVICH_WRAPPER_H
#define ZELDOVICH_WRAPPER_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Opaque pointer type for PowerSpectrum object
typedef void* PowerSpectrumHandle;

// Opaque pointer type for Parameters object
typedef void* ParametersHandle;

// ====================================================================================
// PARAMETERS INTERFACE
// ====================================================================================

// Create Parameters object from parameter file
// Returns NULL on error
ParametersHandle zeldovich_params_create(const char* param_file);

// Create Parameters object from in-memory header bytes
// Buffer must include the ParseHeader expected trailing "\0\0".
// Returns NULL on error.
ParametersHandle zeldovich_params_create_from_buffer(
    const char* header_bytes,
    size_t header_len,
    const char* source_name
);

// Destroy Parameters object
void zeldovich_params_destroy(ParametersHandle params);

// Get fundamental wavenumber (2pi/boxsize)
double zeldovich_params_get_fundamental(ParametersHandle params);

// Get boxsize
double zeldovich_params_get_boxsize(ParametersHandle params);

// Get Pk_scale
double zeldovich_params_get_Pk_scale(ParametersHandle params);

// Get ppd (grid size)
int64_t zeldovich_params_get_ppd(ParametersHandle params);

// Get cpd (coarse particle decomposition; number of slabs for output alignment)
int zeldovich_params_get_cpd(ParametersHandle params);

// Get user-specified number of ranks along z
int zeldovich_params_get_NumZRanks(ParametersHandle params);

// Get seed
int zeldovich_params_get_seed(ParametersHandle params);

// Get Pk_powerlaw_index
double zeldovich_params_get_Pk_powerlaw_index(ParametersHandle params);

// Get Pk_filename (tabular P(k) input). NULL if empty — use power law in that case.
// Valid only while Parameters object exists.
const char* zeldovich_params_get_Pk_filename(ParametersHandle params);

// Get fundamental wavenumber (2pi/BoxSize)
double zeldovich_params_get_fundamental(ParametersHandle params);

// Get f_cluster (fraction of matter that is clustering)
double zeldovich_params_get_f_cluster(ParametersHandle params);

// Get z_initial (initial redshift)
double zeldovich_params_get_z_initial(ParametersHandle params);

// Get qPLT (PLT flag: non-zero if using Particle Linear Theory modes)
int zeldovich_params_get_qPLT(ParametersHandle params);

// Get PLT_filename (file containing PLT eigenmodes)
// Returns: C string (caller should not free). Returns NULL if empty.
// Note: The string is valid only while Parameters object exists.
const char* zeldovich_params_get_PLT_filename(ParametersHandle params);

// Get qPLTrescale (rescaling flag: non-zero to rescale initial amplitudes)
int zeldovich_params_get_qPLTrescale(ParametersHandle params);

// Get PLT_target_z (target redshift for PLT rescaling)
double zeldovich_params_get_PLT_target_z(ParametersHandle params);

// Get ICFormat (output format string, e.g., "RV", "RVDoubleZel", etc.)
// Returns: C string (caller should not free). Returns NULL if empty.
// Note: The string is valid only while Parameters object exists.
const char* zeldovich_params_get_ICFormat(ParametersHandle params);

// Get qdensity (density output mode)
// Returns: 0 = normal mode, 1 = normal + density file, 2 = density only (no displacements)
int zeldovich_params_get_qdensity(ParametersHandle params);

// Get k_cutoff (wavenumber cutoff factor)
// Returns: k_cutoff value (default 1.0, corresponds to k_Nyquist)
double zeldovich_params_get_k_cutoff(ParametersHandle params);

// Get CornerModes (corner mode handling flag)
// Returns: 0 = zero modes with k^2 >= k2_cutoff (default), non-zero = keep corner modes
int zeldovich_params_get_CornerModes(ParametersHandle params);

// ====================================================================================
// POWER SPECTRUM INTERFACE
// ====================================================================================

// Create PowerSpectrum object
PowerSpectrumHandle zeldovich_ps_create(int n, ParametersHandle params);

// Destroy PowerSpectrum object
void zeldovich_ps_destroy(PowerSpectrumHandle ps);

// Initialize from power law
// Returns 0 on success, non-zero on error
int zeldovich_ps_init_powerlaw(PowerSpectrumHandle ps, double powerlaw_index, ParametersHandle params);

// Initialize from file
// Returns 0 on success, non-zero on error
int zeldovich_ps_init_file(PowerSpectrumHandle ps, const char* filename, ParametersHandle params);

// Evaluate power spectrum at given wavenumber
double zeldovich_ps_power(PowerSpectrumHandle ps, double wavenumber);

// Generate random number using zeldovich-PLT's v2rng
double zeldovich_ps_one_rand(PowerSpectrumHandle ps, int64_t rng_index);

// Generate complex Gaussian random number with power spectrum variance
void zeldovich_ps_cgauss(PowerSpectrumHandle ps, double wavenumber, int64_t rng_index, double* real, double* imag);

// Advance RNG for given Y-slice index
void zeldovich_ps_advance_rng(PowerSpectrumHandle ps, ParametersHandle params, int64_t rng_index, int64_t nskip);

// Thread-local RNG support for parallel z-loop
// Get a copy of the RNG for Y-slice rng_index. Caller allocates out_rng (size from zeldovich_ps_rng_buffer_size).
void zeldovich_ps_get_rng_copy(PowerSpectrumHandle ps, int64_t rng_index, void* out_rng);
// Size in bytes for RNG buffer allocation
size_t zeldovich_ps_rng_buffer_size(void);
// Advance RNG in buffer by nskip (units: complex numbers = 2 random numbers each)
void zeldovich_ps_advance_rng_buffer(void* rng_buf, int64_t nskip);
// cgauss using RNG in buffer; needs ps for P(k) and fixed_power
void zeldovich_ps_cgauss_from_buffer(void* rng_buf, PowerSpectrumHandle ps, double wavenumber, double* real, double* imag);

// Get normalization
double zeldovich_ps_get_normalization(PowerSpectrumHandle ps);

// Get powerlaw_index
double zeldovich_ps_get_powerlaw_index(PowerSpectrumHandle ps);

// Get Pk_smooth2
double zeldovich_ps_get_Pk_smooth2(PowerSpectrumHandle ps);

// Get fixed_power flag
int zeldovich_ps_get_fixed_power(PowerSpectrumHandle ps);

#ifdef __cplusplus
}
#endif

#endif

