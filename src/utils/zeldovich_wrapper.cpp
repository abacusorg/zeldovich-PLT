// ====================================================================================
// C++ wrapper for zeldovich-PLT P(k) and Parameters classes
// ====================================================================================

#include "zeldovich_wrapper.h"
#include <power_spectrum.h>
#include <parameters.h>
#include <zeldovich.h>
#include <cmath>
#include <cstdint>

extern "C" {

// --- Parameters ---

ParametersHandle zeldovich_params_create(const char* param_file) {
    try {
        Parameters* params = new Parameters(fs::path(param_file));
        return static_cast<ParametersHandle>(params);
    } catch (...) {
        return NULL;
    }
}

ParametersHandle zeldovich_params_create_from_buffer(
    const char* header_bytes,
    size_t header_len,
    const char* source_name
) {
    if (header_bytes == NULL || header_len == 0 || source_name == NULL) {
        return NULL;
    }
    try {
        Parameters* params = new Parameters(header_bytes, header_len, fs::path(source_name));
        return static_cast<ParametersHandle>(params);
    } catch (...) {
        return NULL;
    }
}

void zeldovich_params_destroy(ParametersHandle params) {
    if (params) {
        delete static_cast<Parameters*>(params);
    }
}

// Pattern: Cast to Parameters object and attribute
double zeldovich_params_get_fundamental(ParametersHandle params) {
    if (!params) return 0.0;
    // Cast to Parameters object and return fundamental
    Parameters* p = static_cast<Parameters*>(params);
    return p->fundamental;
}

double zeldovich_params_get_boxsize(ParametersHandle params) {
    if (!params) return 0.0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->boxsize;
}

double zeldovich_params_get_Pk_scale(ParametersHandle params) {
    if (!params) return 0.0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->Pk_scale;
}

int64_t zeldovich_params_get_ppd(ParametersHandle params) {
    if (!params) return 0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->ppd;
}

int zeldovich_params_get_cpd(ParametersHandle params) {
    if (!params) return 0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->cpd;
}

int zeldovich_params_get_NumZRanks(ParametersHandle params) {
    if (!params) return 0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->num_z_ranks;
}

int zeldovich_params_get_seed(ParametersHandle params) {
    if (!params) return 0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->seed;
}

double zeldovich_params_get_Pk_powerlaw_index(ParametersHandle params) {
    if (!params) return 1000.0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->Pk_powerlaw_index;
}

const char* zeldovich_params_get_Pk_filename(ParametersHandle params) {
    if (!params) return NULL;
    Parameters* p = static_cast<Parameters*>(params);
    if (p->Pk_filename.empty()) return NULL;
    return p->Pk_filename.c_str();
}

double zeldovich_params_get_f_cluster(ParametersHandle params) {
    if (!params) return 0.0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->f_cluster;
}

double zeldovich_params_get_z_initial(ParametersHandle params) {
    if (!params) return 0.0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->z_initial;
}

int zeldovich_params_get_qPLT(ParametersHandle params) {
    if (!params) return 0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->qPLT;
}

const char* zeldovich_params_get_PLT_filename(ParametersHandle params) {
    if (!params) return NULL;
    Parameters* p = static_cast<Parameters*>(params);
    if (p->PLT_filename.empty()) return NULL;
    return p->PLT_filename.c_str();
}

int zeldovich_params_get_qPLTrescale(ParametersHandle params) {
    if (!params) return 0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->qPLTrescale;
}

double zeldovich_params_get_PLT_target_z(ParametersHandle params) {
    if (!params) return 0.0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->PLT_target_z;
}

const char* zeldovich_params_get_ICFormat(ParametersHandle params) {
    if (!params) return NULL;
    Parameters* p = static_cast<Parameters*>(params);
    return p->ICFormat.c_str();
}

int zeldovich_params_get_qdensity(ParametersHandle params) {
    if (!params) return 0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->qdensity;
}

double zeldovich_params_get_k_cutoff(ParametersHandle params) {
    if (!params) return 1.0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->k_cutoff;
}

int zeldovich_params_get_CornerModes(ParametersHandle params) {
    if (!params) return 0;
    Parameters* p = static_cast<Parameters*>(params);
    return p->CornerModes;
}

// --- PowerSpectrum ---

static const int MAX_TRACKED_OBJECTS = 1024;
static PowerSpectrum* destroyed_objects[MAX_TRACKED_OBJECTS];
static int num_destroyed = 0;

PowerSpectrumHandle zeldovich_ps_create(int n, ParametersHandle params) {
    if (!params) return NULL;
    try {
        Parameters* p = static_cast<Parameters*>(params);
        PowerSpectrum* ps = new PowerSpectrum(n, *p);
        return static_cast<PowerSpectrumHandle>(ps);
    } catch (...) {
        return NULL;
    }
}

void zeldovich_ps_destroy(PowerSpectrumHandle ps) {
    if (!ps) return;
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    for (int i = 0; i < num_destroyed; i++) {
        if (destroyed_objects[i] == p)
            return;
    }
    if (num_destroyed < MAX_TRACKED_OBJECTS)
        destroyed_objects[num_destroyed++] = p;
    try {
        delete p;
    } catch (const std::exception& e) {
        fprintf(stderr, "[ERROR] zeldovich_ps_destroy: %s\n", e.what());
    } catch (...) {
        fprintf(stderr, "[ERROR] zeldovich_ps_destroy: unknown exception\n");
    }
}

int zeldovich_ps_init_powerlaw(PowerSpectrumHandle ps, double powerlaw_index, ParametersHandle params) {
    if (!ps || !params) return -1;
    try {
        PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
        Parameters* par = static_cast<Parameters*>(params);
        return p->InitFromPowerLaw(powerlaw_index, *par);
    } catch (...) {
        return -1;
    }
}

int zeldovich_ps_init_file(PowerSpectrumHandle ps, const char* filename, ParametersHandle params) {
    if (!ps || !params) return -1;
    try {
        PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
        Parameters* par = static_cast<Parameters*>(params);
        return p->InitFromFile(fs::path(filename), *par);
    } catch (...) {
        return -1;
    }
}

double zeldovich_ps_power(PowerSpectrumHandle ps, double wavenumber) {
    if (!ps) return 0.0;
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    return p->power(wavenumber);
}

double zeldovich_ps_one_rand(PowerSpectrumHandle ps, int64_t rng_index) {
    if (!ps) return 0.0;
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    return p->one_rand<2>(rng_index);
}

void zeldovich_ps_cgauss(PowerSpectrumHandle ps, double wavenumber, int64_t rng_index, double* real, double* imag) {
    if (!ps || !real || !imag) return;
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    Complx result = p->cgauss<2>(wavenumber, rng_index);
    *real = result.real();
    *imag = result.imag();
}

void zeldovich_ps_advance_rng(PowerSpectrumHandle ps, ParametersHandle params, int64_t rng_index, int64_t nskip) {
    if (!ps || !params || nskip <= 0) return;
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    Parameters* param = static_cast<Parameters*>(params);
    int64_t ppd_half = param->ppd / 2;
    if (rng_index >= 0 && rng_index < ppd_half && p->v2rng) {
        uint64_t advance_amount = (uint64_t)2 * (uint64_t)nskip;
        p->v2rng[rng_index].advance(advance_amount);
    }
}

// All threads read from the same v2rng[global_y] and write into their own buffers.
void zeldovich_ps_get_rng_copy(PowerSpectrumHandle ps, int64_t rng_index, void* out_rng) {
    if (!ps || !out_rng) return;
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    pcg64 copy = p->get_rng_copy(rng_index);
    new (out_rng) pcg64(copy);
}

size_t zeldovich_ps_rng_buffer_size(void) {
    return sizeof(pcg64);
}

void zeldovich_ps_advance_rng_buffer(void* rng_buf, int64_t nskip) {
    if (!rng_buf || nskip <= 0) return;
    pcg64* rng = static_cast<pcg64*>(rng_buf);
    rng->advance((uint64_t)2 * (uint64_t)nskip);
}

// Match zeldovich one_rand<2>: returns (0,1]
static double one_rand_from_pcg64(pcg64* rng) {
    uint64_t r = (*rng)();
    if (r == UINT64_MAX) return 1.0;
    r += (uint64_t)1;
    return ldexp(static_cast<double>(r), -64);
}

void zeldovich_ps_cgauss_from_buffer(void* rng_buf, PowerSpectrumHandle ps, double wavenumber, double* real, double* imag) {
    if (!rng_buf || !ps || !real || !imag) return;
    pcg64* rng = static_cast<pcg64*>(rng_buf);
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    double Pk = p->power(wavenumber);
    double R = one_rand_from_pcg64(rng);
    double theta = one_rand_from_pcg64(rng);
    if (p->fixed_power)
        R = sqrt(Pk);
    else
        R = sqrt(-Pk * log(R));
    theta = 2 * M_PI * theta;
    *real = R * cos(theta);
    *imag = R * sin(theta);
}

double zeldovich_ps_get_normalization(PowerSpectrumHandle ps) {
    if (!ps) return 0.0;
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    return p->normalization;
}

double zeldovich_ps_get_powerlaw_index(PowerSpectrumHandle ps) {
    if (!ps) return 0.0;
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    return p->powerlaw_index;
}

double zeldovich_ps_get_Pk_smooth2(PowerSpectrumHandle ps) {
    if (!ps) return 0.0;
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    return p->Pk_smooth2;
}

int zeldovich_ps_get_fixed_power(PowerSpectrumHandle ps) {
    if (!ps) return 0;
    PowerSpectrum* p = static_cast<PowerSpectrum*>(ps);
    return p->fixed_power;
}

} 

