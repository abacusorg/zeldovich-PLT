#pragma once

#include <stdint.h>

#include <gsl/gsl_rng.h>
#include <pcg-rng/pcg_random.hpp>

#include "parameters.h"
#include "spline_function.h"
#include "zeldovich.h"

class PowerSpectrum : public SplineFunction {
public:
    PowerSpectrum(int n, Parameters &param);
    ~PowerSpectrum();

    int fixed_power;
    int is_powerlaw;
    double powerlaw_index;
    double normalization;
    double Pk_smooth2;  // param.Pk_smooth squared
    double Rnorm;
    double kmax;  // max k in the input PS
    double kmin;  // min (non-zero) k in the input PS
    int block;
    double primordial_norm;
    double n_s;

    gsl_rng **v1rng;  // The random number generators for the deprecated ZD_Version=1
    pcg64 *v2rng;     // The random number generators for version 2 (current version)

    double sigmaR_integrand(double k);
    double sigmaR(double R);

    double Romberg(
       double (PowerSpectrum::*func)(double),
       double a,
       double b,
       double prec,
       double *obtprec
    );

    int InitFromFile(char filename[], Parameters &param);

    int InitFromPowerLaw(double _powerlaw_index, Parameters &param);

    void Normalize(Parameters &param);

    double power(double wavenumber);
    double primordial_power(double wavenumber);
    double infer_Tk(double wavenumber);

    template <int Ver>
    double one_rand(int64_t i);

    template <int Ver>
    Complx cgauss(double wavenumber, int64_t rng);
};
