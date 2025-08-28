#pragma once

#include <filesystem>

#include "ParseHeader.hh"

namespace fs = std::filesystem;

class Parameters : public ParseHeader {
    // This is where we're going to stick all of the control parameters.
    // It is responsible for being able to load from an input param file
    // and to write an output header.
public:
    double boxsize;   // Units for the simulation box
    double Pk_scale;  // A rescaling in case the simulation units
    // are different from the P(k) input file units
    int64_t ppd;  // The size of the simulation grid to generate
    int cpd;
    long long int np;
    int numblock;  // The number of blocks to divide this into.
    // This must be an even divisor!
    double separation;   // boxsize/ppd
    double fundamental;  // 2*PI/boxsize
    double nyquist;      // PI/separation
    double k_cutoff;     // the wavenumber above which to not input any power, expressed
                      // such that k_max = k_nyquist/k_cutoff.  2 = half nyquist, etc.
    int qdensity;     // If 1, output the density. If 2, output *just* the density (no
                      // displacements)
    int qascii;       // If non-zero, output in ASCII
    int qoneslab;     // If >=0, only output this z slab.
    int seed;         // Random number seed
    double Pk_norm;   // The scale to normalize P(k) at, in simulation units!
    double Pk_sigma;  // The normalization at that scale, at the initial redshift!
    double Pk_sigma_ratio;  // An alternative to Pk_sigma: a ratio by which to
                            // normalize, instead of a target amplitude
    double f_cluster;       // The fraction of the matter that is clustering.  Usually 1
                            // without neutrinos, and 1 - Omega_Smooth/Omega_M with.
    double Pk_smooth;       // The scale to smooth P(k) at, in simulation units!
    int qPk_fix_to_mean;    // Don't draw the mode amplitude from a Gaussian; use
                            // sqrt(P(k)) instead.
    fs::path Pk_filename;  // The file name for the P(k) input
    double Pk_powerlaw_index;  // The power law index n for a pure power law P(k) ~ k^n
    fs::path output_dir;     // The file name for the Output
    fs::path density_filename;  // The file name for a density file output
    double z_initial;
    HeaderStream *inputstream;  // Header stream from which the parameters were read.
                                // After instantiation, points to end of header so
                                // binary data could potentially be read

    int qonemode;     // If non-zero, only use the mode given by one_mode
    std::vector<int> one_mode;  // Contains one k-vector to select

    int qPLT;  // If non-zero, use the Particle Linear Theory modes read from a file
    fs::path PLT_filename;  // file containing PLT eigenmodes
    int qPLTrescale;  // If non-zero, rescale the initial amplitudes to match continuum
                      // linear theory at PLT_target_z
    double PLT_target_z;  // The target redshift for the PLT rescaling

    double f_NL;     // phi^2(x) prefactor for local primordial non-Gaussianity
    double n_s;      // Spectral index of the primordial power (only used for f_NL)
    double Omega_M;  // Omega_M at z=0 (only used for f_NL)

    std::string ICFormat;  // Abacus's expected input format (i.e. our output format)

    int AllowDirectIO;  // If -DDIRECTIO, need to know if we're on a AllowDirectIO

    // Version of the algorithm for getting modes from RNG
    // This directly impacts the phases you get out.
    // Use version = 2 unless you need backwards compatibility with old ICs,
    // in which case use version = 1 (but beware the phases will depend
    // on ZD_NumBlock)
    int version;

    int CornerModes;  // fill modes k > k_Ny. Default: 0.

    int setup();
    // Check against the defaults; complain if something is missing.
    // Setup the other variables
    void print(FILE *fp);
    // Write a suitable header into the output file
    //
    Parameters(const fs::path &inputfile);
    ~Parameters();

    void register_vars(void);
};
