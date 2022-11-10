
class Parameters: public Header, public ParseHeader {
    // This is where we're going to stick all of the control parameters.
    // It is responsible for being able to load from an input param file
    // and to write an output header.
public:
    double boxsize;     // Units for the simulation box
    double Pk_scale;    // A rescaling in case the simulation units 
    // are different from the P(k) input file units
    int64_t ppd;        // The size of the simulation grid to generate
    int cpd;
    long long int np;
    int numblock;    // The number of blocks to divide this into.
    // This must be an even divisor!
    double separation;  // boxsize/ppd
    double fundamental; // 2*PI/boxsize
    double nyquist;    // PI/separation
    double k_cutoff; // the wavenumber above which to not input any power, expressed such that k_max = k_nyquist/k_cutoff.  2 = half nyquist, etc.
    int qdensity;    // If 1, output the density. If 2, output *just* the density (no displacements)
    int qascii;        // If non-zero, output in ASCII
    int qoneslab;    // If >=0, only output this z slab.
    int seed;    // Random number seed
    double Pk_norm;    // The scale to normalize P(k) at, in simulation units!
    double Pk_sigma;    // The normalization at that scale, at the initial redshift!
    double Pk_sigma_ratio;    // An alternative to Pk_sigma: a ratio by which to normalize, instead of a target amplitude
    double f_cluster;  // The fraction of the matter that is clustering.  Usually 1 without neutrinos, and 1 - Omega_Smooth/Omega_M with.
    double Pk_smooth;    // The scale to smooth P(k) at, in simulation units!
    int qPk_fix_to_mean;    // Don't draw the mode amplitude from a Gaussian; use sqrt(P(k)) instead.
    char Pk_filename[200];   // The file name for the P(k) input
    double Pk_powerlaw_index;   // The power law index n for a pure power law P(k) ~ k^n
    char output_dir[1024];   // The file name for the Output
    char density_filename[200];   // The file name for a density file output
    double z_initial;
    HeaderStream * inputstream; // Header stream from which the parameters were read. After instantiation, points to end of header so binary data could potentially be read

    int qonemode; // If non-zero, only use the mode given by one_mode
    int one_mode[3]; // Contains one k-vector to select
    
    int qPLT; // If non-zero, use the Particle Linear Theory modes read from a file
    char PLT_filename[1024]; // file containing PLT eigenmodes
    int qPLTrescale; // If non-zero, rescale the initial amplitudes to match continuum linear theory at PLT_target_z
    double PLT_target_z; // The target redshift for the PLT rescaling
    
    char ICFormat[1024]; // Abacus's expected input format (i.e. our output format)
    
    int AllowDirectIO; // If -DDIRECTIO, need to know if we're on a AllowDirectIO

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
    void print(FILE *fp, const char *datatype);
    void print(FILE *fp);
    // Write a suitable header into the output file
    //
    Parameters(char *inputfile): Header() { 
        // Set default values first
        ppd = 0;    // Illegal
        numblock = 2;    // Ok, but you might not want this!
        boxsize = 0;    // Illegal
        Pk_scale = 1;    // Legal default
        qdensity = 0;    // Legal default
        qascii = 0;    // Legal default
        qoneslab = -1;    // Legal default
        Pk_norm = 0;    // Legal default: Don't renormalize the power spectrum
        Pk_sigma = 0;    // Legal default, but you probably don't want this!
        Pk_sigma_ratio = 0;    // Legal default
        f_cluster = 1;  // Legal default
        Pk_smooth = 0;    // Legal default
        qPk_fix_to_mean = 0; // Legal default
        seed = 0;    // Legal default
        strcpy(Pk_filename,"");   // Must specify Pk file or power law
        Pk_powerlaw_index = NAN;  // Must specify Pk file or power law
        strcpy(density_filename,"density%d");  // Legal default
        qonemode = 0; // Legal default
        memset(one_mode, 0, 3*sizeof(int)); // Legal default
        qPLT = 0; // Legal default
        strcpy(PLT_filename,""); // Legal default
        qPLTrescale = 0; // Legal default
        PLT_target_z = 0.; // Legal default, probably don't want!
        k_cutoff = 1.; // Legal default (corresponds to k_nyquist)
        strcpy(ICFormat,""); // Illegal default
        AllowDirectIO = 0;  // Legal default for most cases
        version = -1;  // All new ICs should use verison 2 (default), but version 1 is available for backwards compatibility
        CornerModes = 0;  // Legal default (no corner modes)
        
        // Read the paramater file values
        register_vars();
        inputstream = new HeaderStream(inputfile);
        ReadHeader(*inputstream);
        
        // Check the validity and compute derived quantities
        int status = setup();
        if(status) {
            std::cout <<"Invalid Parameters given \n";
            exit(1);
        }
    }
    ~Parameters() {
        inputstream->Close();

        if(version == 1){
            for(int64_t i = 0; i < ppd/numblock; i++)
                gsl_rng_free(v1rng[i]);
            delete[] v1rng;
        }
        else{
            delete[] v2rng;
        }
    }

    void register_vars(void) {
        installscalar("BoxSize", boxsize,MUST_DEFINE);
        installscalar("ZD_Pk_scale",Pk_scale,MUST_DEFINE);
        installscalar("NP",np,MUST_DEFINE);
        installscalar("ZD_NumBlock",numblock,MUST_DEFINE);
        installscalar("CPD",cpd,MUST_DEFINE);
        installscalar("ZD_qdensity",qdensity,DONT_CARE);
        installscalar("ZD_qoneslab",qoneslab,DONT_CARE);
        installscalar("ZD_Seed",seed,MUST_DEFINE);
        installscalar("ZD_Pk_norm",Pk_norm,MUST_DEFINE);
        installscalar("ZD_Pk_sigma",Pk_sigma,DONT_CARE);
        installscalar("ZD_Pk_sigma_ratio",Pk_sigma_ratio,DONT_CARE);
        installscalar("ZD_f_cluster",f_cluster,DONT_CARE);
        installscalar("ZD_Pk_smooth",Pk_smooth,MUST_DEFINE);
        installscalar("ZD_qPk_fix_to_mean",qPk_fix_to_mean,DONT_CARE);
        installscalar("ZD_Pk_filename",Pk_filename,DONT_CARE);
        installscalar("ZD_Pk_powerlaw_index",Pk_powerlaw_index,DONT_CARE);
        installscalar("InitialConditionsDirectory",output_dir,MUST_DEFINE);
        installscalar("ZD_density_filename",density_filename,DONT_CARE);
        installscalar("InitialRedshift",z_initial,MUST_DEFINE);
        installscalar("ZD_qonemode",qonemode,DONT_CARE);
        installvector("ZD_one_mode",one_mode,3,1,DONT_CARE);
        installscalar("ZD_qPLT",qPLT,DONT_CARE);
        installscalar("ZD_PLT_filename",PLT_filename,DONT_CARE);
        installscalar("ZD_qPLT_rescale",qPLTrescale,DONT_CARE);
        installscalar("ZD_PLT_target_z",PLT_target_z,DONT_CARE);
        installscalar("ZD_k_cutoff",k_cutoff,DONT_CARE);
        installscalar("ICFormat",ICFormat,MUST_DEFINE);
        installscalar("AllowDirectIO",AllowDirectIO,DONT_CARE);
        installscalar("ZD_Version",version,DONT_CARE);
        installscalar("ZD_CornerModes",CornerModes,DONT_CARE);
    }


};

int Parameters::setup() {
    // Compute any derived quantities.  Look for errors.
    // Return 0 if all is well, 1 if this failed.

    if(version == -1){
        fprintf(stderr,R"(
*** ERROR: ZD_Version was not specified for zeldovich-PLT.  New ICs should
    specify ZD_Version = 2; legacy ICs (pre-November 2019) should use
    ZD_Version = 1 to reproduce the old phases.  Please specify one of
    these in the parameter file.
)");
        exit(1);
    }

    assert(version == 1 || version == 2);

    if(version == 1){
        fprintf(stderr,R"(
*** WARNING: zeldovich-PLT is being invoked with ZD_Version = 1.
    This means that the output phases depend on the ZD_NumBlock tuning parameter,
    so version 1 should only be used for backwards compatibility.  Use ZD_Version = 2
    for new ICs.

)");
    }
    
    ppd = (int64_t) round(cbrt(np));
    fprintf(stderr,"Generating ICs for ppd = %lu\n", ppd);
    assert(ppd*ppd*ppd == np);
    assert(ppd <= MAX_PPD);
    
    // NumBlock is only modified in version 1
    if(version == 1){
        // This is critical for random number synchronization among different ppd
        if(k_cutoff != 1.){
            int numblock_old = numblock;
            numblock = numblock*k_cutoff + .5; // Ensure rounding
            fprintf(stderr,"Note: using k_cutoff=%f means that we are using NumBlock=%d instead of the supplied value of NumBlock=%d\n", k_cutoff, numblock, numblock_old);
        }
    }

    // Check for illegal values
    assert(! (boxsize<=0.0) );
    assert(! (ppd<=0) );
    assert(! (numblock<=0) );
    assert(! (Pk_scale<=0.0) );
    assert(! (Pk_norm<0.0) );

    if((bool)(Pk_sigma > 0) == (bool)(Pk_sigma_ratio > 0)){
        fprintf(stderr, "Must specify exactly one of Pk_sigma or Pk_sigma_ratio!\n");
        exit(1);
    }

    assert(f_cluster > 0. && f_cluster <= 1.);  // Anything outside this range is probably a bug

    // Must specify exactly one of Pk file or power law index
    assert( (bool)(strlen(Pk_filename) > 0) != (bool) !std::isnan(Pk_powerlaw_index) );
    if(!std::isnan(Pk_powerlaw_index))
        assert(Pk_powerlaw_index <= 0);  // technically the code supports blue spectra, but it's more likely input error!
    if(qPLT) assert(! (strlen(PLT_filename)==0) );
    assert(k_cutoff >= 1);
    
    // If using PLT, you probably want an output format with velocities
    if(qPLT)
        assert(strncmp(ICFormat, "RV", 2) == 0);
    
    // Compute derived quantities
    separation = boxsize/ppd;    // Length per grid point
    nyquist = M_PI/separation;
    fundamental = 2.0*M_PI/boxsize;    // The k spacing
    
    if(qonemode)
        fprintf(stderr,"one_mode: %d, %d, %d\n",one_mode[0],one_mode[1],one_mode[2]);

    // Set up multiple RNGs to support parallelism and over-/down-sampling
    
    //seed the rng. a seed of zero uses the current time
    unsigned long int longseed = seed;

    if(version == 1){
        int block = ppd/numblock;
        v1rng = new gsl_rng*[block];

        for (int i = 0; i < block; i++){
            v1rng[i] = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(v1rng[i], longseed+i);
        }
    } else {
        // We'll make ppd/2 independent y-planes
        // But we do so by fast-forwarding the base RNG,
        // so logically this is just a single output stream from the RNG
        v2rng = new pcg64[ppd/2];
        v2rng[0] = pcg64(longseed);
        for(int i = 1; i < ppd/2; i++){
            v2rng[i] = v2rng[i-1];
            // Each plane is ppd^2 complexes
            v2rng[i].advance(2*MAX_PPD*MAX_PPD);
        }
    }
    
    return 0;
}

void Parameters::print(FILE *fp) {
    // Print the results
    time_t t = time(0);
    tm * now = localtime(& t);
    
    WriteHStream(fp, *inputstream);
    fprintf(fp,"#modified by ");
    fprintf(fp,VERSION);
    fprintf(fp," TIME:  ");
    fprintf(fp,"%s", asctime(now));
    fprintf(fp,"\n");
    fprintf(fp,"\n");
    double *ainit;
    ainit = new double;
    *ainit = 1.0/(1+z_initial);
    std::cout<<*ainit<<"\n";
    //fwrite(ainit,sizeof(double),1,fp);
    long long int *np;
    np = new long long int;
    *np =  ppd*ppd*ppd;
    //fwrite(np,sizeof(long long int),1,fp);
    
    return;
}