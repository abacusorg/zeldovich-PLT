
class Parameters: public Header, public ParseHeader {
    // This is where we're going to stick all of the control parameters.
    // It is responsible for being able to load from an input param file
    // and to write an output header.
public:
    double boxsize;     // Units for the simulation box
    double Pk_scale;    // A rescaling in case the simulation units 
    // are different from the P(k) input file units
    int ppd;        // The size of the simulation grid to generate
    int cpd;
    long long int np;
    int numblock;    // The number of blocks to divide this into.
    // This must be an even divisor!
    double separation;  // boxsize/ppd
    double fundamental; // 2*PI/boxsize
    double nyquist;    // PI/separation
    double k_cutoff; // the wavenumber above which to not input any power, expressed such that k_max = k_nyquist/k_cutoff.  2 = half nyquist, etc.
    int qdensity;    // If non-zero, output the density
    int qascii;        // If non-zero, output in ASCII
    int qnoheader;    // If non-zero, don't attach a header
    int qvelocity;    // If non-zero, include the velocities in the binary output
    int qoneslab;    // If >=0, only output this z slab.
    int seed;    // Random number seed
    double Pk_norm;    // The scale to normalize P(k) at, in simulation units!
    double Pk_sigma;    // The normalization at that scale, at the initial redshift!
    double Pk_smooth;    // The scale to smooth P(k) at, in simulation units!
    int qPk_fix_to_mean;    // Don't draw the mode amplitude from a Gaussian; use sqrt(P(k)) instead.
    char Pk_filename[200];   // The file name for the P(k) input
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
    
    int ramdisk; // If -DDIRECTIO, need to know if we're on a ramdisk
    

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
        qvelocity = 0;    // Legal default
        qnoheader = 0;    // Legal default
        qoneslab = -1;    // Legal default
        Pk_norm = 0;    // Legal default: Don't renormalize the power spectrum
        Pk_sigma = 0;    // Legal default, but you probably don't want this!
        Pk_smooth = 0;    // Legal default
        qPk_fix_to_mean = 0; // Legal default
        seed = 0;    // Legal default
        strcpy(Pk_filename,"");   // Illegal
        strcpy(density_filename,"output.density");  // Legal default
        qonemode = 0; // Legal default
        memset(one_mode, 0, 3*sizeof(int)); // Legal default
        qPLT = 0; // Legal default
        strcpy(PLT_filename,""); // Legal default
        qPLTrescale = 0; // Legal default
        PLT_target_z = 0.; // Legal default, probably don't want!
        k_cutoff = 1.; // Legal default (corresponds to k_nyquist)
        strcpy(ICFormat,""); // Illegal default
        ramdisk = 0;  // Legal default for most cases
        
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

        // Set up RNGs, one per plane per block, to support parallelism
        const gsl_rng_type * T;
        gsl_rng_env_setup();
        T = gsl_rng_mt19937; //The mersene twister
        int block = ppd/numblock;
        rng = new gsl_rng *[block];
        
        //seed the rng. a seed of zero uses the current time
        unsigned long int longseed = seed;
        if (seed == 0){
            longseed = time(0);
            printf("Seed 0 will use current time %lu\n", longseed);
        }
        for (int i = 0; i < block; i++){
            rng[i] = gsl_rng_alloc(T);
            gsl_rng_set(rng[i], longseed+i);
        }
    }
    ~Parameters() {
        inputstream->Close();
    }

    void register_vars(void) {
        installscalar("BoxSize", boxsize,MUST_DEFINE);
        installscalar("ZD_Pk_scale",Pk_scale,MUST_DEFINE);
        installscalar("NP",np,MUST_DEFINE);
        installscalar("ZD_NumBlock",numblock,MUST_DEFINE);
        installscalar("CPD",cpd,MUST_DEFINE);
        installscalar("ZD_qdensity",qdensity,DONT_CARE);
        installscalar("ZD_qnoheader",qnoheader,DONT_CARE);
        installscalar("ZD_qvelocity",qvelocity,DONT_CARE);
        installscalar("ZD_qoneslab",qoneslab,DONT_CARE);
        installscalar("ZD_Seed",seed,MUST_DEFINE);
        installscalar("ZD_Pk_norm",Pk_norm,MUST_DEFINE);
        installscalar("ZD_Pk_sigma",Pk_sigma,MUST_DEFINE);
        installscalar("ZD_Pk_smooth",Pk_smooth,MUST_DEFINE);
        installscalar("ZD_qPk_fix_to_mean",qPk_fix_to_mean,DONT_CARE);
        installscalar("ZD_Pk_filename",Pk_filename,MUST_DEFINE);
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
        installscalar("RamDisk",ramdisk,DONT_CARE);
    }


};

int Parameters::setup() {
    // Compute any derived quantities.  Look for errors.
    // Return 0 if all is well, 1 if this failed.
    
    double npcr = pow(np,1.0/3.0);
    ppd = (long long int) floor(npcr+0.5);
    if(ppd - npcr > .0001){
        printf("ppd: %d\n",ppd);
        printf("npcr: %5.4f\n", npcr);

        assert(ppd - npcr < .0001);
    }
    
    // This is critical for random number synchronization among different ppd
    if(k_cutoff != 1.){
        int numblock_old = numblock;
        numblock = numblock*k_cutoff + .5; // Ensure rounding
        printf("Note: using k_cutoff=%f means that we are using NumBlock=%d instead of the supplied value of NumBlock=%d\n", k_cutoff, numblock, numblock_old);
    }

    // Check for illegal values
    assert(! (boxsize<=0.0) );
    assert(! (ppd<=0) );
    assert(! (numblock<=0) );
    assert(! (Pk_scale<=0.0) );
    assert(! (Pk_norm<0.0) );
    assert(! (Pk_sigma<0.0) );
    assert(! (strlen(Pk_filename)==0) );
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
        printf("one_mode: %d, %d, %d\n",one_mode[0],one_mode[1],one_mode[2]);
    
    return 0;
}

void Parameters::print(FILE *fp, const char *datatype) {
    print(fp);
    return;
}
void Parameters::print(FILE *fp) {
    // Print the results
    time_t t = time(0);
    tm * now = localtime(& t);
    
    WriteHStream(fp, *inputstream);
    fprintf(fp,"#modified by ");
    fprintf(fp,VERSION);
    fprintf(fp," TIME:  ");
    fprintf(fp,asctime(now));
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