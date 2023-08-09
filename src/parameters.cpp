#include <cstring>
#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>

#include "zeldovich.h"
#include "parameters.h"


// Write a suitable header into the output file
Parameters::Parameters(char *inputfile): Header() { 
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
    Pk_powerlaw_index = 1000;  // Must specify Pk file or power law
    strcpy(density_filename,"density%d");  // Legal default
    qonemode = 0; // Legal default
    memset(one_mode, 0, 3*sizeof(int)); // Legal default
    qPLT = 0; // Legal default
    strcpy(PLT_filename,""); // Legal default
    qPLTrescale = 0; // Legal default
    PLT_target_z = 0.; // Legal default, probably don't want!
    f_NL = 0.;  // Legal default; no primordial non-Gaussianity
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

Parameters::~Parameters() {
    inputstream->Close();
}

void Parameters::register_vars(void) {
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
    installscalar("ZD_f_NL",f_NL,DONT_CARE);
    installscalar("ICFormat",ICFormat,MUST_DEFINE);
    installscalar("AllowDirectIO",AllowDirectIO,DONT_CARE);
    installscalar("ZD_Version",version,DONT_CARE);
    installscalar("ZD_CornerModes",CornerModes,DONT_CARE);
}

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
    assert( (bool)(strlen(Pk_filename) > 0) != (bool) (Pk_powerlaw_index != 1000) );
    if(Pk_powerlaw_index != 1000)
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
