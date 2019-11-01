#include <stdint.h>

class PowerSpectrum:public SplineFunction {
public:
    PowerSpectrum(int n): SplineFunction(n) {
        is_powerlaw = 0;
        powerlaw_index = NAN;
    };
    int fixed_power;
    int is_powerlaw;
    double powerlaw_index;
    double normalization;
    double Pk_smooth2;   // param.Pk_smooth squared
    double Rnorm;
    double kmax;  // max k in the input PS
    double sigmaR_integrand(double k) {
        double x=k*Rnorm;
        double w;
        if (x<=1e-3) w = 1-x*x/10.0;
        else w = 3.0*(sin(x)-x*cos(x))/x/x/x;
        return 0.5/M_PI/M_PI*k*k*w*w*this->power(k);
    }
    double sigmaR(double R) {
        if(!is_powerlaw){
            double target_prec = 1e-6;  // 1e-12 was causing divergence in certain cases
            double precision = 1.0;
            Rnorm = R;
            double retval = sqrt(Romberg(&PowerSpectrum::sigmaR_integrand,0,10.0,target_prec,&precision));
            if(precision > target_prec){
                printf("Error: actual Romberg integration precision (%g) is greater than the target precision (%g); halting.\n", precision, target_prec);
                exit(1);
            }
            return retval;
        }
        else{
            // Use the analytic power law solution
            // See: http://nbviewer.jupyter.org/gist/lgarrison/7e41ee280c57554e256b834ac5c3f753
            double n = powerlaw_index;
            double retval = 9*pow(R,-n-3)/(2*M_PI*sqrt(M_PI)) * tgamma((3+n)/2.)/(tgamma((2-n)/2.)*(n-3)*(n-1));
            retval = sqrt(retval*normalization);
            return retval;
        }
    }
    // Do Romberg integration with up to MAXITER bisections.
    // Give a precision 'prec'.  Also return the estimated precision in 'obtprec'.

    #define MAXITER 32
    double Romberg(double (PowerSpectrum::*func)(double), double a, double b,
    double prec, double *obtprec) {
        int  jj;
        double h, s, fourtokm1, TT[MAXITER+1][MAXITER+1];

        h = 0.5*(b-a);
        TT[0][1] = h * ( (this->*func)(a) + (this->*func)(b) );
        jj = 0;
        do {
            jj++;
            s = 0;
            #pragma omp parallel for reduction(+:s)
            for(uint64_t k=1; k<=(1ULL<<(jj-1)); k++)
                s += (this->*func)(a + (2*k-1)*h);
            TT[jj][1] = 0.5*TT[jj-1][1] + h * s;
            fourtokm1 = 1;
            for(int k=2; k<=jj; k++) {
                fourtokm1 *=4;
                TT[jj][k] = TT[jj][k-1] + ( TT[jj][k-1] - TT[jj-1][k-1])/(fourtokm1 -1);
            }
            h *= 0.5;
            //printf("TT[%d] = %g\n", jj, TT[jj][jj]);
            if(jj>1 && fabs(TT[jj][jj]-TT[jj-1][jj-1])<prec*fabs(TT[jj][jj])) break;
        } while(jj<MAXITER);

        *obtprec = (TT[jj][jj]-TT[jj-1][jj-1])/TT[jj][jj];
        return TT[jj][jj];
    }
    
    int InitFromFile(char filename[], Parameters& param) {
        // Read the file and load/compile the spline function
        // Return 0 if ok
        // Read the file
        // Rescale the given wavenumbers to match to simulation units
        // Renormalize the power spectrum
        // Divide the power spectrum by the box volume, so that the 
        //    inverse FFT is normalized properly.
        char line[200];
        FILE *fp;
        double k,P;
        int nn;
        printf("Loading power spectrum from file \"%s\"\n", filename);
        fp = fopen(filename,"r");
        if(fp == NULL){
            printf("Power spectrum file \"%s\" not found; exiting.\n", filename);
            exit(1);
        }
        nn=0;

        int jj = 0;
        while (fgets(line,200,fp)!=NULL) {
            jj++;
            if (line[0]=='#') continue;
            sscanf(line,"%lf %lf",&k,&P);
            if (k<0.0) continue;
            if (P<0.0) continue;
            k*=param.Pk_scale;
            if (k>0.0) {
                this->load(log(k),log(P));
            } else {
                this->load(-1e3,log(P));
            }
            kmax = std::max(k,kmax);
            nn++;
        }
        this->spline();

        Normalize(param);
        return 0;
    }

    int InitFromPowerLaw(double _powerlaw_index, Parameters& param){
        assert(!std::isnan(_powerlaw_index));
        powerlaw_index = _powerlaw_index;
        is_powerlaw = 1;
        printf("Initializing power spectrum with power law index %g\n", powerlaw_index);

        Normalize(param);
        return 0;
    }

    void Normalize(Parameters& param){
        Pk_smooth2 = 0.0;
        normalization = 1.0;

        // Might still have to normalize things!
        if (param.Pk_norm>0.0) { // Do a normalization 
            printf("Input sigma(%f) = %f\n", param.Pk_norm, sigmaR(param.Pk_norm));
            normalization = param.Pk_sigma/sigmaR(param.Pk_norm);
            normalization *= normalization;
            printf("Final sigma(%f) = %f\n", param.Pk_norm, sigmaR(param.Pk_norm));
        }
        // Might need to normalize to the box volume.  This is appropriate
        // if the iFFT is like FFTW, i.e., not dividing by N.
        // This is because integrals over k bring in a phase space volume
        // per cell of (2*pi/L)^3
        normalization /= param.boxsize*param.boxsize*param.boxsize;
        Pk_smooth2 = param.Pk_smooth*param.Pk_smooth;
        
        fixed_power = param.qPk_fix_to_mean;
        if (fixed_power)
            printf("Fixing density mode amplitudes to sqrt(P(k))\n");
    }

    double power(double wavenumber) {
        // Evaluate the power spectrum and apply a smoothing.
        // Probably we should force P(k=0) to be 0.
        // The smoothing should be exp(-k^2 sig^2/2) for the density,
        // which is exp(-k^2 sig^2) for the power
        
        if (wavenumber<=0.0)
            return 0.0;

        if(is_powerlaw){
            // power law is raw k^n, times normalization and smoothing
            return std::pow(wavenumber,powerlaw_index)*exp(-wavenumber*wavenumber*this->Pk_smooth2)*normalization;
        } else {
            static bool already_warned = false;
            if (wavenumber > kmax && !already_warned) {
                fprintf(stderr, R"(
*** WARNING: power spectrum spline interpolation was requested
    past the maximum k (%f) that was provided in the input power
    spectrum file.  The extrapolation should be well-behaved, but
    make sure that this was expected.  Provide a power spectrum
    that goes to at least k=10 (or higher if your k_Nyquist demands
    it) to get rid of this warning.

)", kmax);
                already_warned = true;
            }
            return exp(this->val(log(wavenumber))-wavenumber*wavenumber*this->Pk_smooth2)*normalization;
        }
    }

    double one_rand(int64_t i) {
        return gsl_rng_uniform(rng[i]);
    }
    Complx cgauss(double wavenumber, int64_t rng) {
        // Return a gaussian complex deviate scaled to the sqrt of the power
        // Box-Muller, adapted from Numerical Recipes
        // If fixed_power is set, the complex deviate always has amplitude sqrt(P(k))
        // rng tells us which of our 2*ppd^2 RNGs to use

        double Pk = this->power(wavenumber);
        double phase1, phase2, r2;
        // printf("P(%f) = %g\n",wavenumber,Pk);
        do { 
            phase1 = one_rand(rng)*2.0-1.0;
            phase2 = one_rand(rng)*2.0-1.0;
            r2 = phase1*phase1+phase2*phase2;
        } while (!(r2<1.0&&r2>0.0));
        if (fixed_power){
            r2 = sqrt(Pk/r2);
        } else {
            r2 = sqrt(-Pk*log(r2)/r2);   // Drop the factor of 2, so these Gaussians
                                         // have variance of 1/2.
        }
        // printf("cgauss: %f %f\n", phase1*r2, phase2*r2);
        return Complx(phase1*r2,phase2*r2);
    }
};
