/* v1.0 -- Initial version

v1.1 -- Fixed critical bug in Box-Muller implementation.
Fixed minor bug in spline_function initialization that was causing 
crashes under Linux.

v1.2 -- Added feature to output only one XY slab (at a chosen z)
so as to keep the files smaller for debugging or plotting.

v1.3-- Changed to use ParseHeader to handle input files

v1.4-- Changed output to match input specification for abacus

v1.5-- Changed standard random call to mersene twister from the GSL

v1.6-- Support for "oversampled" simulations (same modes at different PPD) via the k_cutoff option

v1.7-- Support for PLT eigenmodes and rescaling

v2.0-- Support for oversampling with and without new power.
       N.B. The generation of modes from RNG has changed!
       Use "ZD_Version = 1" to get the old phases
       (but beware version 1 phases depend on ZD_NumBlock).
*/

#define VERSION "zeldovich_v2.0"

#include <cmath>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cctype>
#include <cstring>
#include <complex>
#include <cassert>
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <time.h>
#include "spline_function.h"
#include "header.h"
#include "ParseHeader.hh"
#include <omp.h>
#include "pcg-rng/pcg_random.hpp"

#ifdef DIRECTIO
// DIO libraries
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "file.h"
#include "file.cpp"
#include "iolib.cpp"
#endif

#define Complx std::complex<double>

static double __dcube;
#define CUBE(a) ((__dcube=(a))==0.0?0.0:__dcube*__dcube*__dcube)

gsl_rng **v1rng;  // The random number generators for the deprecated ZD_Version=1
pcg64 *v2rng;  // The random number generators for version 2 (current version)

// The PLT eigenmodes
double* eig_vecs;
int64_t eig_vecs_ppd;

// Global maxima of the particle displacements
double max_disp[3];

#define MAX_PPD ((int64_t) 65536)

#include "parameters.cpp"
#include "power_spectrum.cpp"
#include "block_array.cpp"
#include "output.cpp"

// ===============================================================

// TODO: Replace with our own FFT
#include "fftw3.h"
fftw_plan plan1d, plan2d;
void Setup_FFTW(int n) {
    STimer FFTplanning;

    FFTplanning.Start();
    char wisdom_file[1024];
    int wisdom_exists;
    sprintf(wisdom_file, "fftw_zeldovich.wisdom"); 
        // TODO: Where to put this file?  param.output_dir?  Or in the repository?
        // Currently in the repository
    wisdom_exists = fftw_import_wisdom_from_filename(wisdom_file);
    fprintf(stderr,"FFTW Wisdom import from file \"%s\" returned %d (%s).\n", wisdom_file, wisdom_exists, wisdom_exists == 1 ? "success" : "failure");

    fftw_complex *p;
    // p = new fftw_complex[n*n];
    int ret = posix_memalign((void **)&p, 4096, n*n*sizeof(Complx));
    plan1d = fftw_plan_dft_1d(n, p, p, +1, FFTW_PATIENT);
    plan2d = fftw_plan_dft_2d(n, n, p, p, +1, FFTW_PATIENT);
    free(p);
    // delete []p;
    ret = fftw_export_wisdom_to_filename(wisdom_file);
    fprintf(stderr,"FFTW Wisdom export to file \"%s\" returned %d.\n", wisdom_file, ret);
    FFTplanning.Stop();
    fprintf(stderr,"Creating FFTW plans done in %f sec.\n", FFTplanning.Elapsed());
}

// TODO: This handling of the memory allocation required for plans 
// doesn't seem exceptionally safe.  Is it really true that "new Complx"
// always returns a suitably aligned variable? 
// And we're trusting that every 2-d plane ends up aligned, although
// this has a better chance of being ok because we require an even PPD.

void Inverse1dFFT(Complx *p, int n) {
    // Given a pointer to a 1d complex vector, packed as p[n].
    // Do the 1d inverse FFT in place
    fftw_execute_dft(plan1d, (fftw_complex *)p, (fftw_complex *)p);
}
void Inverse2dFFT(Complx *p, int n) {
    // Given a pointer to a 2d complex array, contiguously packed as p[n][n].
    // Do the 2d inverse FFT in place
    fftw_execute_dft(plan2d, (fftw_complex *)p, (fftw_complex *)p);
}
void InverseFFT_Yonly(Complx *p, int n) {
    // Given a pointer to a 2d complex array, contiguously packed as p[n][n].
    // Do the 1d inverse FFT on the first index (the long stride one)
    // for each value of the second index.
    // Note that the Y in the title doesn't refer to the Y direction in 
    // our 3-d problem!
    Complx *tmp;
    int j,k;
    int ret = posix_memalign((void **)&tmp, 4096, n*sizeof(Complx));
    assert(ret==0);

    // tmp = new Complx[n];

    for (j=0;j<n;j++) {
        // We will load one row at a time
        for (k=0;k<n;k++) tmp[k] = p[k*n+j];
        Inverse1dFFT(tmp, n);
        for (k=0;k<n;k++) p[k*n+j] = tmp[k];
    }
    free(tmp);
    //delete []tmp;
}

//================================================================

// We use a set of X-Z arrays of Complx numbers (ordered by A and Y).
// array.ppd is 64-bit, so this expression should be safe from overflow
#define AYZX(_slab,_a,_y,_z,_x) _slab[(int64_t)(_x)+array.ppd*((_z)+array.ppd*((_a)+array.narray*(_y)))]

typedef struct {
    double vec[3];
    double val;
} eigenmode;

void interp_eigmode(int ikx, int iky, int ikz, int64_t ppd, double *e){
#define EIGMODE(_kx,_ky,_kz,_i) (eig_vecs[(int64_t)(_kx)*eig_vecs_ppd*halfppd*4 + (_ky)*halfppd*4 + (_kz)*4 + (_i)])
    int64_t halfppd = eig_vecs_ppd/2 + 1;
    int64_t ppdhalf = eig_vecs_ppd/2;
    if(eig_vecs_ppd % ppd == 0){
        for(int i = 0; i < 4; i++)
            e[i] = EIGMODE(ikx*eig_vecs_ppd/ppd, iky*eig_vecs_ppd/ppd, ikz*eig_vecs_ppd/ppd, i);
        return;
    }
    
    double fx = ((double) eig_vecs_ppd) / ppd * ikx;
    double fy = ((double) eig_vecs_ppd) / ppd * iky;
    double fz = ((double) eig_vecs_ppd) / ppd * ikz;
    
    // For ppd 64, [0,32] are positive k, [33,63] are negative
    // So don't interpolate between 32-33!  Map upwards instead.
    // if(fx > eig_vecs_ppd/2 && fx < eig_vecs_ppd/2 + 1)
    // if(fy > eig_vecs_ppd/2 && fy < eig_vecs_ppd/2 + 1)
    // if(fz > eig_vecs_ppd/2 && fz < eig_vecs_ppd/2 + 1)
    if(fx > ppdhalf && fx < halfppd)
        fx = floor(fx+1);
    if(fy > ppdhalf && fy < halfppd)
        fy = floor(fy+1);
    if(fz > ppdhalf && fz < halfppd)
        fz = floor(fz+1);
    
    // Build the indices of the nearest grid points
    int ikx_l = (int) fx;
    int ikx_h = ikx_l + 1;  // This is okay when ikx_l == eig_vecs_ppd/2 because fx is an integer so ikx_h is never used
    int iky_l = (int) fy;
    int iky_h = iky_l + 1;
    int ikz_l = (int) fz;
    int ikz_h = ikz_l + 1;
    
    // If ikx = 127, then kx = -1, so we should interpolate
    // between -1 and 0; i.e. between ikx = 63 and 0.
    if(ikx_h == eig_vecs_ppd) ikx_h = 0;
    if(iky_h == eig_vecs_ppd) iky_h = 0;
    if(ikz_h == eig_vecs_ppd) ikz_h = 0;
    
    // The fractional position between the grid points
    fx -= ikx_l;
    fy -= iky_l;
    fz -= ikz_l;
    
    // Trilinear interpolation coefficients
    double f[8];
    f[0] = (1 - fx) * (1 - fy) * (1 - fz);
    f[1] = (1 - fx) * (1 - fy) * (fz);
    f[2] = (1 - fx) * (fy) * (1 - fz);
    f[3] = (1 - fx) * (fy) * (fz);
    f[4] = (fx) * (1 - fy) * (1 - fz);
    f[5] = (fx) * (1 - fy) * (fz); 
    f[6] = (fx) * (fy) * (1 - fz);
    f[7] = (fx) * (fy) * (fz);

    // Treat the eigenmode struct as 4 doubles
    for(int i = 0; i < 4; i++){
        e[i] = f[0]*EIGMODE(ikx_l, iky_l, ikz_l, i) + f[1]*EIGMODE(ikx_l, iky_l, ikz_h, i) +
               f[2]*EIGMODE(ikx_l, iky_h, ikz_l, i) + f[3]*EIGMODE(ikx_l, iky_h, ikz_h, i) +
               f[4]*EIGMODE(ikx_h, iky_l, ikz_l, i) + f[5]*EIGMODE(ikx_h, iky_l, ikz_h, i) +
               f[6]*EIGMODE(ikx_h, iky_h, ikz_l, i) + f[7]*EIGMODE(ikx_h, iky_h, ikz_h, i);
    }
}

eigenmode get_eigenmode(int kx, int ky, int kz, int64_t ppd, int qPLT){
    eigenmode e;
    
    if(qPLT){
        // undo nyquist wrapping
        // These are the necessary array indices
        int ikx = kx < 0 ? ppd + kx : kx;
        int iky = ky < 0 ? ppd + ky : ky;
        int ikz = kz < 0 ? ppd + kz : kz;
        // note: np.fft has the convention of freq[ppd/2] = -ppd/2, instead of +ppd/2
        // This is different from the convention in this code
        // but we normalize to this convention when we generate the eigenmodes.
        // kz is already okay because of rfft.
        ikz = ikz > ppd/2 ? ppd - ikz : ikz; // Use the index from the +k half-space
        double k2 = kx*kx + ky*ky + kz*kz;
        
        // Use interpolation to get the eigenmode
        eigenmode ehat;
        assert(sizeof(ehat) == sizeof(double)*4);
        interp_eigmode(ikx, iky, ikz, ppd, (double *) &ehat);
        // Set the sign of the z component (because the real FFT only gives the +kz half-space)
        ehat.vec[2] *= copysign(1, kz);
        // Linear interpolation might not preserve |ehat| = 1, so enforce this
        double ehatmag = sqrt(ehat.vec[0]*ehat.vec[0] + ehat.vec[1]*ehat.vec[1] + ehat.vec[2]*ehat.vec[2]);
        ehat.vec[0] /= ehatmag; ehat.vec[1] /= ehatmag; ehat.vec[2] /= ehatmag;
        
        // This upweights each mode by 1/(khat*ehat)
        double norm = k2/( kx*ehat.vec[0] + ky*ehat.vec[1] + kz*ehat.vec[2] );
        if(k2 == 0.0 || !std::isfinite(norm)) norm = 0.0;
        e.vec[0] = norm*ehat.vec[0];
        e.vec[1] = norm*ehat.vec[1];
        e.vec[2] = norm*ehat.vec[2];
        e.val = ehat.val;
    } else {
        e.vec[0] = kx;
        e.vec[1] = ky;
        e.vec[2] = kz;
        e.val = 1;
    }

    return e;
}

void LoadPlane(BlockArray& array, Parameters& param, PowerSpectrum& Pk, 
                int yblock, int yres, Complx *slab, Complx *slabHer) {
    // Note that this function is called from within a parallel for-loop over yres
    // STimer cpu, fft;
    // cpu.Start();

    Complx D,F,G,H;
    double f;
    Complx I(0.0,1.0);
    double k2;
    unsigned int a;
    int x,y,z, kx,ky,kz, xHer,yresHer,zHer;
    int64_t ppd = array.ppd;
    int64_t ppdhalf = array.ppdhalf;
    double ifundamental = 1.0/param.fundamental; // store the inverse
    double fundamental2 = param.fundamental*param.fundamental; // store the square
    double ik_cutoff = 1.0/param.k_cutoff;   // store the inverse
    int just_density = param.qdensity == 2;  // Don't generate displacements

    double target_f = (sqrt(1. + 24*param.f_cluster) - 1)/4.;
    double a_NL, a0;
    if(param.qPLTrescale){
        a_NL = 1./(1+param.PLT_target_z);
        a0 = 1./(1+param.z_initial);
    } else { a_NL = a0 = 1.0; }

    // How many RNG calls do we skip?  We'll fast-forward this amount each time we resume
    int64_t nskip = 0;

    double k2_cutoff = param.nyquist*param.nyquist/(param.k_cutoff*param.k_cutoff);
    int ver = param.version;

    y = yres+yblock*array.block;

    pcg64 checkpoint;
    if(ver == 2){
        checkpoint = v2rng[y];
    } 

    ky = y>ppdhalf?y-ppd:y;        // Nyquist wrapping
    yresHer = array.block-1-yres;         // Reflection
    for (z=0;z<ppd;z++) {
        // Just crossed the wrap, skip MAX_PPD-ppd rows
        if(z == ppdhalf+1 && ver == 2)
            nskip += (MAX_PPD - ppd)*MAX_PPD;
        kz = z>ppdhalf?z-ppd:z;        // Nyquist wrapping
        zHer = ppd-z; if (z==0) zHer=0;     // Reflection
        for (x=0;x<ppd;x++) {
            // Just cross the wrap, skip MAX_PPD-ppd particles
            if(x == ppdhalf+1 && ver == 2)
                nskip += MAX_PPD - ppd;
            kx = x>ppdhalf?x-ppd:x;        // Nyquist wrapping
            xHer = ppd-x; if (x==0) xHer=0;    // Reflection
            // We will pack two complex arrays
            k2 = (kx*kx+ky*ky+kz*kz)*fundamental2;
            
            // Force Nyquist elements to zero, being extra careful with rounding
            int kmax = (double)ppdhalf*ik_cutoff+.5;
            if ( (abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax)
                    // Force all elements with wavenumber above k_cutoff (nominally k_Nyquist) to zero
                    || (!param.CornerModes && k2>=k2_cutoff)
                    // Pick out one mode
                    || (param.qonemode && !(kx==param.one_mode[0] && ky==param.one_mode[1] && kz==param.one_mode[2])) ) {
                D=0.0;

                if (ver == 2)
                    nskip++;
            }
            else if (ver != 1){
                if(nskip){
                    v2rng[y].advance(2*nskip);
                    nskip = 0;
                }
                D = Pk.cgauss<2>(sqrt(k2), y);
            } else {
                // We deliberately only call cgauss() if we are inside the k_cutoff region
                // to get the same phase for a given k and cutoff region, no matter the ppd
                D = Pk.cgauss<1>(sqrt(k2),yres);
            }
            // D = 0.1;    // If we need a known level
            
            k2 *= ifundamental; // Get units of F,G,H right
            if (k2==0.0) k2 = 1.0;  // Avoid divide by zero
            // if (!(ky==5)) D=0.0;    // Pick out one plane
            double ik2 = 1.0/k2;
            
            // No-op this math if we aren't going to use it
            if(D != 0.){
                eigenmode e = get_eigenmode(kx, ky, kz, ppd, param.qPLT);

                double rescale = 1.;
                f = 1.0;
                if(param.qPLT){
                    // This is f_growth, the logarithmic derivative of the growth factor that scales the velocities
                    // The corrections are sourced from:
                    // 1) PLT growth rate
                    // 2) Addition of a smooth, non-clustering component to the background (<= NOT A PLT EFFECT)
                    // If PLT is turned on, we have to combine the effects here.  If not, we apply f_cluster during output.
                    f = (sqrt(1. + 24*e.val*param.f_cluster) - 1)*.25; // 1/4 instead of 1/6 because v = alpha*u/t0 = 3/2*H*alpha*u

                    if(param.qPLTrescale){
                        /// double a_NL = 1./(1+param.PLT_target_z);
                        /// double a0 = 1./(1+param.z_initial);
                        // First is continuum linear theory growth rate, possibly including f_smooth
                        // Second is PLT growth rate, also including f_smooth
                        /// double target_f = (sqrt(1. + 24*param.f_cluster) - 1)/4.;
                        // double plt_f = (sqrt(1. + 24*e.val*param.f_cluster) - 1)/4.;
                        double plt_f = f;
                        rescale = pow(a_NL/a0, target_f - plt_f);
                    }
                }
        
                F = rescale*I*e.vec[0]*ik2*D;
                G = rescale*I*e.vec[1]*ik2*D;
                H = rescale*I*e.vec[2]*ik2*D;
            } else {
                F = G = H = 0.0;
                f = 0.;
            }
            
            if(!just_density){
                // fprintf(stderr,"%d %d %d   %d %d %d   %f   %f %f\n",
                // x,y,z, kx,ky,kz, k2, real(D), imag(D));
                // H = F = D = 0.0;   // Test that the Hermitian aspects work
                // Now A = D+iF and B = G+iH.  
                // A is in array 0; B is in array 1
                AYZX(slab,0,yres,z,x) = D+I*F;
                AYZX(slab,1,yres,z,x) = G+I*H;
                if(param.qPLT){
                    AYZX(slab,2,yres,z,x) = Complx(0,0) + I*F*f;
                    AYZX(slab,3,yres,z,x) = G*f + I*H*f;
                }
                // And we need to store the complex conjugate
                // in the reflected entry.  We are reflecting
                // each element.  Note that we are storing one element
                // displaced in y; we will need to fix this when loading
                // for the y transform.  For now, we want the two block
                // boundaries to match.  This means that the conjugates 
                // for y=0 are being saved, which will be used below.
                AYZX(slabHer,0,yresHer,zHer,xHer) = conj(D)+I*conj(F);
                AYZX(slabHer,1,yresHer,zHer,xHer) = conj(G)+I*conj(H);
                if(param.qPLT){
                    AYZX(slabHer,2,yresHer,zHer,xHer) = 0. + I*conj(F*f);
                    AYZX(slabHer,3,yresHer,zHer,xHer) = conj(G*f) + I*conj(H*f);
                }
            } else {
                AYZX(slab,0,yres,z,x) = D;
                AYZX(slabHer,0,yresHer,zHer,xHer) = conj(D);
            }
        }
    } // End the x-z loops

    if(ver != 1){
        //printf("Made %lu calls (nskip at end %lu)\n", (uint64_t) (v2rng[y]-checkpoint), nskip);
        v2rng[y].advance(2*nskip);
        assert(v2rng[y]-checkpoint == 2*MAX_PPD*MAX_PPD);
    }

    // Need to do something special for ky=0 to enforce the 
    // Hermitian structure.  Recall that this whole plane was
    // stored in reflection and conjugate; we just need to copy
    // half of it back.
    if (yblock==0&&yres==0) {
        // Copy the first half plane onto the second
        for (z=0;z<ppdhalf;z++) {
            zHer = ppd-z; if (z==0) zHer=0;
            // Treat y=z=0 as a half line
            int xmax = (z==0?ppdhalf:ppd);
            for (x=0;x<xmax;x++) {
                xHer = ppd-x;if (x==0) xHer=0;
                for (a=0;a<array.narray;a++) {
                    AYZX(slab,a,yres,zHer,xHer) =
                    (AYZX(slabHer,a,yresHer,zHer,xHer));
                }
            }
        }
        // And the origin must be zero
        for (a=0;a<array.narray;a++) AYZX(slab,a,0,0,0) = 0.0; 
    }

    // Now do the Z FFTs, since those data are contiguous
    // cpu.Stop();
    // fft.Start();
    for (a=0;a<array.narray;a++) {
        InverseFFT_Yonly(&(AYZX(slab,a,yres,0,0)),ppd);
        InverseFFT_Yonly(&(AYZX(slabHer,a,yresHer,0,0)),ppd);
    }
    // fft.Stop();
    // printf("FFT fraction %f\n", fft.Elapsed()/(cpu.Elapsed()+fft.Elapsed()));
    return;
}

void ZeldovichZ(BlockArray& array, Parameters& param, PowerSpectrum& Pk) {
    // Generate the Fourier space density field, one Y block at a time
    // Use it to generate all arrays (density, qx, qy, qz) in Fourier space,
    // Do Z direction inverse FFTs.
    // Pack the result into 'array'.
    Complx *slab, *slabHer;
    int64_t len = (int64_t)array.block*array.ppd*array.ppd*array.narray;
    
    //slab    = new Complx[len];  // avoid hidden single-threaded zeroing
    //slabHer = new Complx[len];

    int ret = posix_memalign((void **) &slab, 4096, sizeof(Complx)*len);
    assert(ret == 0);
    ret = posix_memalign((void **) &slabHer, 4096, sizeof(Complx)*len);
    assert(ret == 0);
    #pragma omp parallel for schedule(static)
    for(int64_t i = 0; i < len; i++){
        slab[i] = Complx(0.,0.);
        slabHer[i] = Complx(0.,0.);
    }

    STimer compute_planes, storing;
    //
    fprintf(stderr,"Looping over Y: ");
    for (int yblock=0;yblock<array.numblock/2;yblock++) {
        // We're going to do each pair of Y slabs separately.
        // Load the deltas and do the FFTs for each pair of planes
        fprintf(stderr,".."); fflush(stderr);
        compute_planes.Start();
        #pragma omp parallel for schedule(dynamic,1)
        for (int yres=0;yres<array.block;yres++) {     
            LoadPlane(array,param,Pk,yblock,yres,slab,slabHer);
        }
        compute_planes.Stop();

        // Now store it into the primary BlockArray.  
        // Can't openMP an I/O loop.
        storing.Start();
        #pragma omp parallel for schedule(dynamic,1)
        for (int zblock=0;zblock<array.numblock;zblock++) {
            array.StoreBlock(yblock,zblock,slab);
            array.StoreBlock(array.numblock-1-yblock,zblock,slabHer);
        }
        storing.Stop();
    }  // End yblock for loop
    free(slabHer);
    free(slab);
    fprintf(stderr,"\n");
    fprintf(stderr,"Computing, Saving the Planes took %f %f sec\n", compute_planes.Elapsed(), storing.Elapsed());
    return;
}

// ===============================================================

// We use a set of X-Y arrays of Complx numbers (ordered by A and Z).
#define AZYX(_slab,_a,_z,_y,_x) _slab[(int64_t)(_x)+array.ppd*((_y)+array.ppd*((_a)+array.narray*(_z)))]


void ZeldovichXY(BlockArray& array, Parameters& param, FILE *output) {
    // Do the Y & X inverse FFT and output the results.
    // Do this one Z slab at a time; try to load the data in order.
    // Try to write the output file in z order
    
    Complx *slab;
    int64_t len = (int64_t)array.block*array.ppd*array.ppd*array.narray;
    //slab = new Complx[len];

    int ret = posix_memalign((void **) &slab, 4096, sizeof(Complx)*len);
    assert(ret == 0);
    #pragma omp parallel for schedule(static)
    for(int64_t i = 0; i < len; i++){
        slab[i] = Complx(0.,0.);
    }

    fprintf(stderr,"Looping over Z: ");
    STimer loading, writing, fft;

    for (int zblock=0;zblock<array.numblock;zblock++) {
        // We'll do one Z slab at a time
        // Load the slab back in.  
        // Can't openMP an I/O loop.
        loading.Start();
        fprintf(stderr,".");
        #pragma omp parallel for schedule(dynamic,1)
        for (int yblock=0;yblock<array.numblock;yblock++) {
            array.LoadBlock(yblock, zblock, slab);
        } 
        loading.Stop();

        // The Nyquist frequency y=array.ppd/2 must now be set to 0
        // because we shifted the data by one location.
        // FLAW: this assumes PPD is even.
        fft.Start();
        int y = array.ppd/2;
        #pragma omp parallel for schedule(static)
        for (int zres=0;zres<array.block;zres++) {
            for (int a=0;a<array.narray;a++) {
                for (int x=0;x<array.ppd;x++) AZYX(slab,a,zres,y,x) = 0.0;
            }
        }

        // Now we want to do the Y & X inverse FFT.
        for (int a=0;a<array.narray;a++) {
            #pragma omp parallel for schedule(dynamic,1)
            for (int zres=0;zres<array.block;zres++) {
                Inverse2dFFT(&(AZYX(slab,a,zres,0,0)),array.ppd);
            }
        }
        fft.Stop();

        // Now write out these rows of [z][y][x] positions
        // TODO: For now, we can't openMP an I/O loop.
        // There are internal buffers in the output routines that need to be made thread-safe.

        writing.Start();
        for (int zres=0;zres<array.block;zres++) {
            int z = zres+array.block*zblock;
            if (param.qoneslab<0||z==param.qoneslab) {
                // We have the option to output only one z slab.
                WriteParticlesSlab(output,z,
                &(AZYX(slab,0,zres,0,0)), &(AZYX(slab,1,zres,0,0)),
                &(AZYX(slab,2,zres,0,0)), &(AZYX(slab,3,zres,0,0)),
                array, param);
            }
        }
        writing.Stop();
    } // End zblock for loop
    free(slab);
    fprintf(stderr,"\n");
    fprintf(stderr,"Loading, FFTs, Writing took %f %f %f seconds\n", loading.Elapsed(), fft.Elapsed(), writing.Elapsed());
    return;
}

// ===============================================================

void load_eigmodes(Parameters &param){
    fprintf(stderr,"Using PLT eigenmodes.\n");
    // The eigvecs file consists of the ppd (32-bit int)
    // followed by PPDxPPDx(PPD/2+1)*4 doubles
    std::ifstream eigf;
    eigf.open(param.PLT_filename, std::ios::in|std::ios::binary|std::ios::ate);  //opens to end of file
    if (!eigf){
        std::cerr << "[Error] Could not open eigenmode file \"" << param.PLT_filename << "\".\n";
        exit(1);
    }
    std::streampos size;
    size = eigf.tellg();
    eigf.seekg (0, std::ios::beg);
    // Read as a 4-byte int, but store as a 8-byte int
    int eig_vecs_ppd_32;
    eigf.read((char*) &eig_vecs_ppd_32, sizeof(eig_vecs_ppd_32));
    eig_vecs_ppd = (int64_t) eig_vecs_ppd_32;
    
    size_t nelem = eig_vecs_ppd*eig_vecs_ppd*(eig_vecs_ppd/2 + 1)*4;
    size_t nbytes = nelem*sizeof(double);
    if((size_t) size != nbytes + sizeof(eig_vecs_ppd_32)){
        std::cerr << "[Error] Eigenmode file \"" << param.PLT_filename << "\" of size " << size
            << " did not match expected size " << nbytes << " from eig_vecs_ppd " << eig_vecs_ppd << ".\n";
        exit(1);
    }
    
    // eig_vecs = new double[nelem];
    int ret = posix_memalign((void **)&eig_vecs, 4096, sizeof(double)*nelem);
    assert(ret==0);
    eigf.read((char*) eig_vecs, nbytes);
    
    eigf.close();
}

// ===============================================================

// Determine which parts to emit
#if defined(PART1)
    #define NOPART2
    const int PART = 1;
    #ifdef PART2
    #error "Define at most one of PART1 and PART2"
    #endif
#elif defined(PART2)
    #define NOPART1
    const int PART = 2;
#else
    const int PART = -1;
#endif

int main(int argc, char *argv[]) {
    if (argc != 2){
        fprintf(stderr,"Usage: %s param_file\n", argv[0]);
        exit(1);
    }

    if(PART > 0){
        fprintf(stderr, "This is zeldovich part %d\n", PART);
    }

    STimer totaltime;
    totaltime.Start();
    
    FILE *output;
    double memory;
    density_variance = 0.0;
    Parameters param(argv[1]);

    PowerSpectrum Pk(10000);
    if(strlen(param.Pk_filename) > 0){
        if (Pk.InitFromFile(param.Pk_filename,param)!=0) return 1;
    } else {
        if (Pk.InitFromPowerLaw(param.Pk_powerlaw_index,param)!=0) return 1;
    }
    
    if(!Pk.is_powerlaw)
        param.append_file_to_comments(param.Pk_filename);

    //param.print(stdout);   // Inform the command line user
    // Two arrays for dens,x,y,z, two more for vx,vy,vz
    int narray;
    if(param.qdensity == 2){
        narray = 1;
    } else {
        narray = param.qPLT ? 4 : 2;
    }   
        
    memory = CUBE(param.ppd/1024.0)*narray*sizeof(Complx);

#ifndef NOPART1
    // Remove any existing IC files
    SetupOutputDir(param);
#endif

double _outbufferGiB = 0;
#ifndef NOPART2
    _outbufferGiB = InitOutputBuffers(param);
#endif

#ifdef DISK
    fprintf(stderr,"Compiled with -DDISK; blocks will be buffered on disk.\n");
    fprintf(stderr,"Total (out-of-core) memory usage: %5.3f GiB\n", memory);
    fprintf(stderr,"Two slab (in-core) memory usage: %5.3f GiB\n",
        memory/param.numblock*2.0 + memory/param.numblock/param.numblock + _outbufferGiB);  // extra from StoreBlock_tmp
    fprintf(stderr,"Block file size: %5.3f GiB\n", memory/param.numblock/param.numblock);
#else
    fprintf(stderr,"Not compiled with -DDISK; whole problem will reside in memory.\n");
    fprintf(stderr,"Total memory usage: %5.3f GiB\n", memory + memory/param.numblock*2.0 + _outbufferGiB);  // extra is from 2 blocks in ZeldovichZ
    fprintf(stderr,"Two slab memory usage: %5.3f GiB\n", memory/param.numblock*2.0);
    fprintf(stderr,"Block size: %5.3f GiB\n", memory/param.numblock/param.numblock);
#endif

    if(param.qPLT){
        load_eigmodes(param);
    }
    
    if(param.k_cutoff != 1){
        fprintf(stderr,"Using k_cutoff = %f (effective ppd = %d)\n", param.k_cutoff, (int)(param.ppd/param.k_cutoff+.5));
    }

    int quickdelete = 1;
    BlockArray array(param.ppd,param.numblock,narray,param.output_dir,param.AllowDirectIO, quickdelete, PART);

    totaltime.Stop();
    fprintf(stderr,"Preamble took %f seconds\n", totaltime.Elapsed());
    totaltime.Start();
    Setup_FFTW(param.ppd);

    totaltime.Stop();
    fprintf(stderr,"Time so far: %f seconds\n", totaltime.Elapsed());
    totaltime.Start();

#ifndef NOPART1
    ZeldovichZ(array, param, Pk);
    fprintf(stderr,"Wrote %d files\n", array.files_written);
    totaltime.Stop();
    fprintf(stderr,"Time so far: %f seconds\n", totaltime.Elapsed());
    totaltime.Start();
#endif

#ifndef NOPART2
    output = 0; // Current implementation doesn't use user-provided output
    ZeldovichXY(array, param, output);
    totaltime.Stop();
    fprintf(stderr,"Time so far: %f seconds\n", totaltime.Elapsed());
    totaltime.Start();

    fprintf(stderr,"The rms density variation of the pixels is %f\n", sqrt(density_variance/CUBE(param.ppd)));
    fprintf(stderr,"This could be compared to the P(k) prediction of %f\n",
    Pk.sigmaR(param.separation/4.0)*pow(param.boxsize,1.5));
    
    if(param.qdensity != 2){
        fprintf(stderr,"The maximum component-wise displacements are (%g, %g, %g), same units as BoxSize.\n", max_disp[0], max_disp[1], max_disp[2]);
        fprintf(stderr,"For Abacus' 2LPT implementation to work (assuming FINISH_WAIT_RADIUS = 1),\n\tthis implies a maximum CPD of %d\n", (int) (param.boxsize/(2*max_disp[2])));  // The slab direction is z in this code
    }
    // fclose(output);
    totaltime.Stop();
    fprintf(stderr,"Time so far: %f seconds\n", totaltime.Elapsed());
    totaltime.Start();

    TeardownOutput();
#endif
    
    if(param.qPLT)
        free(eig_vecs);
        // delete[] eig_vecs;

    totaltime.Stop();
    fprintf(stderr, "zeldovich took %.4g sec for ppd %lu ==> %.3g Mpart/sec\n",
        totaltime.Elapsed(), param.ppd, param.np/1e6/totaltime.Elapsed());
    
    return 0;
}
