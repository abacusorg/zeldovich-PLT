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

#include "../include/STimer.cc"

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

    // For big n, this is slow enough to notice
    if(n >= 512)
        fprintf(stderr,"Creating FFTW plans...");

    fftw_complex *p;
    p = new fftw_complex[n*n];
    plan1d = fftw_plan_dft_1d(n, p, p, +1, FFTW_PATIENT);
    plan2d = fftw_plan_dft_2d(n, n, p, p, +1, FFTW_PATIENT);
    delete []p;

    if(n >= 512)
        fprintf(stderr," done.\n");
}

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
    tmp = new Complx[n];
    for (j=0;j<n;j++) {
        // We will load one row at a time
        for (k=0;k<n;k++) tmp[k] = p[k*n+j];
        Inverse1dFFT(tmp, n);
        for (k=0;k<n;k++) p[k*n+j] = tmp[k];
    }
    delete []tmp;
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
    if(fx > eig_vecs_ppd/2 && fx < eig_vecs_ppd/2 + 1)
        fx = floor(fx+1);
    if(fy > eig_vecs_ppd/2 && fy < eig_vecs_ppd/2 + 1)
        fy = floor(fy+1);
    if(fz > eig_vecs_ppd/2 && fz < eig_vecs_ppd/2 + 1)
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

    Complx D,F,G,H,f;
    Complx I(0.0,1.0);
    double k2;
    unsigned int a;
    int x,y,z, kx,ky,kz, xHer,yresHer,zHer;
    int64_t ppd = array.ppd;

    // How many RNG calls do we skip?  We'll fast-forward this amount each time we resume
    int64_t nskip = 0;

    double k2_cutoff = param.nyquist*param.nyquist/(param.k_cutoff*param.k_cutoff);
    int ver = param.version;

    y = yres+yblock*array.block;

    pcg64 checkpoint;
    if(ver == 2){
        checkpoint = v2rng[y];
    } 

    ky = y>ppd/2?y-ppd:y;        // Nyquist wrapping
    yresHer = array.block-1-yres;         // Reflection
    for (z=0;z<ppd;z++) {
        // Just crossed the wrap, skip MAX_PPD-ppd rows
        if(z == ppd/2+1 && ver == 2)
            nskip += (MAX_PPD - ppd)*MAX_PPD;
        kz = z>ppd/2?z-ppd:z;        // Nyquist wrapping
        zHer = ppd-z; if (z==0) zHer=0;     // Reflection
        for (x=0;x<ppd;x++) {
            // Just cross the wrap, skip MAX_PPD-ppd particles
            if(x == ppd/2+1 && ver == 2)
                nskip += MAX_PPD - ppd;
            kx = x>ppd/2?x-ppd:x;        // Nyquist wrapping
            xHer = ppd-x; if (x==0) xHer=0;    // Reflection
            // We will pack two complex arrays
            k2 = (kx*kx+ky*ky+kz*kz)*param.fundamental*param.fundamental;
            
            // Force Nyquist elements to zero, being extra careful with rounding
            int kmax = ppd/2./param.k_cutoff+.5;
            if ( (abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax)
                    // Force all elements with wavenumber above k_cutoff (nominally k_Nyquist) to zero
                    || (k2>=k2_cutoff)
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
            
            k2 /= param.fundamental; // Get units of F,G,H right
            if (k2==0.0) k2 = 1.0;  // Avoid divide by zero
            // if (!(ky==5)) D=0.0;    // Pick out one plane
            
            // No-op this math if we aren't going to use it
            if(D != 0.){
                eigenmode e = get_eigenmode(kx, ky, kz, ppd, param.qPLT);
                double rescale = 1.;
                if(param.qPLTrescale){
                    double a_NL = 1./(1+param.PLT_target_z);
                    double a0 = 1./(1+param.z_initial);
                    double alpha_m = (sqrt(1. + 24*e.val) - 1)/6.;
                    rescale = pow(a_NL/a0, 1 - 1.5*alpha_m);
                }
                F = rescale*I*e.vec[0]/k2*D;
                G = rescale*I*e.vec[1]/k2*D;
                H = rescale*I*e.vec[2]/k2*D;
                
                if(param.qPLT)
                    f = (sqrt(1. + 24*e.val) - 1)*.25; // 1/4 instead of 1/6 because v = alpha*u/t0 = 3/2*H*alpha*u
            } else {
                F = G = H = f = 0.;
            }
            
            // printf("%d %d %d   %d %d %d   %f   %f %f\n",
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
        for (z=0;z<ppd/2;z++) {
            zHer = ppd-z; if (z==0) zHer=0;
            // Treat y=z=0 as a half line
            int xmax = (z==0?ppd/2:ppd);
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
    for (a=0;a<array.narray;a++) {
        InverseFFT_Yonly(&(AYZX(slab,a,yres,0,0)),ppd);
        InverseFFT_Yonly(&(AYZX(slabHer,a,yresHer,0,0)),ppd);
    }
    return;
}

void StoreBlock(BlockArray& array, int yblock, int zblock, Complx *slab) {
    // We must be sure to store the block sequentially.
    // data[zblock=0..NB-1][yblock=0..NB-1]
    //     [array=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
    // Can't openMP an I/O loop.
    unsigned int a;
    int yres,zres,z;
    array.bopen(yblock,zblock,"w");
    for (a=0;a<array.narray;a++) 
    for (zres=0;zres<array.block;zres++) 
    for (yres=0;yres<array.block;yres++) {
        z = zres+array.block*zblock;
        //y = yres+array.block*yblock;  // the slab is in the y coordinate, so we never need the absolute y index
        // Copy the whole X skewer
        array.bwrite(&(AYZX(slab,a,yres,z,0)),array.ppd);
    }
    array.bclose();
    return;
}

void ZeldovichZ(BlockArray& array, Parameters& param, PowerSpectrum& Pk) {
    // Generate the Fourier space density field, one Y block at a time
    // Use it to generate all arrays (density, qx, qy, qz) in Fourier space,
    // Do Z direction inverse FFTs.
    // Pack the result into 'array'.
    Complx *slab, *slabHer;
    int yblock,zblock;
    int64_t len = (int64_t)array.block*array.ppd*array.ppd*array.narray;
    slab    = new Complx[len];
    slabHer = new Complx[len];
    //
    printf("Looping over Y: ");
    for (yblock=0;yblock<array.numblock/2;yblock++) {
        // We're going to do each pair of Y slabs separately.
        // Load the deltas and do the FFTs for each pair of planes
        printf(".."); fflush(stdout);
        #pragma omp parallel for schedule(static)
        for (int yres=0;yres<array.block;yres++) {     
            LoadPlane(array,param,Pk,yblock,yres,slab,slabHer);
        }

        // Now store it into the primary BlockArray.  
        // Can't openMP an I/O loop.
        for (zblock=0;zblock<array.numblock;zblock++) {
            StoreBlock(array,yblock,zblock,slab);
            StoreBlock(array,array.numblock-1-yblock,zblock,slabHer);
        }
    }  // End yblock for loop
    delete []slabHer;
    delete []slab;
    printf("\n"); fflush(stdout);
    return;
}

// ===============================================================

// We use a set of X-Y arrays of Complx numbers (ordered by A and Z).
#define AZYX(_slab,_a,_z,_y,_x) _slab[(int64_t)(_x)+array.ppd*((_y)+array.ppd*((_a)+array.narray*(_z)))]

void LoadBlock(BlockArray& array, int yblock, int zblock, Complx *slab) {
    // We must be sure to access the block sequentially.
    // data[zblock=0..NB-1][yblock=0..NB-1]
    //     [array=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
    // Can't openMP an I/O loop.
    unsigned int a;
    int yres,y,zres,yshift;
    array.bopen(yblock,zblock,"r");
    for (a=0;a<array.narray;a++)
    for (zres=0;zres<array.block;zres++) 
    for (yres=0;yres<array.block;yres++) {
        //z = zres+array.block*zblock;  // slabs are in z; never need the absolute z coord
        y = yres+array.block*yblock;
        // Copy the whole X skewer.  However, we want to
        // shift the y frequencies in the reflected half
        // by one.
        // FLAW: Assumes array.ppd is even.
        if (y>=array.ppd/2) yshift=y+1; else yshift=y;
        if (yshift==array.ppd) yshift=array.ppd/2;
        // Put it somewhere; this is about to be overwritten
        array.bread(&(AZYX(slab,a,zres,yshift,0)),array.ppd);
    }
    array.bclose();
    return;
}

void ZeldovichXY(BlockArray& array, Parameters& param, FILE *output, FILE *densoutput) {
    // Do the Y & X inverse FFT and output the results.
    // Do this one Z slab at a time; try to load the data in order.
    // Try to write the output file in z order
    
    Complx *slab;
    int64_t len = (int64_t)array.block*array.ppd*array.ppd*array.narray;
    slab = new Complx[len];
    unsigned int a;
    int x,yblock,y,zblock,z;

    printf("Looping over Z: ");
    for (zblock=0;zblock<array.numblock;zblock++) {
        // We'll do one Z slab at a time
        // Load the slab back in.  
        // Can't openMP an I/O loop.
        printf("."); fflush(stdout);
        for (yblock=0;yblock<array.numblock;yblock++) {
            LoadBlock(array, yblock, zblock, slab);
        } 

        // The Nyquist frequency y=array.ppd/2 must now be set to 0
        // because we shifted the data by one location.
        // FLAW: this assumes PPD is even.
        y = array.ppd/2;
        for (int zres=0;zres<array.block;zres++) {
            for (a=0;a<array.narray;a++) {
                for (x=0;x<array.ppd;x++) AZYX(slab,a,zres,y,x) = 0.0;
            }
        }

        // Now we want to do the Y & X inverse FFT.
        for (a=0;a<array.narray;a++) {
            #pragma omp parallel for schedule(static)
            for (int zres=0;zres<array.block;zres++) {
                Inverse2dFFT(&(AZYX(slab,a,zres,0,0)),array.ppd);
            }
        }

        // Now write out these rows of [z][y][x] positions
        // Can't openMP an I/O loop.
        

        for (int zres=0;zres<array.block;zres++) {
            z = zres+array.block*zblock;
            if (param.qoneslab<0||z==param.qoneslab) {
                // We have the option to output only one z slab.



                WriteParticlesSlab(output,densoutput,z,
                &(AZYX(slab,0,zres,0,0)), &(AZYX(slab,1,zres,0,0)),
                &(AZYX(slab,2,zres,0,0)), &(AZYX(slab,3,zres,0,0)),
                array, param);
            }
        }
    } // End zblock for loop
    delete []slab;
    printf("\n"); fflush(stdout);
    return;
}

// ===============================================================

void load_eigmodes(Parameters &param){
    printf("Using PLT eigenmodes.\n");
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
    
    eig_vecs = new double[nelem];
    eigf.read((char*) eig_vecs, nbytes);
    
    eigf.close();
}

int main(int argc, char *argv[]) {
    if (argc != 2){
        printf("Usage: %s param_file\n", argv[0]);
        exit(1);
    }
    
    FILE *output, *densoutput;
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
    int narray = param.qPLT ? 4 : 2;
    memory = CUBE(param.ppd/1024.0)*narray*sizeof(Complx);

#ifdef DISK
    printf("Compiled with -DDISK; blocks will be buffered on disk.\n");
    printf("Total (out-of-core) memory usage (GB): %5.3f\n", memory);
    printf("Two slab (in-core) memory usage (GB): %5.3f\n", memory/param.numblock*2.0);
    printf("Block file size (GB): %5.3f\n", memory/param.numblock/param.numblock);
#else
    printf("Not compiled with -DDISK; whole problem will reside in memory.\n");
    printf("Total memory usage (GB): %5.3f\n", memory + memory/param.numblock*2.0);
    printf("Two slab memory usage (GB): %5.3f\n", memory/param.numblock*2.0);
    printf("Block size (GB): %5.3f\n", memory/param.numblock/param.numblock);
#endif

    if (param.qdensity>0) {
        densoutput = fopen(param.density_filename,"w");
        assert(densoutput!=NULL);
    } else densoutput = NULL;

    if(param.qPLT){
        load_eigmodes(param);
    }
    
    if(param.k_cutoff != 1){
        printf("Using k_cutoff = %f (effective ppd = %d)\n", param.k_cutoff, (int)(param.ppd/param.k_cutoff+.5));
    }

    Setup_FFTW(param.ppd);

    // We're about to start the BlockArray, thus using the disk
    // Remove any existing IC files
    CleanDirectory(param.output_dir);
    CreateDirectories(param.output_dir);

    BlockArray array(param.ppd,param.numblock,narray,param.output_dir,param.ramdisk);    
    ZeldovichZ(array, param, Pk);

    output = 0; // Current implementation doesn't use user-provided output
    ZeldovichXY(array, param, output, densoutput);

    printf("The rms density variation of the pixels is %f\n", sqrt(density_variance/CUBE(param.ppd)));
    printf("This could be compared to the P(k) prediction of %f\n",
    Pk.sigmaR(param.separation/4.0)*pow(param.boxsize,1.5));
    
    printf("The maximum component-wise displacements are (%g, %g, %g).\n", max_disp[0], max_disp[1], max_disp[2]);
    printf("For Abacus' 2LPT implementation to work (assuming FINISH_WAIT_RADIUS = 1),\nthis implies a maximum CPD of %d\n", (int) (param.boxsize/(2*max_disp[2])));  // The slab direction is z in this code
    // fclose(output);
    
    if(param.qPLT)
        delete[] eig_vecs;
    
    return 0;
}
