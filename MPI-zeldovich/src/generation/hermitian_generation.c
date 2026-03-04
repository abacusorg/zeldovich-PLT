#include "hermitian_generation.h"
#include "../utils/verification.h"
#include "../utils/zeldovich_wrapper.h" 
#include "../utils/plt_eigenmodes.h" 
#include "../config.h"  
#include "../precision.h"  
#include <stdio.h>
#include <stdlib.h>  
#include <stdint.h> 
#include <math.h>    
#include <omp.h>

// ====================================================================================
// Thread-local RNG helpers for parallel z-loop
// ====================================================================================

#if PARALLELIZE_Z_LOOP
// Map (z, x) to virtual position in MAX_PPD x MAX_PPD grid (stream position in rand counts)
static inline int64_t compute_virtual_position(int z, int x, int N, int Nhalf) {
    int z_v = (z <= Nhalf) ? z : (MAX_PPD - N + z);
    int x_v = (x <= Nhalf) ? x : (MAX_PPD - N + x);
    return 2 * ((int64_t)z_v * MAX_PPD + x_v); // divide by 2 later
}

// manual static chunking:
// split z = 0..N-1 range into ~equal chunks for each thread 
// last thread may get slightly shorter chunk if N not divisible by nthreads.
static inline void get_thread_z_range(int tid, int nthreads, int N, int* z_start, int* z_end) {
    int chunk_size = (N + nthreads - 1) / nthreads;
    *z_start = tid * chunk_size;
    *z_end = (*z_start + chunk_size > N) ? N : *z_start + chunk_size;
}
#endif

// Choose which RNG to use based on local_rng_buf
// local_rng_buf == NULL: use shared RNG (zeldovich_ps_advance_rng, zeldovich_ps_cgauss)
// local_rng_buf != NULL: use buffer RNG (zeldovich_ps_advance_rng_buffer, zeldovich_ps_cgauss_from_buffer)
static inline void get_cgauss(PowerSpectrumHandle ps_handle, ParametersHandle params_handle,
    int64_t rng_index, double kmag, int64_t* nskip, void* local_rng_buf,
    double* D_real, double* D_imag) {
    if (*nskip > 0) {
        if (local_rng_buf)
            zeldovich_ps_advance_rng_buffer(local_rng_buf, *nskip);
        else
            zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, *nskip);
        *nskip = 0;
    }
    if (local_rng_buf)
        zeldovich_ps_cgauss_from_buffer(local_rng_buf, ps_handle, kmag, D_real, D_imag);
    else
        zeldovich_ps_cgauss(ps_handle, kmag, rng_index, D_real, D_imag);
}

// ====================================================================================
// Generates one pair of Hermitian Y-slices (primary + conjugate) in Fourier space
// Gaussian w/ power spectrum weighting (w/ zeldovich-PLT ps_handle & params_handle)
// 2D FFT: Fourier --> real space (X,Z)
// Handles RNG nskip tracking
// ====================================================================================

void generate_hermitian_slice_pair_local(
    int N,
    int global_y,
    int y_mirror,
    fftw_complex_t *primary_slices,   
    fftw_complex_t *conjugate_slices, 
    int narray,                       // 4 arrays per slice, 7 C numbers
    fftw_plan_t plan_2d,
    int rank,
    PowerSpectrumHandle ps_handle,   // zeldovich-PLT PowerSpectrum handle
    ParametersHandle params_handle,  // zeldovich-PLT Parameters handle
    void** thread_rng_buffers)        // Pre-allocated RNG buffers [nthreads] (NULL = use malloc)
{
    // ========== DIAGNOSTIC TIMING: Function-level ==========
    double t_func_start = omp_get_wtime();
    double t_setup_end, t_zloop_end, t_verify_end, t_fft_end;
    
    // Debug: Log entry for seg fault error
    // #if DEBUG_PRINTS
    // fprintf(stderr,
    //         "[Rank %d] ENTER generate_hermitian_slice_pair_local: "
    //         "Y_primary=%d, Y_mirror=%d, ps_handle=%p, params_handle=%p\n",
    //         rank, global_y, y_mirror, (void*)ps_handle, (void*)params_handle);
    // fflush(stderr);
    // #endif

    // Precompute N/2 to avoid repeated division in loops
    int Nhalf = N / 2;
    
    // ========== Check for density-only mode (qdensity == 2) ==========
    // In density-only mode, skip compute F, G, H (displacements), only store D (density)
    int just_density = 0;
    if (params_handle != NULL) {
        int qdensity = zeldovich_params_get_qdensity(params_handle);
        just_density = (qdensity == 2);
    }
    
    // ========== PLT Rescaling factors (computed once per function call) ==========
    // These are used to rescale F, G, H when qPLTrescale is enabled
    // target_f: continuum linear theory growth rate
    // a_NL: scale factor at target redshift
    // a0: scale factor at initial redshift
    double target_f = 1.0;
    double a_NL = 1.0;
    double a0 = 1.0;
    int qPLTrescale = 0;
    
    // Only compute rescaling factors if we're computing displacements
    if (params_handle != NULL && !just_density) {
        double f_cluster = zeldovich_params_get_f_cluster(params_handle);
        // target_f is the continuum linear theory growth rate (fluid limit)
        target_f = (sqrt(1. + 24. * f_cluster) - 1.) * 0.25;
        
        qPLTrescale = zeldovich_params_get_qPLTrescale(params_handle);
        if (qPLTrescale) {
            double z_initial = zeldovich_params_get_z_initial(params_handle);
            double PLT_target_z = zeldovich_params_get_PLT_target_z(params_handle);
            a_NL = 1. / (1. + PLT_target_z);  // Scale factor at target redshift
            a0 = 1. / (1. + z_initial);       // Scale factor at initial redshift
        }
    }

    // ========== k_cutoff filtering parameters (computed once per function call) ==========
    // Calculate k2_cutoff for filtering high-wavenumber modes
    // See zeldovich.cpp: k2_cutoff = nyquist^2 / (k_cutoff^2)
    double k_cutoff = 1.0;  // Default value
    int CornerModes = 0;    // Default value
    double k2_cutoff = 0.0; // Will be calculated
    
    if (params_handle != NULL) {
        k_cutoff = zeldovich_params_get_k_cutoff(params_handle);
        CornerModes = zeldovich_params_get_CornerModes(params_handle);
    }
    
    // Calculate k2_cutoff: k2_cutoff = (N/2)^2 / (k_cutoff^2)
    // Ex. N=16, k_cutoff=1.0: k2_cutoff = 8^2 / 1.0^2 = 64
    double Nhalf_dbl = (double)Nhalf;
    k2_cutoff = (Nhalf_dbl * Nhalf_dbl) / (k_cutoff * k_cutoff);

    // Fundamental wavenumber (same for conjugate-pair and self-conjugate; set once per function call)
    double fundamental = 1.0;
    if (params_handle != NULL) {
        fundamental = zeldovich_params_get_fundamental(params_handle);
    }

    // Access primary and conj slices
    #define PRIM_SLICE(array_idx, x, z) \
        primary_slices[(int64_t)(x) + (N) * ((z) + (N) * (array_idx))]
     
    #define CONJ_SLICE(array_idx, x, z) \
        conjugate_slices[(int64_t)(x) + (N) * ((z) + (N) * (array_idx))]
    
    // RNG verification: Track total RNG calls and skips
    #if VERIFY_RNG_CALLS
    int64_t total_rng_calls = 0;  // Count of cgauss() calls (each uses 2 random numbers)
    int64_t total_rng_skips = 0;  // Count of skipped numbers (from nskip advances, including D=0 skips)
    #endif
    
    // RNG consistency: Skip random numbers for missing grid points when N < MAX_PPD
    // so we can reproduce sim boxes w/ up to MAX_PPD x MAX_PPD grid generators
    int64_t nskip = 0;
    if (N < MAX_PPD) {
        // Sequential: track skip incrementally during iteration (like zeldovich.cpp)
        nskip = 0;  // Will be accumulated during loops
    }
    
    t_setup_end = omp_get_wtime();
    
    if (y_mirror != global_y) {
        // ========== CONJUGATE PAIR: Y=i and Y=N-i ==========
#if PARALLELIZE_Z_LOOP
        #pragma omp parallel
        {
            int64_t rng_index_cp = global_y;
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
            int z_start, z_end;
            get_thread_z_range(tid, nthreads, N, &z_start, &z_end);
            
            // Diagnostic timing
            double t_start = omp_get_wtime();
            
            // Use persistent buffer if provided, otherwise allocate
            void* local_rng_buf;
            int need_free = 0;
            if (thread_rng_buffers != NULL && tid < omp_get_max_threads()) {
                local_rng_buf = thread_rng_buffers[tid];
            } else {
                fprintf(stderr, "[Rank %d] WARNING: No persistent RNG buffer for thread %d, allocating new one (EXPENSIVE!) \n", rank, tid);
                size_t rng_size = zeldovich_ps_rng_buffer_size();
                local_rng_buf = malloc(rng_size);
                need_free = 1;
            }
            
            double t_after_alloc = omp_get_wtime();
            
            zeldovich_ps_get_rng_copy(ps_handle, rng_index_cp, local_rng_buf);
            double t_after_copy = omp_get_wtime();
            
            int64_t virtual_start = compute_virtual_position(z_start, 0, N, Nhalf) / 2;
            if (virtual_start > 0)
                zeldovich_ps_advance_rng_buffer(local_rng_buf, virtual_start);
            double t_after_advance = omp_get_wtime();
            
            int64_t nskip = 0;
            for (int z = z_start; z < z_end; z++) {
#else
        void* local_rng_buf = NULL;
        for (int z = 0; z < N; z++) {
#endif
            // RNG consistency: When crossing Nyquist boundary (z == Nhalf + 1),
            // skip ALL missing z-rows (z = N to MAX_PPD-1, each containing MAX_PPD x-values)
            // See zeldovich.cpp: skip at Nyquist boundary before processing negative kz region
            // The missing frequencies are in the MIDDLE of the MAX_PPD array (high positive and negative k),
            // due to FFT ordering: [0, 1, ..., N/2, -N/2+1, ..., -1]
            // High frequencies (both positive and negative) are located in the middle,
            // so we skip them all at once when we first enter the mirrored region
            if (z == Nhalf + 1 && N < MAX_PPD) {
                // Skip ALL missing z-rows: from z=N to z=MAX_PPD-1
                // This accounts for the "gap" in frequency space (high frequencies in the middle of array)
                int64_t skip_amount = (int64_t)(MAX_PPD - N) * (int64_t)MAX_PPD;
                nskip += skip_amount; // skip missing z-rows
                
                #if DEBUG_RNG_SKIP
                // Debug: Log skip accumulation for test coordinates
                int log_skip = (global_y <= MAX_DEBUG_COORD) || (global_y == Nhalf - 1);
                if (log_skip) {
                    fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d z=%d: ACCUMULATE skip at Nyquist boundary: +%lld (missing z-rows), total nskip=%lld\n",
                            N, global_y, z, (long long)skip_amount, (long long)nskip);
                    fflush(stderr);
                }
                #endif
            }
            
            // ky, kz, abs_ky, abs_kz are constant over x for this (z, global_y); compute once per z-row
            int ky = (global_y > Nhalf) ? global_y - N : global_y;
            int kz = (z > Nhalf) ? z - N : z;
            int abs_ky = (ky < 0) ? -ky : ky;
            int abs_kz = (kz < 0) ? -kz : kz;
            
            for (int x = 0; x < N; x++) {
                int x_mirror = (x == 0) ? 0 : N - x;
                int z_mirror = (z == 0) ? 0 : N - z;
                
                // ========== STEP 1: Calculate k-vector components ==========
                // RNG consistency: When crossing Nyquist boundary (x == Nhalf + 1),
                if (x == Nhalf + 1 && N < MAX_PPD) {
                    // Skip ALL missing x-values: from x=N to x=MAX_PPD-1
                    int64_t skip_amount = MAX_PPD - N;
                    nskip += skip_amount; // skip missing x-values in this z-row
                    
                    #if DEBUG_RNG_SKIP
                    // Debug: Log skip accumulation for test coordinates
                    int log_skip = (z <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD) ||
                                   (z == Nhalf - 1 && global_y == Nhalf - 1);
                    if (log_skip) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d z=%d x=%d: ACCUMULATE skip at Nyquist boundary: +%lld (missing x-values), total nskip=%lld\n",
                                N, global_y, z, x, (long long)skip_amount, (long long)nskip);
                        fflush(stderr);
                    }
                    #endif
                }
                
                // kx depends on x; ky, kz, abs_ky, abs_kz already set above
                int kx = (x > Nhalf) ? x - N : x;
                int k2_int = kx*kx + ky*ky + kz*kz;  // int for k_cutoff comparison
                double k2 = (double)k2_int;  // float for p(k)
                
                // ========== STEP 2: Generate D using RNG or cgauss() ==========
                // Nyquist frequency zeroing: Force Nyquist elements to zero for all three axes
                // This matches zeldovich.cpp line 354: abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax
                // The Nyquist frequency (k = N/2) is self-conjugate and doesn't have a separate partner.
                // Due to the Y-shift in the reflected shell, we need to zero these to align mirroring expectations.
                int abs_kx = (kx < 0) ? -kx : kx;
                int is_nyquist = (abs_kx == Nhalf || abs_ky == Nhalf || abs_kz == Nhalf);
                
                fftw_complex D;
                // Zero D for: DC mode, Nyquist frequency, or k_cutoff filtering
                if ((k2 == 0.0)
                     || (is_nyquist)
                    // Force all elements with wavenumber above k_cutoff (nominally k_Nyquist) to zero
                     || (!CornerModes && (double)k2_int >= k2_cutoff)) {
                    D[0] = D[1] = 0.0;
                    // RNG consistency: When D=0, we skip the RNG call - accumulate nskip, dont advance immediately
                    nskip++;  // Accumulate skip, will be applied before next cgauss() call
                    #if DEBUG_RNG_SKIP
                    int log_skip = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                   (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_skip) {
                        if (k2 == 0.0) {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (DC mode), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, (long long)nskip);
                        } else if (is_nyquist) {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (Nyquist: kx=%d ky=%d kz=%d), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, kx, ky, kz, (long long)nskip);
                        } else {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (k_cutoff: k2=%d >= %.1f), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, k2_int, k2_cutoff, (long long)nskip);
                        }
                        fflush(stderr);
                    }
                    #endif
                } 
                else if (ps_handle != NULL && params_handle != NULL) {
                    // Convert k indices to physical wavenumber
                    double k2_phys = k2 * fundamental * fundamental;
                    double kmag = sqrt(k2_phys);
                    
                    // zeldovich_ps_cgauss returns double precision, convert to real_t
                    // v2rng array is sized to ppd/2, so valid indices are 0 to (N/2 - 1)
                    // In conjugate pair branch, global_y is in range [1, N/2-1] or [N/2+1, N-1]
                    int64_t rng_index = global_y;
                    
                    double D_real, D_imag;
                    #if VERIFY_RNG_CALLS
                    if (nskip > 0 && local_rng_buf == NULL) total_rng_skips += nskip;
                    #endif
                    // Phase 4: use get_cgauss (shared RNG when local_rng_buf==NULL, buffer RNG when parallel)
                    get_cgauss(ps_handle, params_handle, rng_index, kmag, &nskip, local_rng_buf, &D_real, &D_imag);

                    #if VERIFY_RNG_CALLS
                    total_rng_calls++;  // Each cgauss() call uses 2 random numbers
                    #endif
                    D[0] = (real_t)D_real;
                    D[1] = (real_t)D_imag;
                    // Note: For VERIFY_HERMITIAN_SYMMETRY==2, we set D=0 AFTER computing F and H
                }
                
                // ========== STEP 3: Compute F, G, H from D ==========
                fftw_complex F, G, H;
                int use_plt = 0; 
                eigenmode e;      
                double f = 1.0; // for velocity arrays in all modes

                if (D[0] == 0.0 && D[1] == 0.0) {
                    f = 0.0;
                }
                
                if (just_density) {
                    // Density-only mode (qdensity == 2): Zero F, G, H
                    F[0] = F[1] = G[0] = G[1] = H[0] = H[1] = 0.0;
                } else {
                double ik2 = 1.0 / k2; // later: define this at top of code and use it instead of dividing
                double rescale = 1.0;
                double factor;

                #if defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 1
                // Verification mode 1: Set F=0 and H=0 to test Hermitian symmetry
                // With F=0 and H=0, conj slices should be true conjugates of primary slices (result = pure real)
                factor = rescale / (k2 * fundamental);
                
                F[0] = 0.0;
                F[1] = 0.0;
                
                // Compute G using the same formula as normal mode
                G[0] = -ky * factor * D[1];
                G[1] =  ky * factor * D[0];
                
                H[0] = 0.0;
                H[1] = 0.0;
                #elif defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 2
                // Verification mode 2: Compute F and H from D, then set D=0 and G=0
                // Result after 3D FFT: Array 0 and Array 1 are purely imag
                factor = rescale / (k2 * fundamental);
                
                // Now compute F and H using the same formula as normal mode
                F[0] = -kx * factor * D[1];
                F[1] =  kx * factor * D[0];
                
                G[0] = 0.0;  // Set G=0
                G[1] = 0.0;
                
                H[0] = -kz * factor * D[1];
                H[1] =  kz * factor * D[0];
                
                // Now set D=0 (after F and H are computed)
                D[0] = 0.0;
                D[1] = 0.0;

                #else
                // Normal operation: Compute F, G, H from D
                // Check if PLT is enabled (fundamental already set at top of else block)
                int qPLT = 0;
                if (params_handle != NULL) {
                    qPLT = zeldovich_params_get_qPLT(params_handle);
                    if (qPLT) {
                        // Get PLT eigenmode for this k-vector
                        // Convert kx, ky, kz to array indices for plt_get_eigenmode
                        int ikx = (kx < 0) ? N + kx : kx;
                        int iky = (ky < 0) ? N + ky : ky;
                        int ikz = (kz < 0) ? N + kz : kz;
                        // Handle z index (only positive half-space stored)
                        if (ikz > N / 2) ikz = N - ikz;
                        
                        if (plt_get_eigenmode(ikx, iky, ikz, (int64_t)N, &e) == 0) {
                            use_plt = 1;
                            double e_mag, k_dot_e, norm;

                            // Apply the same normalization as zeldovich.cpp
                            // 1. Set sign of z component (real FFT only gives +kz half-space)
                            // 2. Ensure |e| = 1 (normalize after interpolation, like zeldovich)
                            // 3. Apply norm = k2 / (k * e) to get the final eigenvector
                            
                            // Set the sign of the z component (because the real FFT only gives the +kz half-space)
                            // see zeldovich.cpp line 255: ehat.vec[2] *= copysign(1, kz);
                            if (kz < 0) {
                                e.vec[2] = -e.vec[2];
                            }
                            
                            // Normalize eigenvector (interpolation might not preserve |e| = 1)
                            // see zeldovich.cpp lines 257-263
                            e_mag = sqrt(e.vec[0] * e.vec[0] + e.vec[1] * e.vec[1] + e.vec[2] * e.vec[2]);
                            if (e_mag > 0.0) {
                                e.vec[0] /= e_mag;
                                e.vec[1] /= e_mag;
                                e.vec[2] /= e_mag;
                            }
                            
                            // Apply normalization: norm = k2 / (k * e)
                            // This upweights each mode by 1/(khat*ehat),see zeldovich.cpp line 266
                            k_dot_e = kx * e.vec[0] + ky * e.vec[1] + kz * e.vec[2];
                            norm = (k2 > 0.0 && k_dot_e != 0.0) ? k2 / k_dot_e : 0.0;
                            if (!isfinite(norm)) norm = 0.0;
                            
                            // Scale eigenvector by norm (see zeldovich.cpp lines 268-270)
                            e.vec[0] *= norm;
                            e.vec[1] *= norm;
                            e.vec[2] *= norm;
                        } else {
                            // If eigenmode lookup fails, fall back to normal computation
                            if (rank == 0 && x == 0 && z == 0) {
                                fprintf(stderr, "[WARNING] Failed to get PLT eigenmode for (kx=%d, ky=%d, kz=%d), using normal computation\n",
                                        kx, ky, kz);
                            }
                        }
                    }
                }
                
                // ========== STEP 3.5: Compute PLT growth rate f and rescale (before computing F, G, H) ==========
                // f is the logarithmic derivative of the growth factor that scales velocities
                // When PLT is enabled: f = (sqrt(1. + 24 * e.val * f_cluster) - 1) / 4.
                // When PLT is not enabled: f = 1.0 (default)
                // Skip in density-only mode (qdensity == 2)
                // NOTE: Use outer f variable (declared on 397), not a new inner f!
                // Otherwise the PLT f value goes out of scope before velocity storage.
                if (!just_density && use_plt && params_handle != NULL) {
                    double f_cluster = zeldovich_params_get_f_cluster(params_handle);
                    // f scales the velocities. The corrections are sourced from:
                    // 1) PLT growth rate 2) Addition of a smooth, non-clustering component
                    // to the background (<= NOT A PLT EFFECT)
                    // If PLT true, combine the effects here. Else: apply f_cluster during output.
                    f = (sqrt(1. + 24. * e.val * f_cluster) - 1.) * 0.25;
                    
                    // Compute rescaling if qPLTrescale true
                    // rescale = pow(a_NL/a0, target_f - plt_f)
                    // where plt_f = f (PLT growth rate) and target_f is the continuum growth rate
                    if (qPLTrescale) {
                        double plt_f = f;  // PLT growth rate for this mode
                        rescale = pow(a_NL / a0, target_f - plt_f);
                    }
                }
                
                // Compute factor (used identically in both PLT and non-PLT cases)
                // In zeldovich.cpp: k2 includes fundamental^2, so ik2 = 1/(k2_index * fundamental^2)
                // F = rescale * I * vec * fundamental * ik2 * D
                //   = rescale * I * vec * fundamental / (k2_index * fundamental^2) * D
                //   = rescale * I * vec / (k2_index * fundamental) * D
                // where vec is either e.vec[i] (PLT) or k[i] (non-PLT)
                // Note: our ik2 = 1/k2_index (no fundamental), so we use 1/(k2*fundamental)
                factor = rescale / (k2 * fundamental);
                
                if (use_plt) {
                    // PLT mode: Use eigenvector instead of k-vector
                    #if DEBUG_EIGENVECTOR
                    // Debug: Print eigenvector components for test coordinates
                    int debug_eigen = (rank == 0 && x <= 2 && global_y <= 2 && z <= 2);
                    if (debug_eigen) {
                        fprintf(stderr, "[EIGEN-DEBUG] N=%d Y=%d (x,z)=(%d,%d) k=(%d,%d,%d): e.vec=[%.6f, %.6f, %.6f] e.val=%.6f\n",
                                N, global_y, x, z, kx, ky, kz, e.vec[0], e.vec[1], e.vec[2], e.val);
                        fflush(stderr);
                    }
                    #endif
                    // F = rescale * i * e.vec[0] * fundamental * ik2 * D
                    //   = rescale * i * e.vec[0] * fundamental * ik2 * (D_re + i*D_im)
                    //   = rescale * (-e.vec[0] * fundamental * ik2 * D_im + i * e.vec[0] * fundamental * ik2 * D_re)
                    F[0] = -e.vec[0] * factor * D[1];  // Real part
                    F[1] =  e.vec[0] * factor * D[0];  // Imaginary part
                    
                    G[0] = -e.vec[1] * factor * D[1];
                    G[1] =  e.vec[1] * factor * D[0];
                    
                    H[0] = -e.vec[2] * factor * D[1];
                    H[1] =  e.vec[2] * factor * D[0];
                } else {
                    // Normal operation: Compute F, G, H from D using k-vector
                    // In zeldovich.cpp: k2 includes fundamental^2, so F = I * kx * fundamental * ik2 * D
                    // In our code: k2 doesn't include fundamental^2, so we need to multiply by fundamental
                    // F = rescale * i * kx * fundamental * ik2 * D
                    //   = rescale * i * kx * fundamental * ik2 * (D_re + i*D_im)
                    //   = rescale * (-kx * fundamental * ik2 * D_im + i * kx * fundamental * ik2 * D_re)
                    F[0] = -kx * factor * D[1];  // Real part
                    F[1] =  kx * factor * D[0];  // Imaginary part
                    
                    G[0] = -ky * factor * D[1];
                    G[1] =  ky * factor * D[0];
                    
                    H[0] = -kz * factor * D[1];
                    H[1] =  kz * factor * D[0];
                }
                #endif  // VERIFY_HERMITIAN_SYMMETRY
                }  // End of else block for !just_density
                
                // ========== STEP 4: Store in arrays ==========
                if (just_density) {
                    // Density-only mode: Only store D (density) in Array 0
                    // Array 0: D (density only, no displacement)
                    PRIM_SLICE(0, x, z)[0] = D[0];  // Real = D_re
                    PRIM_SLICE(0, x, z)[1] = D[1];  // Imag = D_im

                    // Test D+iF calculation
                    // PRIM_SLICE(0, x, z)[0] = D[0] - F[1];  // Real = D_re - F_im
                    // PRIM_SLICE(0, x, z)[1] = D[1] + F[0];  // Imag = D_im + F_re
                } else {
                    // Normal mode: Store D+iF, G+iH, and optionally velocities
                    // Array 0: D + i*F (density + X-displacement)
                    // D + i*F = (D[0] + i*D[1]) + i*(F[0] + i*F[1]) = (D[0] - F[1]) + i*(D[1] + F[0])
                    PRIM_SLICE(0, x, z)[0] = D[0] - F[1];  // Real = D_re - F_im
                    PRIM_SLICE(0, x, z)[1] = D[1] + F[0];  // Imag = D_im + F_re
                    
                    // Array 1: G + i*H (Y-displacement + Z-displacement)
                    // G + i*H = (G[0] + i*G[1]) + i*(H[0] + i*H[1]) = (G[0] - H[1]) + i*(G[1] + H[0])
                    PRIM_SLICE(1, x, z)[0] = G[0] - H[1];  // Real = G_re - H_im
                    PRIM_SLICE(1, x, z)[1] = G[1] + H[0];  // Imag = G_im + H_re
                    
                    if (narray >= 4) { // qPLT == true
                        // Array 2: 0 + i*F*f (X-velocity)
                        // 0 + i*(F*f) = 0 + i*((F[0] + i*F[1])*f) = -F[1]*f + i*(F[0]*f)
                        // f is computed above (PLT growth rate if PLT enabled, else 1.0)
                        PRIM_SLICE(2, x, z)[0] = -F[1] * f;  // Real = -F_im * f
                        PRIM_SLICE(2, x, z)[1] = F[0] * f;   // Imag = F_re * f
                        
                        // Array 3: G*f + i*H*f (Y-velocity + Z-velocity)
                        // (G*f) + i*(H*f) = (G[0] - H[1])*f + i*((G[1] + H[0])*f)
                        PRIM_SLICE(3, x, z)[0] = (G[0] - H[1]) * f;  // Real = (G_re - H_im) * f
                        PRIM_SLICE(3, x, z)[1] = (G[1] + H[0]) * f;  // Imag = (G_im + H_re) * f
                    }
                }
                
                // ========== STEP 5: Store conjugates (Zeldovich scheme: conj(D) + i*conj(F)) ==========
                if (just_density) {
                    // Density-only mode: Only store D (density) in Array 0
                    // For conjugate, store conj(D) = (D[0], -D[1])
                    CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0];   // Real = D_re
                    CONJ_SLICE(0, x_mirror, z_mirror)[1] = -D[1];  // Imag = -D_im (conjugate)
                    
                    // Test D+iF calculation
                    // CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0] + F[1];  // Real = D_re + F_im
                    // CONJ_SLICE(0, x_mirror, z_mirror)[1] = F[0] - D[1];  // Imag = F_re - D_im
                } else {
                    // Normal mode: Store conj(D)+i*conj(F), conj(G)+i*conj(H), and optionally velocities
                    // For mode -k, store conj(D) + i*conj(F), NOT the conjugate of (D + i*F)!
                    // conj(D) + i*conj(F) = (D[0] - i*D[1]) + i*(F[0] - i*F[1]) = (D[0] + F[1]) + i*(F[0] - D[1])
                    CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0] + F[1];  // Real = D_re + F_im
                    CONJ_SLICE(0, x_mirror, z_mirror)[1] = F[0] - D[1];  // Imag = F_re - D_im
                    
                    // conj(G) + i*conj(H) = (G[0] - i*G[1]) + i*(H[0] - i*H[1]) = (G[0] + H[1]) + i*(H[0] - G[1])
                    CONJ_SLICE(1, x_mirror, z_mirror)[0] = G[0] + H[1];  // Real = G_re + H_im
                    CONJ_SLICE(1, x_mirror, z_mirror)[1] = H[0] - G[1];  // Imag = H_re - G_im
                    
                    if (narray >= 4) {
                        // Array 2: 0 + i*conj(F*f) = i*conj(F*f)
                        // conj(F*f) = F_re*f - i*F_im*f
                        // i*conj(F*f) = i*(F_re*f - i*F_im*f) = F_im*f + i*F_re*f
                        // f is computed above (PLT growth rate if PLT enabled, else 1.0)
                        CONJ_SLICE(2, x_mirror, z_mirror)[0] = F[1] * f;   // Real = F_im * f
                        CONJ_SLICE(2, x_mirror, z_mirror)[1] = F[0] * f;   // Imag = F_re * f (FIXED: was -F[0])
                        
                        // Array 3: conj(G*f) + i*conj(H*f) = (G[0] + H[1])*f + i*((H[0] - G[1])*f)
                        CONJ_SLICE(3, x_mirror, z_mirror)[0] = (G[0] + H[1]) * f;  // Real = (G_re + H_im) * f
                        CONJ_SLICE(3, x_mirror, z_mirror)[1] = (H[0] - G[1]) * f;  // Imag = (H_re - G_im) * f
                    }
                }
                
            }
            
        }
#if PARALLELIZE_Z_LOOP
            double t_after_loop = omp_get_wtime();
            
            // Print timing diagnostics (only for first Y-slice to reduce output)
            if (global_y == 1 || global_y == 2) {
                #pragma omp critical
                {
                    fprintf(stderr, "[OMP-TIMING] Rank=%d Y=%d Thread=%d/%d z=[%d,%d) vstart=%lld | "
                            "alloc=%.3fms copy=%.3fms advance=%.3fms loop=%.3fms total=%.3fms\n",
                            rank, global_y, tid, nthreads, z_start, z_end, (long long)virtual_start,
                            (t_after_alloc - t_start) * 1000.0,
                            (t_after_copy - t_after_alloc) * 1000.0,
                            (t_after_advance - t_after_copy) * 1000.0,
                            (t_after_loop - t_after_advance) * 1000.0,
                            (t_after_loop - t_start) * 1000.0);
                    fflush(stderr);
                }
            }
            
            if (need_free) free(local_rng_buf);
        }
#endif
    } else {
        // ========== SELF-CONJUGATE: Y=0 or Y=N/2 ==========
        // Reset nskip for self-conjugate case
        // For self-conjugate slices, we process full NxN plane (like zeldovich.cpp)
        if (N < MAX_PPD) {
            nskip = 0;  // Will be accumulated during loops
        } else {
            nskip = 0;
        }
        
        // Zeldovich method: Fill half the plane, mirror the rest
        // Process FULL plane first: z = 0 to N-1, x = 0 to N-1
#if PARALLELIZE_Z_LOOP
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
            int z_start, z_end;
            get_thread_z_range(tid, nthreads, N, &z_start, &z_end);
            
            // Diagnostic timing
            double t_start = omp_get_wtime();
            
            // Use persistent buffer if provided, otherwise allocate
            void* local_rng_buf;
            int need_free = 0;
            if (thread_rng_buffers != NULL && tid < omp_get_max_threads()) {
                local_rng_buf = thread_rng_buffers[tid];
            } else {
                size_t rng_size = zeldovich_ps_rng_buffer_size();
                local_rng_buf = malloc(rng_size);
                need_free = 1;
            }
            
            double t_after_alloc = omp_get_wtime();
            
            zeldovich_ps_get_rng_copy(ps_handle, global_y, local_rng_buf);
            double t_after_copy = omp_get_wtime();
            
            int64_t virtual_start = compute_virtual_position(z_start, 0, N, Nhalf) / 2;
            if (virtual_start > 0)
                zeldovich_ps_advance_rng_buffer(local_rng_buf, virtual_start);
            double t_after_advance = omp_get_wtime();
            
            int64_t nskip = 0;
            for (int z = z_start; z < z_end; z++) {
#else
        // Non-omp-parallel: Process z range sequentially
        for (int z = 0; z < N; z++) {
            void* local_rng_buf = NULL;
#endif
            // RNG skipping: match zeldovich.cpp - skip at Nyquist boundary (z == Nhalf + 1)
            // This applies to ALL Y slices including self-conjugate (Y=0 and Y=N/2)
            if (z == Nhalf + 1 && N < MAX_PPD) {
                // Skip ALL missing z-rows: from z=N to z=MAX_PPD-1
                // Use explicit int64_t casts to avoid integer overflow
                int64_t skip_amount = (int64_t)(MAX_PPD - N) * (int64_t)MAX_PPD;
                nskip += skip_amount; // skip missing z-rows
                
                #if DEBUG_RNG_SKIP
                int log_skip = (global_y <= MAX_DEBUG_COORD) || (global_y == Nhalf - 1);
                if (log_skip) {
                    fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d z=%d: ACCUMULATE skip at Nyquist boundary (self-conj, full plane): +%lld (missing z-rows), total nskip=%lld\n",
                            N, global_y, z, (long long)skip_amount, (long long)nskip);
                    fflush(stderr);
                }
                #endif
            }
            
            int ky = (global_y > Nhalf) ? global_y - N : global_y;
            int kz = (z > Nhalf) ? z - N : z;
            int abs_ky = (ky < 0) ? -ky : ky;
            int abs_kz = (kz < 0) ? -kz : kz;
            
            for (int x = 0; x < N; x++) {
                // RNG skipping: match zeldovich.cpp - skip at Nyquist boundary (x == Nhalf + 1)
                // This applies to ALL Y slices including self-conjugate (Y=0 and Y=N/2)
                if (x == Nhalf + 1 && N < MAX_PPD) {
                    // Skip ALL missing x-values: from x=N to x=MAX_PPD-1
                    int64_t skip_amount = (int64_t)(MAX_PPD - N);
                    nskip += skip_amount; // skip missing x-values in this z-row
                    
                    #if DEBUG_RNG_SKIP
                    int log_skip = (z <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD) ||
                                   (z == Nhalf - 1 && global_y == Nhalf - 1);
                    if (log_skip) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d z=%d x=%d: ACCUMULATE skip at Nyquist boundary (self-conj, full plane): +%lld (missing x-values), total nskip=%lld\n",
                                N, global_y, z, x, (long long)skip_amount, (long long)nskip);
                        fflush(stderr);
                    }
                    #endif
                }
                
                int kx = (x > Nhalf) ? x - N : x;
                int k2_int = kx*kx + ky*ky + kz*kz;  // Integer k^2 for k_cutoff comparison
                double k2 = (double)k2_int;  // Floating-point k^2 for power spectrum
                
                // Nyquist frequency zeroing: Force Nyquist elements to zero for all three axes
                // This matches zeldovich.cpp line 354: abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax
                // The Nyquist frequency (k = N/2) is self-conjugate and doesnt have a separate partner
                // Due to the Y-shift in the reflected shell, we need to zero these to align mirroring expectations
                int abs_kx = (kx < 0) ? -kx : kx;
                int is_nyquist = (abs_kx == Nhalf || abs_ky == Nhalf || abs_kz == Nhalf);
                
                // Generate D using RNG or cgauss()
                // For self-conjugate slices, we still need to check for power spectrum mode
                fftw_complex D;
                if ((k2 == 0.0) || (is_nyquist) || (!CornerModes && (double)k2_int >= k2_cutoff)) {
                    // Zero D for: DC mode, Nyquist frequency, or k_cutoff filtering
                    // This matches zeldovich.cpp line 360-364: zeroing conditions
                    D[0] = D[1] = 0.0;
                    // RNG consistency: When D=0, we skip the RNG call, so accumulate skip
                    nskip++;  // Accumulate skip, will be applied before next cgauss() call
                    #if DEBUG_RNG_SKIP
                    int log_skip = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                   (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_skip) {
                        if (k2 == 0.0) {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (DC mode, self-conj), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, (long long)nskip);
                        } else if (is_nyquist) {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (Nyquist: kx=%d ky=%d kz=%d, self-conj), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, kx, ky, kz, (long long)nskip);
                        } else {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (k_cutoff: k2=%d >= %.1f, self-conj), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, k2_int, k2_cutoff, (long long)nskip);
                        }
                        fflush(stderr);
                    }
                    #endif
                } else if (ps_handle != NULL && params_handle != NULL) {
                    // v15.2: Use zeldovich-PLT power spectrum-weighted Gaussian
                    // Convert k indices to physical wavenumber: k_phys = k_index * fundamental
                    double k2_phys = k2 * fundamental * fundamental;
                    double kmag = sqrt(k2_phys);
                    
                    // zeldovich_ps_cgauss returns double precision, convert to real_t
                    // zeldovich-PLT's v2rng array is sized to ppd/2, so valid indices are 0 to (N/2 - 1)
                    // Since we've already handled global_y == N/2 above, global_y is now < N/2
                    int64_t rng_index = global_y;
                    
                    double D_real, D_imag;
                    #if VERIFY_RNG_CALLS
                    if (nskip > 0 && local_rng_buf == NULL) total_rng_skips += nskip;
                    #endif
                    // Phase 4: use get_cgauss (shared RNG when local_rng_buf==NULL, buffer RNG when parallel)
                    get_cgauss(ps_handle, params_handle, rng_index, kmag, &nskip, local_rng_buf, &D_real, &D_imag);
                    #if VERIFY_RNG_CALLS
                    total_rng_calls++;  // Each cgauss() call uses 2 random numbers
                    #endif
                    D[0] = (real_t)D_real;
                    D[1] = (real_t)D_imag;
                    // Note: For VERIFY_HERMITIAN_SYMMETRY==2, we set D=0 AFTER computing F and H
                } 
                // Compute F, G, H from D
                // Skip F, G, H computation in density-only mode (qdensity == 2)
                fftw_complex F, G, H;
                
                // Declare f_sc before if/else blocks (needed for velocity arrays in all modes)
                double f_sc = 1.0;
                
                // Match zeldovich.cpp: set f=0 when D==0 (line 443)
                if (D[0] == 0.0 && D[1] == 0.0) {
                    f_sc = 0.0;
                }
                
                if (just_density) {
                    // Density-only mode: Set F, G, H to zero (not used)
                    F[0] = F[1] = G[0] = G[1] = H[0] = H[1] = 0.0;

                } else if (k2 == 0.0) {
                    F[0] = F[1] = G[0] = G[1] = H[0] = H[1] = 0.0;
                    f_sc = 0.0;  // Match zeldovich.cpp: set f=0 when k2==0 (D==0)
                } else {
                    double ik2 = 1.0 / k2;
                    
                    #if defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 1
                    // Verification mode 1: Set F=0 and H=0 to test Hermitian symmetry
                    // With F=0 and H=0, only D and G remain, and result should be purely real
                    F[0] = 0.0;
                    F[1] = 0.0;
                    G[0] = -ky * ik2 * D[1];
                    G[1] =  ky * ik2 * D[0];
                    H[0] = 0.0;
                    H[1] = 0.0;
                    #elif defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 2
                    // Verification mode 2: Compute F and H from D, then set D=0 and G=0
                    // This makes Array 0 and Array 1 purely imaginary (real parts = 0)
                    // After 3D FFT, the result should be purely imaginary (real parts = 0)
                    // IMPORTANT: Must use the same factor computation as normal mode to ensure
                    // F is computed identically. Otherwise the test is invalid.
                    // Get fundamental wavenumber (needed for factor computation)
                    double fundamental_sc = 1.0;
                    if (params_handle != NULL) {
                        fundamental_sc = zeldovich_params_get_fundamental(params_handle);
                    }
                    // For VERIFY_HERMITIAN_SYMMETRY=2, we don't use PLT rescaling, so rescale=1.0
                    double rescale_sc = 1.0;
                    // Compute factor the same way as normal mode
                    double factor_sc = rescale_sc / (k2 * fundamental_sc);
                    
                    // Now compute F and H using the same formula as normal mode
                    F[0] = -kx * factor_sc * D[1];
                    F[1] =  kx * factor_sc * D[0];
                    
                    G[0] = 0.0;  // Set G=0
                    G[1] = 0.0;
                    
                    H[0] = -kz * factor_sc * D[1];
                    H[1] =  kz * factor_sc * D[0];
                    
                    // Now set D=0 (after F and H are computed)
                    D[0] = 0.0;
                    D[1] = 0.0;
                    #else
                    // Normal operation: Compute F, G, H from D
                    // double ik2 = 1.0 / k2; --> UNCOMMENT LATER! AFTER D_0 TEST
                    // Check if PLT is enabled
                    int qPLT_sc = 0;
                    // Default fundamental = 1.0 
                    double fundamental_sc = 1.0;
                    eigenmode e_sc;
                    int use_plt_sc = 0;
                    
                    // ========== STEP 3.5: Compute PLT growth rate f and rescale (before computing F, G, H) ==========
                    // f is the logarithmic derivative of the growth factor that scales velocities
                    // When PLT is enabled: f = (sqrt(1. + 24 * e.val * f_cluster) - 1) / 4.
                    // When PLT is not enabled: f = 1.0 (default)
                    // Skip in density-only mode (qdensity == 2)
                    // f_sc already declared above, update it here if PLT is enabled
                    double rescale_sc = 1.0;
                    
                    // Get fundamental wavenumber (needed for both PLT and non-PLT cases)
                    if (!just_density && params_handle != NULL) {
                        fundamental_sc = zeldovich_params_get_fundamental(params_handle);
                        qPLT_sc = zeldovich_params_get_qPLT(params_handle);
                        if (qPLT_sc) {
                            // Get PLT eigenmode for this k-vector
                            // Convert kx, ky, kz to array indices for plt_get_eigenmode
                            int ikx = (kx < 0) ? N + kx : kx;
                            int iky = (ky < 0) ? N + ky : ky;
                            int ikz = (kz < 0) ? N + kz : kz;
                            // Handle z index (only positive half-space stored)
                            if (ikz > N / 2) ikz = N - ikz;
                            
                            if (plt_get_eigenmode(ikx, iky, ikz, (int64_t)N, &e_sc) == 0) {
                                use_plt_sc = 1;
        
                                // Eigenval normalization
                                if (kz < 0) {
                                    e_sc.vec[2] = -e_sc.vec[2];
                                }
                                
                                // Normalize eigenvector (interpolation might not preserve |e_sc| = 1)
                                double e_sc_mag = sqrt(e_sc.vec[0] * e_sc.vec[0] + e_sc.vec[1] * e_sc.vec[1] + e_sc.vec[2] * e_sc.vec[2]);
                                if (e_sc_mag > 0.0) {
                                    e_sc.vec[0] /= e_sc_mag;
                                    e_sc.vec[1] /= e_sc_mag;
                                    e_sc.vec[2] /= e_sc_mag;
                                }
                                
                                // Apply normalization: norm = k2 / (k * e_sc)
                                double k_dot_e_sc = kx * e_sc.vec[0] + ky * e_sc.vec[1] + kz * e_sc.vec[2];
                                double norm_sc = (k2 > 0.0 && k_dot_e_sc != 0.0) ? k2 / k_dot_e_sc : 0.0;
                                if (!isfinite(norm_sc)) norm_sc = 0.0;
                                
                                // Scale eigenvector by norm
                                e_sc.vec[0] *= norm_sc;
                                e_sc.vec[1] *= norm_sc;
                                e_sc.vec[2] *= norm_sc;
                                
                                // Compute f and rescale before computing F, G, H
                                double f_cluster = zeldovich_params_get_f_cluster(params_handle);
                                // PLT growth rate: f = (sqrt(1. + 24 * e.val * f_cluster) - 1) / 4.
                                f_sc = (sqrt(1. + 24. * e_sc.val * f_cluster) - 1.) * 0.25;
                                
                                // Compute rescaling if qPLTrescale is enabled
                                if (qPLTrescale) {
                                    double plt_f = f_sc;  // PLT growth rate for this mode
                                    rescale_sc = pow(a_NL / a0, target_f - plt_f);
                                }
                            } else {
                                // If eigenmode lookup fails, fall back to normal computation
                                if (rank == 0 && x == 0 && z == 0) {
                                    fprintf(stderr, "[WARNING] Failed to get PLT eigenmode for (kx=%d, ky=%d, kz=%d), using normal computation\n",
                                            kx, ky, kz);
                                }
                            }
                        }
                    }
                    
                    // Compute factor (used identically in both PLT and non-PLT cases)
                    // In zeldovich.cpp: k2 includes fundamental^2, so factor = rescale / (k2 * fundamental)
                    // where vec is either e.vec[i] (PLT) or k[i] (non-PLT)
                    double factor = rescale_sc / (k2 * fundamental_sc);
                    
                    if (use_plt_sc) {
                        // PLT mode: Use eigenvector instead of k-vector
                        // F = rescale * i * e.vec[0] * fundamental * ik2 * D
                        F[0] = -e_sc.vec[0] * factor * D[1];
                        F[1] =  e_sc.vec[0] * factor * D[0];
                        
                        G[0] = -e_sc.vec[1] * factor * D[1];
                        G[1] =  e_sc.vec[1] * factor * D[0];
                        
                        H[0] = -e_sc.vec[2] * factor * D[1];
                        H[1] =  e_sc.vec[2] * factor * D[0];
                    } else {
                        // Normal operation: Compute F, G, H from D using k-vector
                        // In zeldovich.cpp: k2 includes fundamental^2, so F = I * kx * fundamental * ik2 * D
                        // In our code: k2 doesn't include fundamental^2, so we need to multiply by fundamental
                        F[0] = -kx * factor * D[1];
                        F[1] =  kx * factor * D[0];
                        G[0] = -ky * factor * D[1];
                        G[1] =  ky * factor * D[0];
                        H[0] = -kz * factor * D[1];
                        H[1] =  kz * factor * D[0];
                    }
                    #endif
                }
                
                // Compute mirror indices for self-conjugate slices (like zeldovich's zHer, xHer)
                int x_mirror = (x == 0) ? 0 : N - x;
                int z_mirror = (z == 0) ? 0 : N - z;
                
                // Store in arrays (Zeldovich packing scheme)
                // For self-conjugate slices: store D+i*F in primary slice,
                // and conj(D)+i*conj(F) in conjugate slice at mirrored positions
                // (like zeldovich stores in slab and slabHer)
                if (just_density) {
                    // Density-only mode: Only store D (density) in Array 0
                    PRIM_SLICE(0, x, z)[0] = D[0];  // Real = D_re
                    PRIM_SLICE(0, x, z)[1] = D[1];  // Imag = D_im

                    CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0];  // Real = D_re
                    CONJ_SLICE(0, x_mirror, z_mirror)[1] = -D[1];  // Imag = -D_im (conjugate)

                    // Test D+iF calculation
                    // PRIM_SLICE(0, x, z)[0] = D[0] - F[1];  // Real = D_re - F_im
                    // PRIM_SLICE(0, x, z)[1] = D[1] + F[0];  // Imag = D_im + F_re
                    
                    // Conjugate slice (at mirrored position): conj(D) + i*conj(F) = (D[0] + F[1]) + i*(F[0] - D[1])
                    // CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0] + F[1];  // Real = D_re + F_im
                    // CONJ_SLICE(0, x_mirror, z_mirror)[1] = F[0] - D[1];  // Imag = F_re - D_im
                } else {
                    // Normal mode: Store D+iF, G+iH, and optionally velocities
                    // Primary slice: Array 0: D + i*F = (D[0] - F[1]) + i*(D[1] + F[0])
                    PRIM_SLICE(0, x, z)[0] = D[0] - F[1];
                    PRIM_SLICE(0, x, z)[1] = D[1] + F[0];
                    
                    // Array 1: G + i*H = (G[0] - H[1]) + i*(G[1] + H[0])
                    PRIM_SLICE(1, x, z)[0] = G[0] - H[1];
                    PRIM_SLICE(1, x, z)[1] = G[1] + H[0];
                    
                    // Conjugate slice (at mirrored position): conj(D)+i*conj(F), conj(G)+i*conj(H)
                    // conj(D) + i*conj(F) = (D[0] + F[1]) + i*(F[0] - D[1])
                    CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0] + F[1];  // Real = D_re + F_im
                    CONJ_SLICE(0, x_mirror, z_mirror)[1] = F[0] - D[1];  // Imag = F_re - D_im
                    
                    // conj(G) + i*conj(H) = (G[0] + H[1]) + i*(H[0] - G[1])
                    CONJ_SLICE(1, x_mirror, z_mirror)[0] = G[0] + H[1];  // Real = G_re + H_im
                    CONJ_SLICE(1, x_mirror, z_mirror)[1] = H[0] - G[1];  // Imag = H_re - G_im
                    
                    if (narray >= 4) {
                        // Primary slice: Array 2: 0 + i*F*f = -F[1]*f + i*(F[0]*f)
                        // f_sc is computed above (PLT growth rate if PLT enabled, else 1.0)
                        PRIM_SLICE(2, x, z)[0] = -F[1] * f_sc;
                        PRIM_SLICE(2, x, z)[1] = F[0] * f_sc;
                        // Array 3: G*f + i*H*f = (G[0] - H[1])*f + i*((G[1] + H[0])*f)
                        PRIM_SLICE(3, x, z)[0] = (G[0] - H[1]) * f_sc;
                        PRIM_SLICE(3, x, z)[1] = (G[1] + H[0]) * f_sc;
                        
                        // Conjugate slice (at mirrored position): i*conj(F*f), conj(G*f)+i*conj(H*f)
                        // i*conj(F*f) = i*(F_re*f - i*F_im*f) = F_im*f + i*F_re*f
                        CONJ_SLICE(2, x_mirror, z_mirror)[0] = F[1] * f_sc;   // Real = F_im * f
                        CONJ_SLICE(2, x_mirror, z_mirror)[1] = F[0] * f_sc;   // Imag = F_re * f
                        
                        // conj(G*f) + i*conj(H*f) = (G[0] + H[1])*f + i*((H[0] - G[1])*f)
                        CONJ_SLICE(3, x_mirror, z_mirror)[0] = (G[0] + H[1]) * f_sc;  // Real = (G_re + H_im) * f
                        CONJ_SLICE(3, x_mirror, z_mirror)[1] = (H[0] - G[1]) * f_sc;  // Imag = (H_re - G_im) * f
                    }
                }
                
            }
        }
#if PARALLELIZE_Z_LOOP
            double t_after_loop = omp_get_wtime();
            
            // Print timing diagnostics for self-conjugate slices (Y=0 and Y=N/2)
            if (global_y == 0) {
                #pragma omp critical
                {
                    fprintf(stderr, "[OMP-TIMING] Rank=%d Y=%d (self-conj) Thread=%d/%d z=[%d,%d) vstart=%lld | "
                            "alloc=%.3fms copy=%.3fms advance=%.3fms loop=%.3fms total=%.3fms\n",
                            rank, global_y, tid, nthreads, z_start, z_end, (long long)virtual_start,
                            (t_after_alloc - t_start) * 1000.0,
                            (t_after_copy - t_after_alloc) * 1000.0,
                            (t_after_advance - t_after_copy) * 1000.0,
                            (t_after_loop - t_after_advance) * 1000.0,
                            (t_after_loop - t_start) * 1000.0);
                    fflush(stderr);
                }
            }
            
            if (need_free) free(local_rng_buf);
        }
#endif
        
        // Post-processing: Mirror first half to second half (match zeldovich.cpp lines 555-573)
        // This is done AFTER processing the full plane, matching zeldovich.cpp behavior
        // zeldovich.cpp: for (z = 0; z < ppdhalf; z++) where ppdhalf = N/2
        // Apply to ALL self-conjugate slices (Y=0 and Y=N/2) since they use the same processing method
        if (global_y == 0 || global_y == Nhalf) {
            for (int z = 0; z < Nhalf; z++) {
                int z_mirror = (z == 0) ? 0 : N - z;
                // Match zeldovich.cpp: xmax = ppdhalf for z=0, ppd for z>0
                // zeldovich.cpp: int xmax = (z == 0 ? ppdhalf : ppd);
                int x_max = (z == 0 ? Nhalf : N);
                for (int x = 0; x < x_max; x++) {
                    int x_mirror = (x == 0) ? 0 : N - x;
                    // Mirror by taking complex conjugate (match zeldovich.cpp line 566-567)
                    // zeldovich.cpp copies from slabHer (which contains conjugates) to slab
                    // For self-conjugate slices, we need to conjugate when mirroring to preserve Hermitian symmetry
                    // f(kx, 0, kz) = conj(f(-kx, 0, -kz)) for self-conjugate slice at Y=0 (normal operation and Mode 1)
                    // For Mode 2 (purely imaginary result), we need anti-Hermitian: f(kx, 0, kz) = -conj(f(-kx, 0, -kz))
                    #if defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 2
                    // Mode 2: Anti-Hermitian symmetry for purely imaginary result
                    // f(-k) = -conj(f(k)) means we negate the conjugate
                    if (x == 0 && z == 0) {
                        fprintf(stderr, "[MIRROR-DEBUG] Y=%d: Mode 2 branch executing (VERIFY_HERMITIAN_SYMMETRY=2)\n", global_y);
                        fflush(stderr);
                    }
                    for (int a = 0; a < narray; a++) {
                        // Negated complex conjugate: -(a + i*b)* = -(a - i*b) = -a + i*b
                        PRIM_SLICE(a, x_mirror, z_mirror)[0] = CONJ_SLICE(a, x_mirror, z_mirror)[0];  // -PRIM_SLICE(a, x, z)[0];  // Real part (negated)
                        PRIM_SLICE(a, x_mirror, z_mirror)[1] = CONJ_SLICE(a, x_mirror, z_mirror)[1]; // PRIM_SLICE(a, x, z)[1];   // Imaginary part (same)
                    }
                    #else
                    // Normal operation and Mode 1: Hermitian symmetry for purely real result
                    // Match zeldovich's scheme: copy conj(D)+i*conj(F) from conjugate_slices to primary_slices
                    // This copies the conjugate values that were stored in conjugate_slices during main loop
                    // (like zeldovich copies from slabHer to slab at lines 560-561)
                    if (x == 0 && z == 0) {
                        #if defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 1
                        fprintf(stderr, "[MIRROR-DEBUG] Y=%d: Mode 1 branch executing (VERIFY_HERMITIAN_SYMMETRY=1)\n", global_y);
                        #else
                        fprintf(stderr, "[MIRROR-DEBUG] Y=%d: Normal operation branch executing (VERIFY_HERMITIAN_SYMMETRY not defined or != 1,2)\n", global_y);
                        #endif
                        fflush(stderr);
                    }
                    for (int a = 0; a < narray; a++) {
                        // Copy from conjugate slice (which contains conj(D)+i*conj(F) at mirror positions)
                        // to primary slice at mirror positions
                        PRIM_SLICE(a, x_mirror, z_mirror)[0] = CONJ_SLICE(a, x_mirror, z_mirror)[0];  // Real part
                        PRIM_SLICE(a, x_mirror, z_mirror)[1] = CONJ_SLICE(a, x_mirror, z_mirror)[1];  // Imaginary part
                    }
                    #endif
                }
            }
            // Set origin to zero (match zeldovich.cpp line 572)
            // This applies to both Y=0 and Y=N/2 self-conjugate slices
            for (int a = 0; a < narray; a++) {
                PRIM_SLICE(a, 0, 0)[0] = 0.0;
                PRIM_SLICE(a, 0, 0)[1] = 0.0;
            }
        }
        // Set imaginary parts to 0 at special points
        for (int a = 0; a < narray; a++) {
            PRIM_SLICE(a, 0, 0)[1] = 0.0;
            PRIM_SLICE(a, Nhalf, 0)[1] = 0.0;
            PRIM_SLICE(a, 0, Nhalf)[1] = 0.0;
            PRIM_SLICE(a, Nhalf, Nhalf)[1] = 0.0;
        }
    }
    
    t_zloop_end = omp_get_wtime();
    
    // Verify Hermitian symmetry BEFORE 2D FFT
    // STAGE 7: Updated to check all arrays independently
    #if VERIFY_HERMITIAN_SYMMETRY
    for (int a = 0; a < narray; a++) {
        fftw_complex_t *prim_array = &PRIM_SLICE(a, 0, 0);
        fftw_complex_t *conj_array = (y_mirror != global_y) ? &CONJ_SLICE(a, 0, 0) : NULL;
        verify_initial_fourier_hermitian_symmetry(N, prim_array, global_y, y_mirror, conj_array, a);
    }
    #endif
    
    t_verify_end = omp_get_wtime();
    
    // Apply 2D FFT to all arrays independently (hybrid: outer parallelism over arrays, inner FFTW threading)
    // Limit outer threads to narray to avoid oversubscription (each FFT uses multiple inner threads)
    // Note: With OMP_MAX_ACTIVE_LEVELS=1 (default), inner FFTW threads are serialized,
    // resulting in narray concurrent single-threaded FFTs - which is actually optimal!
    #pragma omp parallel for num_threads(narray)
    for (int a = 0; a < narray; a++) {
        // Primary slice
        fftw_complex_t *prim_array_start = &PRIM_SLICE(a, 0, 0);
        
        FFTW_EXECUTE_DFT(plan_2d, prim_array_start, prim_array_start);
        
        // Conjugate slice (if not self-conjugate)
        if (y_mirror != global_y) {
            fftw_complex_t *conj_array_start = &CONJ_SLICE(a, 0, 0);

            FFTW_EXECUTE_DFT(plan_2d, conj_array_start, conj_array_start);
            
            // // DEBUG: Memory guard after conjugate FFT
            // #if DEBUG_PRINTS
            // if (rank < 4 && global_y <= 2) {
            //     // Verify buffer is still valid after FFT
            //     size_t conj_array_size = (size_t)N * N * sizeof(fftw_complex_t);
            //     volatile char checksum = 0;
            //     for (size_t i = 0; i < conj_array_size && i < 1024*1024; i += 4096) {
            //         checksum ^= ((char*)conj_array_start)[i];
            //     }
            //     fprintf(stderr, "[Rank %d] After conjugate FFT (array %d, Y=%d): buffer OK (checksum=%d)\n",
            //            rank, a, global_y, (int)checksum);
            //     fflush(stderr);
            // }
            // #endif
        }
    }
    
    t_fft_end = omp_get_wtime();
    
    // ========== DIAGNOSTIC: Print per-Y timing breakdown ==========
    // Print for first few Y-slices and occasionally thereafter to avoid flooding
    if (global_y <= 3 || (global_y % 64 == 0)) {
        fprintf(stderr, "[SLICE-TIMING] Rank=%d Y=%d/%d | setup=%.3fms zloop=%.3fms verify=%.3fms fft=%.3fms total=%.3fms\n",
                rank, global_y, y_mirror,
                (t_setup_end - t_func_start) * 1000.0,
                (t_zloop_end - t_setup_end) * 1000.0,
                (t_verify_end - t_zloop_end) * 1000.0,
                (t_fft_end - t_verify_end) * 1000.0,
                (t_fft_end - t_func_start) * 1000.0);
        fflush(stderr);
    }
    
    // ========== Apply any remaining nskip at end of function ==========
    // See zeldovich.cpp: Pk.v2rng[y].advance(2 * nskip)
    // Ensures that exactly MAX_PPD * MAX_PPD complex numbers (2 * MAX_PPD * MAX_PPD real numbers)
    // were consumed/skipped for this y-row
    if (nskip > 0) {
        if (ps_handle != NULL && params_handle != NULL) {
            int64_t rng_index = global_y;
            zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, nskip);
            #if VERIFY_RNG_CALLS
            total_rng_skips += nskip;
            #endif
        }
        nskip = 0;  // Reset after applying
    }
    
    // // ========== RNG VERIFICATION: Verify total RNG calls match expected ==========

    // #if VERIFY_RNG_CALLS
    // {
    //     int64_t max_ppd_val = (int64_t)MAX_PPD;
    //     int64_t expected_total = 2LL * max_ppd_val * max_ppd_val;
    //     int64_t actual_total = 2LL * total_rng_calls + 2LL * total_rng_skips;
        
    //     if (actual_total != expected_total) {
    //         fprintf(stderr, 
    //                 "[RNG-VERIFY ERROR] Rank %d, Y=%d: Expected %ld random numbers (MAX_PPD=%ld), got %ld "
    //                 "(calls=%ld * 2 = %ld, skips=%ld * 2 = %ld)\n",
    //                 rank, global_y, (long)expected_total, (long)max_ppd_val, (long)actual_total,
    //                 (long)total_rng_calls, (long)(2LL * total_rng_calls),
    //                 (long)total_rng_skips, (long)(2LL * total_rng_skips));
    //         fflush(stderr);
    //         // Don't abort in production, but warn
    //         #ifdef DEBUG
    //         assert(actual_total == expected_total);
    //         #endif
    //     } else if (rank == 0 && global_y <= 2) {
    //         // Print verification for first few slices in debug mode
    //         fprintf(stderr,
    //                 "[RNG-VERIFY OK] Rank %d, Y=%d: Total=%ld (calls=%ld, skips=%ld)\n",
    //                 rank, global_y, (long)actual_total,
    //                 (long)total_rng_calls, (long)total_rng_skips);
    //         fflush(stderr);
    //     }
    // }
    // #endif
    
    // Clean up local macros
    #undef PRIM_SLICE
    #undef CONJ_SLICE
}

