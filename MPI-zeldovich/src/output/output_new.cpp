#include <dirent.h>
#include <libgen.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>
#include <mutex>
#include <string.h>
#include <filesystem>
#include <vector>
#include <string>

#include <fmt/base.h>
#include <fmt/format.h>

// Include zeldovich-PLT headers to get Complx, Parameters, OutputType, and particle types
#include <zeldovich.h>
#include <parameters.h>
#include <output.h>  // For OutputType enum and particle struct definitions
#include <complex> 
#include <STimer.h>  // zeldovich-PLT's STimer (angle brackets so reassembly build gets zeldovich's, not src/STimer.h)
#include "output_new.h" 
#include "utils/decomposition.h"

namespace fs = std::filesystem;

// Unified indexing macro - works for both full range and local range
// Takes stride as parameter: use param.ppd for full range, k_extent for local range
// Formula: _slab[(_k) + (_stride) * (_j)] where:
//   - _j = Y coordinate (j index)
//   - _k = X coordinate (k_local index)
//   - _stride = stride in X direction (ppd for full range, k_extent for local range)
#define JK(_slab, _j, _k, _stride) _slab[(_k) + (_stride) * (_j)]

// Legacy macros for backward compatibility
#define YX(_slab, _y, _x, _ppd) JK(_slab, _y, _x, _ppd)
#define YX_LOCAL(_slab, _j, _k_local, _k_extent) JK(_slab, _j, _k_local, _k_extent)

// Unified slab indexing macro - handles both [x][y] and [y][x] layouts
// Parameters:
//   _slab: pointer to slab data
//   _j: Y coordinate (j index)
//   _k_local: X coordinate local index (k_local)
//   _use_xy: true for [x][y] layout (ZSLAB), false for [y][x] layout (JK)
//   _stride: stride for [y][x] layout (param.ppd for full range, k_extent for local range)
//   _N: grid size (param.ppd) for [x][y] layout
#define INDEX_SLAB(_slab, _j, _k_local, _use_xy, _stride, _N) \
    ((_use_xy) ? &(_slab)[(_k_local) * (_N) + (_j)] : &JK(_slab, _j, _k_local, _stride))

// Global maxima of the particle displacements
double max_disp[3];

double density_variance;
void *output_tmp;       // output buffer
float *densoutput_tmp;  // dens output buffer
FILE *densfp;           // dens output fp
STimer outtimer;
size_t output_bytes_written = 0;
std::mutex output_mutex;

OutputType param_icformat;
size_t sizeof_outputtype;

// ====================================================================================
// INTERNAL: Unified particle writing function
// ====================================================================================
// Handles both full-range (WriteParticlesSlab_new) and local-range (WriteParticlesSlab_range)
// cases by using appropriate indexing macro and output strategy.
//
// Parameters:
//   - is_full_range: true for full range (ppd x ppd), false for local range (ppd x k_extent)
//   - use_global_buffers: true to use global output_tmp/densoutput_tmp, false to allocate local
//   - rank: MPI rank (used for per-rank file naming, -1 for full-range mode)
//   - k_start_global: Global X start (0 for full range, rank's X start for local range)
//   - k_extent: Number of X values (ppd for full range, rank's X extent for local range)
//   - use_x_y_layout: true for [x][y] layout (ZSLAB format), false for [y][x] layout (JK format)
//                      For [x][y] layout, uses param.ppd for indexing: slab[k_local * param.ppd + j]
// ====================================================================================

static void WriteParticlesSlab_unified(
   int i,                  // i = Z index (legacy z)
   Complx *slab1,
   Complx *slab2,
   Complx *slab3,
   Complx *slab4,
   Parameters &param,
   bool is_full_range,     // true: full range (ppd x ppd), false: local range
   bool use_global_buffers, // true: use global buffers, false: allocate local
   int rank,               // MPI rank (-1 for full-range mode)
   int k_start_global,     // Global X start
   int k_extent,           // Number of X values
   bool use_x_y_layout,    // true: [x][y] layout (ZSLAB), false: [y][x] layout (JK)
   int N                   // Grid size (for ZSLAB format indexing - must match data layout)
) {
    STimer thisouttimer;
    thisouttimer.Start();
    double thisdensity_variance = 0.0;
    int just_density = param.qdensity == 2;

    int j;
    int k_local;
    double pos[3], vel[3], dens;
    double norm, densitynorm, vnorm;

    norm = 1.0;
    densitynorm = 1.0;
    dens = 0.;

    if (param.qPLT) {
        vnorm = 1.0;
    } else {
        vnorm = (sqrt(1. + 24 * param.f_cluster) - 1) * .25;
    }

    // Determine loop bounds and indexing
    // For ZSLAB format ([x][y] layout), use N (grid size) for Y bounds
    // For JK format ([y][x] layout), use param.ppd for Y bounds
    int j_max = use_x_y_layout ? N : param.ppd; // change later: N = param.ppd
    int k_max = is_full_range ? param.ppd : k_extent;
    int64_t num_particles = (int64_t)j_max * (int64_t)k_max;

    // Buffer management
    void *particle_buffer = NULL;
    float *density_buffer = NULL;
    bool owns_buffers = false;

    if (use_global_buffers) {
        // Full-range mode: use global buffers
        particle_buffer = output_tmp;
        density_buffer = densoutput_tmp;
        owns_buffers = false;
    } else {
        // Local-range mode: allocate local buffers
        if (!just_density) {
            switch (param_icformat) {
                case OUTPUT_RVDOUBLEZEL:
                    particle_buffer = new RVdoubleZelParticle[num_particles];
                    break;
                case OUTPUT_RVZEL:
                    particle_buffer = new RVZelParticle[num_particles];
                    break;
                case OUTPUT_ZEL:
                    particle_buffer = new ZelParticle[num_particles];
                    break;
                case OUTPUT_ZEL_SIMPLE:
                    particle_buffer = new ZelSimpleParticle[num_particles];
                    break;
                default:
                    fmt::print(stderr, "Error: unknown ICFormat \"{:s}\". Aborting.\n", param.ICFormat);
                    exit(1);
            }
        }
        if (param.qdensity) {
            density_buffer = new float[num_particles];
        }
        owns_buffers = true;
    }

    int64_t count = 0;

    // Main loop: iterate over all j and k values
    for (j = 0; j < j_max; j++) {
        for (k_local = 0; k_local < k_max; k_local++) {
            int k_value;  // The k value to use (local for full range, global for local range)
            if (is_full_range) {
                k_value = k_local;  // Full range: k is just the loop index
            } else {
                k_value = k_start_global + k_local;  // Local range: calculate global k
            }

            // Access data using layout-appropriate indexing macro
            // [x][y] layout: slab[k_local * N + j] (ZSLAB format - stride N in x direction)
            // [y][x] layout: slab[k_local + stride * j] (JK format)
            int stride = is_full_range ? param.ppd : k_extent;
            // For ZSLAB format, use N (grid size) for stride; for JK format, use param.ppd
            int index_stride = use_x_y_layout ? N : param.ppd;
            Complx *slab1_val = INDEX_SLAB(slab1, j, k_local, use_x_y_layout, stride, index_stride);
            Complx *slab2_val = INDEX_SLAB(slab2, j, k_local, use_x_y_layout, stride, index_stride);
            Complx *slab3_val = slab3 ? INDEX_SLAB(slab3, j, k_local, use_x_y_layout, stride, index_stride) : NULL;
            Complx *slab4_val = slab4 ? INDEX_SLAB(slab4, j, k_local, use_x_y_layout, stride, index_stride) : NULL;

            // DEBUG: Print indexing for first few values (N=4, i=0 only)
            if (i == 0 && param.ppd <= 4 && j < 2 && k_local < 2) {
                int idx_used = use_x_y_layout ? (k_local * index_stride + j) : (k_local + stride * j);
                fprintf(stderr, "[DEBUG-WRITE] i=%d j=%d(Y) k_local=%d(X) use_xy=%d stride=%d index_stride=%d idx=%d val=%.10e\n",
                       i, j, k_local, use_x_y_layout, stride, index_stride, idx_used, (double)real(*slab1_val));
            }

            dens = real(*slab1_val) * densitynorm;
            if (!just_density) {
                pos[0] = imag(*slab2_val) * norm;
                pos[1] = real(*slab2_val) * norm;
                pos[2] = imag(*slab1_val) * norm;
                if (param.qPLT) {
                    vel[0] = slab4_val ? imag(*slab4_val) * vnorm : 0.0;
                    vel[1] = slab4_val ? real(*slab4_val) * vnorm : 0.0;
                    vel[2] = slab3_val ? imag(*slab3_val) * vnorm : 0.0;
                } else {
                    vel[0] = imag(*slab2_val) * vnorm;
                    vel[1] = real(*slab2_val) * vnorm;
                    vel[2] = imag(*slab1_val) * vnorm;
                }

                if (param.qascii) {
                    // ASCII output: full-range uses provided FILE*, local-range uses stdout
                    if (is_full_range) {
                        // Note: output FILE* is handled by caller in WriteParticlesSlab_new
                        // For now, use stdout for both (can be enhanced if needed)
                        fmt::print(stdout,
                           "{:d} {:d} {:d} {:f} {:f} {:f} {:f} {:f} {:f} {:f}\n",
                           i, j, k_value, pos[0], pos[1], pos[2], dens, vel[0], vel[1], vel[2]);
                    } else {
                        fmt::print(stdout,
                           "{:d} {:d} {:d} {:f} {:f} {:f} {:f} {:f} {:f} {:f}\n",
                           i, j, k_value, pos[0], pos[1], pos[2], dens, vel[0], vel[1], vel[2]);
                    }
                } else {
                    switch (param_icformat) {
                        case OUTPUT_RVDOUBLEZEL: {
                            RVdoubleZelParticle out;
                            out.i = i;
                            out.j = j;
                            out.k = k_value;
                            out.displ[0] = pos[0];
                            out.displ[1] = pos[1];
                            out.displ[2] = pos[2];
                            out.vel[0] = vel[0];
                            out.vel[1] = vel[1];
                            out.vel[2] = vel[2];
                            ((RVdoubleZelParticle *) particle_buffer)[count] = out;
                            break;
                        }
                        case OUTPUT_RVZEL: {
                            RVZelParticle out;
                            out.i = i;
                            out.j = j;
                            out.k = k_value;
                            out.displ[0] = pos[0];
                            out.displ[1] = pos[1];
                            out.displ[2] = pos[2];
                            out.vel[0] = vel[0];
                            out.vel[1] = vel[1];
                            out.vel[2] = vel[2];
                            ((RVZelParticle *) particle_buffer)[count] = out;
                            break;
                        }
                        case OUTPUT_ZEL: {
                            ZelParticle out;
                            out.i = i;
                            out.j = j;
                            out.k = k_value;
                            out.displ[0] = pos[0];
                            out.displ[1] = pos[1];
                            out.displ[2] = pos[2];
                            ((ZelParticle *) particle_buffer)[count] = out;
                            break;
                        }
                        case OUTPUT_ZEL_SIMPLE: {
                            ZelSimpleParticle out;
                            out.displ[0] = pos[0];
                            out.displ[1] = pos[1];
                            out.displ[2] = pos[2];
                            ((ZelSimpleParticle *) particle_buffer)[count] = out;
                            break;
                        }
                        default:
                            fmt::print(stderr, "Error: unknown ICFormat \"{:s}\". Aborting.\n", param.ICFormat);
                            exit(1);
                    }
                }

                // Track global max displacement
                for (int idx = 0; idx < 3; idx++) {
                    max_disp[idx] = fabs(pos[idx]) > fabs(max_disp[idx]) ? pos[idx] : max_disp[idx];
                }
            }

            // Store density
            if (param.qdensity && density_buffer != NULL) {
                density_buffer[count] = (float)dens;
            }
            thisdensity_variance += dens * dens;

            count++;
        }
    }

    assert(count == num_particles);

    // File writing
    int64_t totsize = 0;
    int ic_index = (i * param.cpd) / param.ppd; // rounds up (z-slab index * cpd ) / ppd -> 0, 1, ...

    if (!just_density && particle_buffer != NULL) {
        if (is_full_range) {
            // Full range: write to ic_{ic_index} file
            fs::path fn = param.output_dir / fmt::format("ic_{:d}", ic_index);
            FILE *fp = fopen(fn.c_str(), "ab");
            if (fp) {
                fwrite(particle_buffer, sizeof_outputtype, num_particles, fp);
                fclose(fp);
                totsize += num_particles * (int64_t)sizeof_outputtype;
            } else {
                fmt::print(stderr, "ERROR: Failed to open file {} for writing (ic_index={}, i={}): {}\n",
                          fn.c_str(), ic_index, i, strerror(errno));
            }
        } else {
            // Local range: write to per-rank file
            int k_end_global = k_start_global + k_extent;
            fs::path fn_ic = param.output_dir /
                             fmt::format("ic_rank{:d}_i{:d}_k{:d}_{:d}",
                                         rank, i, k_start_global, k_end_global);
            FILE *fp_ic = fopen(fn_ic.c_str(), "wb");
            if (fp_ic) {
                size_t written = fwrite(particle_buffer, sizeof_outputtype, (size_t)num_particles, fp_ic);
                if (written != (size_t)num_particles) {
                    fmt::print(stderr, "Error: short write to {} (rank={}, wrote {}, expected {}).\n",
                               fn_ic.string(), rank, written, (size_t)num_particles);
                }
                fclose(fp_ic);
                totsize += num_particles * (int64_t)sizeof_outputtype;
            } else {
                fmt::print(stderr, "Error: cannot open {} for writing (rank={}, errno={}).\n",
                           fn_ic.string(), rank, errno);
            }
        }
    }

    // Density file writing
    if (param.qdensity && density_buffer != NULL) {
        if (is_full_range) {
            // Full range: append to global density file
            fwrite(density_buffer, sizeof(*density_buffer) * num_particles, 1, densfp);
            totsize += sizeof(*density_buffer) * num_particles;
        } else {
            // Local range: write to per-rank density file
            int k_end_global = k_start_global + k_extent;
            fs::path fn_dens = param.output_dir /
                               fmt::format("dens_rank{:d}_i{:d}_k{:d}_{:d}",
                                           rank, i, k_start_global, k_end_global);
            FILE *fp_d = fopen(fn_dens.c_str(), "wb");
            if (fp_d) {
                size_t written = fwrite(density_buffer, sizeof(float), (size_t)num_particles, fp_d);
                if (written != (size_t)num_particles) {
                    fmt::print(stderr, "Warning: short density write to {} (rank={}, wrote {}, expected {}).\n",
                               fn_dens.string(), rank, written, (size_t)num_particles);
                }
                fclose(fp_d);
                totsize += num_particles * (int64_t)sizeof(float);
            } else {
                fmt::print(stderr, "Warning: cannot open {} for density write (rank={}, errno={}).\n",
                           fn_dens.string(), rank, errno);
            }
        }
    }

    // Cleanup local buffers if we allocated them
    if (owns_buffers) {
        if (particle_buffer != NULL) {
            switch (param_icformat) {
                case OUTPUT_RVDOUBLEZEL:
                    delete[] (RVdoubleZelParticle *) particle_buffer;
                    break;
                case OUTPUT_RVZEL:
                    delete[] (RVZelParticle *) particle_buffer;
                    break;
                case OUTPUT_ZEL:
                    delete[] (ZelParticle *) particle_buffer;
                    break;
                case OUTPUT_ZEL_SIMPLE:
                    delete[] (ZelSimpleParticle *) particle_buffer;
                    break;
            }
        }
        if (density_buffer != NULL) {
            delete[] density_buffer;
        }
    }

    thisouttimer.Stop();
    output_mutex.lock();
    output_bytes_written += totsize;
    outtimer.increment(thisouttimer.timer);
    density_variance += thisdensity_variance;
    output_mutex.unlock();
}

// ====================================================================================
// WriteParticlesSlab_new
// ====================================================================================
// Writes full-range particle ICs (all X, all Y) for one i-slab.
// Uses global output buffers and writes to ic_{ic_index} files.
// ====================================================================================

void WriteParticlesSlab_new(
   FILE *output,
   int i,   // i = Z index (legacy z)
   Complx *slab1,
   Complx *slab2,
   Complx *slab3,
   Complx *slab4,
   Parameters &param
) {
    // Call unified function with full-range parameters
    // NOTE: Data passed to WriteParticlesSlab_new is in [Y][X] format from reassembled .bin files
    // The reassembled data has layout: full_slab[array * N * N + j * N + k] where j=Y, k=X
    // So element (Y, X) is at: slab1[Y * N + X]
    // 
    // With use_x_y_layout=false (JK format): accesses slab1[X + stride * Y] = slab1[X + N * Y] = slab1[Y * N + X] (correct?)
    // With use_x_y_layout=true (ZSLAB format): accesses slab1[X * N + Y] ≠ slab1[Y * N + X] X
    // 
    // However, if this doesn't work, the issue might be in the reassembly or .bin file format
    WriteParticlesSlab_unified(
        i, slab1, slab2, slab3, slab4, param,
        true,   // is_full_range = true (ppd x ppd)
        true,   // use_global_buffers = true
        -1,     // rank = -1 (not used in full-range mode)
        0,      // k_start_global = 0 (full range starts at 0)
        param.ppd,  // k_extent = ppd (full range)
        true,   // use_x_y_layout = true ([x][y] layout, ZSLAB format) - REVERTED: test if original was correct
        param.ppd  // N = ppd (for full-range mode, ppd should equal grid size)
    );
    (void)output;  // Unused in unified function (handled internally)
}

// ====================================================================================
// WriteParticlesSlab_range
// ====================================================================================
// Same physics as WriteParticlesSlab_new, but restricted to a local X-range:
//   - slab1..4 store local data with layout [j][k_local], where:
//       * j runs over all Y (0 .. ppd-1)
//       * k_local runs over this rank's X-extent (0 .. k_extent-1)
//   - k_start_global specifies the global X-index of k_local = 0
//   - out.k records the global k index (k_global = k_start_global + k_local)
//   - Each rank writes its own files:
//       * Particle ICs: ic_rank{rank}_i{i}_k{k_start_global}_{k_end_global}
//       * Density (if qdensity=1): dens_rank{rank}_i{i}_k{k_start_global}_{k_end_global}
//     where:
//       - i = Z-slab index (z value)
//       - k_start_global = global X start index
//       - k_end_global = k_start_global + k_extent
//
// This function does NOT use the global output_tmp buffer; it allocates a
// per-call buffer sized for ppd * k_extent particles.
// Density planes are written to per-rank files when param.qdensity is enabled
// (consistent with per-rank file strategy used in the MPI module).
// ====================================================================================

// ====================================================================================
// WriteParticlesSlab_range - Overloaded versions
// ====================================================================================
// Version 1: Accepts data in [array][x][y] layout (ZSLAB format) - NO TRANSPOSE NEEDED
// Version 2: Accepts data in [y][x] layout (JK format) - for backward compatibility
// ====================================================================================
// version 1 takes as input slab_data and 4 slab pointers are made from it, while version 2 takes as input slab1, slab2, slab3, slab4. [y][x] vs [x][y]
// Version 1: [array][x][y] layout (ZSLAB format) - eliminates transpose
void WriteParticlesSlab_range(
   int rank,
   int i,                  // i = Z index (legacy z)
   int k_start_global,     // global X start for this rank's extent
   int k_extent,           // number of X values for this rank
   Complx *slab_data,      // Data in [array][x_local][y] layout (ZSLAB format)
   int N,                  // Grid size (ppd)
   int narray,             // Number of arrays (typically 4)
   Parameters &param
) {
    // Extract pointers to each array's data (each array is [x][y] layout)
    Complx *slab1 = &slab_data[0 * k_extent * N];
    Complx *slab2 = &slab_data[1 * k_extent * N];
    Complx *slab3 = (narray > 2) ? &slab_data[2 * k_extent * N] : NULL;
    Complx *slab4 = (narray > 3) ? &slab_data[3 * k_extent * N] : NULL;
    
    // Call unified function with [x][y] layout
    WriteParticlesSlab_unified(
        i, slab1, slab2, slab3, slab4, param,
        false,  // is_full_range = false (local range: ppd x k_extent)
        false,  // use_global_buffers = false (allocate local buffers)
        rank,   // rank = MPI rank (for per-rank file naming)
        k_start_global,  // k_start_global = rank's X start
        k_extent,  // k_extent = rank's X extent
        true,   // use_x_y_layout = true ([x][y] layout, ZSLAB format)
        N       // N = grid size (must match data layout stride)
    );
}

// Version 2: [y][x] layout (JK format) - for backward compatibility
void WriteParticlesSlab_range(
   int rank,
   int i,                  // i = Z index (legacy z)
   int k_start_global,     // global X start for this rank's extent
   int k_extent,           // number of X values for this rank
   Complx *slab1,
   Complx *slab2,
   Complx *slab3,
   Complx *slab4,
   Parameters &param
) {
    // Call unified function with [y][x] layout
    // For JK format, use param.ppd as grid size (Y dimension)
    WriteParticlesSlab_unified(
        i, slab1, slab2, slab3, slab4, param,
        false,  // is_full_range = false (local range: ppd x k_extent)
        false,  // use_global_buffers = false (allocate local buffers)
        rank,   // rank = MPI rank (for per-rank file naming)
        k_start_global,  // k_start_global = rank's X start
        k_extent,  // k_extent = rank's X extent
        false,  // use_x_y_layout = false ([y][x] layout, JK format)
        param.ppd  // N = ppd (for JK format, ppd is the grid size)
    );
}

// ====================================================================================
// WriteParticlesSlab_range_from_zslab
// ====================================================================================
// Variant of WriteParticlesSlab_range that accepts data in [array][x_local][y] layout
// (matching main.cpp's ZSLAB macro format) instead of [y][x_local].
//
// This allows direct particle writing from main.cpp without transposing the z-slab.
//
// Parameters:
//   - slab_data: Pointer to contiguous memory with layout [narray][x_count][N]
//                where data[array_idx][x_idx][y] is at:
//                slab_data[array_idx * x_count * N + x_idx * N + y]
//   - Other parameters same as WriteParticlesSlab_range
//
// Example call from main.cpp:
//   WriteParticlesSlab_range_from_zslab(
//       rank, i, k_start_global, k_extent, local_z_slab, 
//       N, narray, param
//   );
// ====================================================================================

void WriteParticlesSlab_range_from_zslab(
   int rank,
   int i,                  // i = Z index (Zeldovich i, legacy z)
   int k_start_global,     // global X start for this rank's extent (Zeldovich k)
   int k_extent,           // number of X values for this rank (Zeldovich k range)
   Complx *slab_data,      // Data in [array][k_local][j] layout (ZSLAB format: [array][x_local][y])
   int N,                  // Grid size (ppd)
   int narray,             // Number of arrays (typically 4)
   Parameters &param
) {
    STimer thisouttimer;
    thisouttimer.Start();
    double thisdensity_variance = 0.0;
    int just_density            = param.qdensity == 2;  // no displacements

    // Write out one slab of particles for this rank's k-range (X-range)
    // Using Zeldovich notation: i=Z, j=Y, k=X
    int j;       // j = Y coordinate (Zeldovich j, legacy y)
    int k_local; // k = X coordinate local index (Zeldovich k, legacy x_local)
    double pos[3], vel[3], dens;
    double norm, densitynorm, vnorm;

    norm        = 1.0;
    densitynorm = 1.0;
    dens        = 0.;

    if (param.qPLT) {
        vnorm = 1.0;
    } else {
        vnorm = (sqrt(1. + 24 * param.f_cluster) - 1) * .25;
    }

    // Number of particles written by this rank for this slab
    int64_t num_particles = (int64_t)N * (int64_t)k_extent;

    // Allocate local output buffers (per-call, k-range only)
    void *local_output_tmp = NULL;
    float *local_dens_tmp  = NULL;

    if (!just_density) {
        switch (param_icformat) {
            case OUTPUT_RVDOUBLEZEL:
                local_output_tmp = new RVdoubleZelParticle[num_particles];
                break;
            case OUTPUT_RVZEL:
                local_output_tmp = new RVZelParticle[num_particles];
                break;
            case OUTPUT_ZEL:
                local_output_tmp = new ZelParticle[num_particles];
                break;
            case OUTPUT_ZEL_SIMPLE:
                local_output_tmp = new ZelSimpleParticle[num_particles];
                break;
            default:
                fmt::print(
                   stderr, "Error: unknown ICFormat \"{:s}\" in WriteParticlesSlab_range_from_zslab. Aborting.\n",
                   param.ICFormat
                );
                exit(1);
        }
    }

    if (param.qdensity) {
        local_dens_tmp = new float[num_particles];
    }

    int64_t count = 0;

    // Indexing macro for ZSLAB format: [array][k_local][j] (equiv: [array][x_local][y])
    // slab_data[array_idx][k_local][j] = slab_data[array_idx * k_extent * N + k_local * N + j]
    // This matches main.cpp's ZSLAB macro: [array][x_local][y]
#define ZSLAB_LOCAL(_array_idx, _k_local, _j) \
    slab_data[(_array_idx) * k_extent * N + (_k_local) * N + (_j)]

    // Extract array pointers for easier access
    Complx *slab1 = &slab_data[0 * k_extent * N];
    Complx *slab2 = &slab_data[1 * k_extent * N];
    Complx *slab3 = (narray > 2) ? &slab_data[2 * k_extent * N] : NULL;
    Complx *slab4 = (narray > 3) ? &slab_data[3 * k_extent * N] : NULL;

    for (j = 0; j < N; j++) {               // j = Y coordinate (all Y values)
        for (k_local = 0; k_local < k_extent; k_local++) {  // k = X coordinate (local X range)
            int k_global = k_start_global + k_local;

            // Access using ZSLAB format: [k_local][j] (equiv: [x_local][y])
            // slab1[k_local][j] = slab1[k_local * N + j]
            dens = real(slab1[k_local * N + j]) * densitynorm;
            
            if (!just_density) {
                pos[0] = imag(slab2[k_local * N + j]) * norm;
                pos[1] = real(slab2[k_local * N + j]) * norm;
                pos[2] = imag(slab1[k_local * N + j]) * norm;
                
                if (param.qPLT && slab3 && slab4) {
                    vel[0] = imag(slab4[k_local * N + j]) * vnorm;
                    vel[1] = real(slab4[k_local * N + j]) * vnorm;
                    vel[2] = imag(slab3[k_local * N + j]) * vnorm;
                } else {
                    vel[0] = imag(slab2[k_local * N + j]) * vnorm;
                    vel[1] = real(slab2[k_local * N + j]) * vnorm;
                    vel[2] = imag(slab1[k_local * N + j]) * vnorm;
                }

                if (param.qascii) {
                    fmt::print(
                       stdout,
                       "{:d} {:d} {:d} {:f} {:f} {:f} {:f} {:f} {:f} {:f}\n",
                       i,
                       j,
                       k_global,
                       pos[0],
                       pos[1],
                       pos[2],
                       dens,
                       vel[0],
                       vel[1],
                       vel[2]
                    );
                } else {
                    switch (param_icformat) {
                        case OUTPUT_RVDOUBLEZEL: {
                            RVdoubleZelParticle out;
                            out.i        = i;
                            out.j        = j;
                            out.k        = k_global;
                            out.displ[0] = pos[0];
                            out.displ[1] = pos[1];
                            out.displ[2] = pos[2];
                            out.vel[0]   = vel[0];
                            out.vel[1]   = vel[1];
                            out.vel[2]   = vel[2];
                            ((RVdoubleZelParticle *) local_output_tmp)[count] = out;
                            break;
                        }

                        case OUTPUT_RVZEL: {
                            RVZelParticle out;
                            out.i        = i;
                            out.j        = j;
                            out.k        = k_global;
                            out.displ[0] = pos[0];
                            out.displ[1] = pos[1];
                            out.displ[2] = pos[2];
                            out.vel[0]   = vel[0];
                            out.vel[1]   = vel[1];
                            out.vel[2]   = vel[2];
                            ((RVZelParticle *) local_output_tmp)[count] = out;
                            break;
                        }

                        case OUTPUT_ZEL: {
                            ZelParticle out;
                            out.i        = i;
                            out.j        = j;
                            out.k        = k_global;
                            out.displ[0] = pos[0];
                            out.displ[1] = pos[1];
                            out.displ[2] = pos[2];
                            ((ZelParticle *) local_output_tmp)[count] = out;
                            break;
                        }

                        case OUTPUT_ZEL_SIMPLE: {
                            ZelSimpleParticle out;
                            out.displ[0] = pos[0];
                            out.displ[1] = pos[1];
                            out.displ[2] = pos[2];
                            ((ZelSimpleParticle *) local_output_tmp)[count] = out;
                            break;
                        }

                        default:
                            fmt::print(
                               stderr,
                               "Error: unknown ICFormat \"{:s}\" in WriteParticlesSlab_range_from_zslab switch. Aborting.\n",
                               param.ICFormat
                            );
                            exit(1);
                    }
                }

                // Track the global max displacement
                for (int idx = 0; idx < 3; idx++) {
                    max_disp[idx] =
                       fabs(pos[idx]) > fabs(max_disp[idx]) ? pos[idx] : max_disp[idx];
                }
            }

            // Store local density (if requested) and track variance
            if (param.qdensity && local_dens_tmp != NULL) {
                local_dens_tmp[count] = (float)dens;
            }
            thisdensity_variance += dens * dens;

            count++;
        }
    }

#undef ZSLAB_LOCAL

    assert(count == num_particles);

    // Per-rank particle file:
    //   ic_rank{rank}_i{i}_k{k_start_global}_{k_end_global}
    //   where i = Z-slab index (z value), k = X extent
    int k_end_global = k_start_global + k_extent;

    int64_t totsize = 0;

    if (!just_density && local_output_tmp != NULL) {
        fs::path fn_ic = param.output_dir /
                         fmt::format("ic_rank{:d}_i{:d}_k{:d}_{:d}",
                                     rank, i, k_start_global, k_end_global);

        FILE *fp_ic = fopen(fn_ic.c_str(), "wb");
        if (!fp_ic) {
            fmt::print(stderr,
                       "Error: cannot open {} for writing (rank={}, errno={}).\n",
                       fn_ic.string(), rank, errno);
        } else {
            size_t written = fwrite(local_output_tmp, sizeof_outputtype,
                                    (size_t)num_particles, fp_ic);
            if (written != (size_t)num_particles) {
                fmt::print(stderr,
                           "Error: short write to {} (rank={}, wrote {}, expected {}).\n",
                           fn_ic.string(), rank, written, (size_t)num_particles);
            }
            fclose(fp_ic);
            totsize += num_particles * (int64_t)sizeof_outputtype;
        }
    }

    // Optional per-rank density file:
    //   dens_rank{rank}_i{i}_k{k_start_global}_{k_end_global}
    //   where i = Z-slab index (z value), k = X extent
    if (param.qdensity && local_dens_tmp != NULL) {
        fs::path fn_dens = param.output_dir /
                           fmt::format("dens_rank{:d}_i{:d}_k{:d}_{:d}",
                                       rank, i, k_start_global, k_end_global);

        FILE *fp_d = fopen(fn_dens.c_str(), "wb");
        if (!fp_d) {
            fmt::print(stderr,
                       "Warning: cannot open {} for density write (rank={}, errno={}).\n",
                       fn_dens.string(), rank, errno);
        } else {
            size_t written = fwrite(local_dens_tmp, sizeof(float),
                                    (size_t)num_particles, fp_d);
            if (written != (size_t)num_particles) {
                fmt::print(stderr,
                           "Warning: short density write to {} (rank={}, wrote {}, expected {}).\n",
                           fn_dens.string(), rank, written, (size_t)num_particles);
            }
            fclose(fp_d);
            totsize += num_particles * (int64_t)sizeof(float);
        }
    }

    thisouttimer.Stop();
    output_mutex.lock();
    output_bytes_written += totsize;
    outtimer.increment(thisouttimer.timer);
    density_variance += thisdensity_variance;
    output_mutex.unlock();

    // Free local buffer
    if (local_output_tmp != NULL) {
        switch (param_icformat) {
            case OUTPUT_RVDOUBLEZEL:
                delete[] (RVdoubleZelParticle *) local_output_tmp;
                break;
            case OUTPUT_RVZEL:
                delete[] (RVZelParticle *) local_output_tmp;
                break;
            case OUTPUT_ZEL:
                delete[] (ZelParticle *) local_output_tmp;
                break;
            case OUTPUT_ZEL_SIMPLE:
                delete[] (ZelSimpleParticle *) local_output_tmp;
                break;
        }
    }

    if (local_dens_tmp != NULL) {
        delete[] local_dens_tmp;
    }

    return;
}

void SetupOutputDir(Parameters &param) {
    // remove files named ic_* and zeldovich.* from param.output_dir
    if (fs::exists(param.output_dir)) {
        for (const auto& entry : fs::directory_iterator(param.output_dir)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                if ((filename.compare(0, 3, "ic_") == 0) ||
                    (filename.compare(0, 10, "zeldovich.") == 0)) {
                    fs::remove(entry.path());
                }
            }
        }
    }

    fs::create_directories(param.output_dir);
}

// Returns GiB size of allocated buffer
double InitOutputBuffers(Parameters &param) {
    if (param.qdensity != 2) {
        if (param.ICFormat == "RVdoubleZel") {
            param_icformat    = OUTPUT_RVDOUBLEZEL;
            output_tmp        = new RVdoubleZelParticle[param.ppd * param.ppd];
            sizeof_outputtype = sizeof(RVdoubleZelParticle);
        } else if (param.ICFormat == "RVZel") {
            param_icformat    = OUTPUT_RVZEL;
            output_tmp        = new RVZelParticle[param.ppd * param.ppd];
            sizeof_outputtype = sizeof(RVZelParticle);
        } else if (param.ICFormat == "Zeldovich") {
            param_icformat    = OUTPUT_ZEL;
            output_tmp        = new ZelParticle[param.ppd * param.ppd];
            sizeof_outputtype = sizeof(ZelParticle);
        } else if (param.ICFormat == "ZelSimple") {
            param_icformat    = OUTPUT_ZEL_SIMPLE;
            output_tmp        = new ZelSimpleParticle[param.ppd * param.ppd];
            sizeof_outputtype = sizeof(ZelSimpleParticle);
        } else {
            fmt::print(
               stderr, "Error: unknown ICFormat \"{:s}\". Aborting.\n", param.ICFormat
            );
            exit(1);
        }
    } else {
        output_tmp = NULL;
    }

    if (param.qdensity) {
        fs::path path = param.output_dir / fmt::format(fmt::runtime(param.density_filename.string()), param.ppd);

        densfp = fopen(path.c_str(), "wb");
        assert(densfp != NULL);
        densoutput_tmp = new float[param.ppd * param.ppd];
    }

    return sizeof_outputtype * param.ppd * param.ppd / 1024. / 1024. / 1024.;
}

void TeardownOutput() {
    if (output_tmp != NULL) {
        switch (param_icformat) {
            case OUTPUT_RVDOUBLEZEL:
                delete[] (RVdoubleZelParticle *) output_tmp;
                break;

            case OUTPUT_RVZEL:
                delete[] (RVZelParticle *) output_tmp;
                break;

            case OUTPUT_ZEL:
                delete[] (ZelParticle *) output_tmp;
                break;

            case OUTPUT_ZEL_SIMPLE:
                delete[] (ZelSimpleParticle *) output_tmp;
                break;
        }
    }

    if (densoutput_tmp != NULL) {
        delete[] densoutput_tmp;
        fclose(densfp);
    }

    fmt::print(
       stderr,
       "WriteParticlesSlab took {:.3g} sec to write {:.3g} MB ==> {:.3g} MB/sec\n",
       outtimer.Elapsed(),
       output_bytes_written / 1e6,
       output_bytes_written / 1e6 / outtimer.Elapsed()
    );
}
