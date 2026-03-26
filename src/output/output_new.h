#ifndef OUTPUT_NEW_H
#define OUTPUT_NEW_H

// Header for output_new.cpp functions
// Provides WriteParticlesSlab_range (overloaded for both [array][x][y] and [y][x] layouts),
// WriteParticlesSlab_range_from_zslab (legacy, now equivalent to WriteParticlesSlab_range),
// and WriteParticlesSlab_new (Path 2: full-range reassembly)

#include <stdio.h>
#include "../precision.h"
// Undefine MAX_PPD macro from config.h to avoid conflict with zeldovich-PLT's const definition
#ifdef MAX_PPD
#undef MAX_PPD
#endif
// Include zeldovich-PLT headers to get Complx and Parameters definitions
#include <zeldovich.h>
#include <parameters.h>
#include "output_types.h"

// WriteParticlesSlab_range - C++ function overloads (same name, different signatures)
// The compiler selects the appropriate version based on the arguments provided.

// Overload 1: [array][x][y] layout (ZSLAB format) - NO TRANSPOSE NEEDED
// This is the preferred version for new code as it works directly with main.cpp's ZSLAB format
void WriteParticlesSlab_range(
   int rank,
   int i,                  // i = Z index (legacy z)
   int k_start_global,     // global X start for this rank's extent
   int k_extent,           // number of X values for this rank
   Complx *slab_data,      // Data in [array][x_local][y] layout (ZSLAB format)
   int N,                  // Grid size (ppd)
   int narray,             // Number of arrays (typically 4)
   Parameters &param
);

// Overload 2: [y][x] layout (JK format) - for backward compatibility
// This version accepts separate pointers for each array in [y][x] layout
void WriteParticlesSlab_range(
   int rank,
   int i,                  // i = Z index (legacy z)
   int k_start_global,     // global X start for this rank's extent
   int k_extent,           // number of X values for this rank
   Complx *slab1,          // Array 0 in [y][x] layout
   Complx *slab2,          // Array 1 in [y][x] layout
   Complx *slab3,          // Array 2 in [y][x] layout
   Complx *slab4,          // Array 3 in [y][x] layout
   Parameters &param
);

// WriteParticlesSlab_range_from_zslab - Legacy function (now equivalent to WriteParticlesSlab_range)
// This function is kept for backward compatibility but is functionally identical to
// WriteParticlesSlab_range overload 1. New code should use WriteParticlesSlab_range instead.
void WriteParticlesSlab_range_from_zslab(
   int rank,
   int i,                  // i = Z index (legacy z)
   int k_start_global,     // global X start for this rank's extent
   int k_extent,           // number of X values for this rank (x_count)
   Complx *slab_data,      // Data in [array][x_local][y] layout (ZSLAB format)
   int N,                  // Grid size (ppd)
   int narray,             // Number of arrays (typically 4)
   Parameters &param
);

// Write full-range particle ICs (all X, all Y) for one i-slab
// Used by Path 2 reassembly tool
void WriteParticlesSlab_new(
   FILE *output,
   int i,                  // i = Z index (legacy z)
   Complx *slab1,          // Array 0 in [y][x] layout
   Complx *slab2,          // Array 1 in [y][x] layout
   Complx *slab3,          // Array 2 in [y][x] layout
   Complx *slab4,          // Array 3 in [y][x] layout
   Parameters &param
);

// Setup and teardown functions
void SetupOutputDir(Parameters &param);
double InitOutputBuffers(Parameters &param);
void TeardownOutput();

// ====================================================================================
// MODE 3: One file per x-slab (ic2D_{xslab}_z{rz}.bin)
// ====================================================================================
//
// CPD-aligned; layout must match Abacus RVZel_2D reader.
// Each file contains [z0 segment][z1 segment]... for one x-slab (sequential, no offset).
// Indices: i, j, k are all global (same convention as grid_x==1 / zeldovich).

// Write one x-slab segment for one z to that slab's file.
void AppendSlabZSegment(
    FILE *fp,                 // This slab's particle file (caller owns)
    FILE *fp_dens,            // This slab's density file, or NULL
    int slab_s,               // CPD slab index (for firstx/lastx)
    int cpd,                  // Cells per dimension
    int z,                    // Global z index for this segment
    int k_start_global,
    int k_extent,
    fftw_complex_t *slab_data,
    int N,
    int narray,
    Parameters &param
);

// ====================================================================================
// MODE 3 (grid_x==1): One file per z-group, matching zeldovich output format
// ====================================================================================
//
// When grid_x==1 each rank owns all N x-values.  Each call writes one full
// N×N z-plane to the z-group file.  Indices are all global (i=z, j=y, k=x).
// Particle ordering: y-outer, x-inner (matching zeldovich WriteParticlesSlab).

void AppendZSlabFull(
    FILE *fp,                 // z-group file (caller owns, opened for append)
    FILE *fp_dens,            // density file, or NULL
    int z,                    // global z index
    int k_start_global,       // == 0 when grid_x==1
    int k_extent,             // == N when grid_x==1
    fftw_complex_t *slab_data,
    int N,
    int narray,
    Parameters &param
);

#endif 

