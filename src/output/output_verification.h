// ====================================================================================
// OUTPUT VERIFICATION UTILITIES
// ====================================================================================
//
// Functions to verify output compatibility with legacy Zeldovich-PLT format
//
// ====================================================================================

#ifndef OUTPUT_VERIFICATION_H
#define OUTPUT_VERIFICATION_H

#include "output_MPI.h"

// ====================================================================================
// VERIFICATION FUNCTIONS
// ====================================================================================

/**
 * Compare two particle files for bit-identical match
 * 
 * @param file1     Path to first file
 * @param file2     Path to second file
 * @param format    Particle format
 * @param verbose   Print detailed differences
 * @return          0 if identical, 1 if different, -1 on error
 */
int CompareParticleFiles(
    const char *file1,
    const char *file2,
    OutputFormat format,
    int verbose
);

/**
 * Verify particle data integrity within a file
 * 
 * Checks:
 * - All i,j,k indices are within valid range [0, N)
 * - Displacements are finite (not NaN or Inf)
 * - Velocities are finite (if present)
 * 
 * @param filename  Path to particle file
 * @param format    Particle format
 * @param N         Expected grid size
 * @param verbose   Print details
 * @return          0 if valid, number of errors found, or -1 on error
 */
int VerifyParticleFile(
    const char *filename,
    OutputFormat format,
    int N,
    int verbose
);

/**
 * Verify axis convention mapping
 * 
 * Checks that particle indices and displacements follow the
 * axis0/axis1/axis2 convention correctly by sampling particles
 * and verifying their positions.
 * 
 * @param filename  Path to particle file
 * @param format    Particle format
 * @param N         Grid size
 * @param verbose   Print details
 * @return          0 if correct, number of errors, or -1 on error
 */
int VerifyAxisConvention(
    const char *filename,
    OutputFormat format,
    int N,
    int verbose
);

/**
 * Compute statistics of particle data in a file
 * 
 * Computes min/max/mean/stddev of displacements and velocities.
 * Useful for sanity checking output.
 * 
 * @param filename  Path to particle file
 * @param format    Particle format
 * @param stats     Output array: [min_displ, max_displ, mean_displ, std_displ, 
 *                                 min_vel, max_vel, mean_vel, std_vel]
 *                  (length 8, or 4 if no velocities)
 * @return          0 on success, -1 on error
 */
int ComputeParticleStats(
    const char *filename,
    OutputFormat format,
    double *stats
);

/**
 * Print summary of all Z-slabs written by a rank
 * 
 * Lists all Z-slab files and their sizes.
 * 
 * @param output_dir  Output directory
 * @param rank        MPI rank number
 * @param N           Grid size
 * @param format      Particle format
 */
void PrintOutputSummary(
    const char *output_dir,
    int rank,
    int N,
    OutputFormat format
);

#endif // OUTPUT_VERIFICATION_H


