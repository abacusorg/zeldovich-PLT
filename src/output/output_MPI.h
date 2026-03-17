// ====================================================================================
// MPI OUTPUT MODULE - HEADER
// ====================================================================================
// 
// This module provides output functions for the MPI-parallelized IC generation
// code, with backward compatibility to the Zeldovich-PLT output format.
//
// KEY FEATURES:
//   - Clear axis0/axis1/axis2 naming convention (no legacy confusion)
//   - Z-distributed slab output (axis0-distributed)
//   - Each MPI rank writes its X-range (axis2-range) to separate files
//   - Backward compatible with existing Abacus loadIC.cpp
//   - Supports all legacy particle formats (RVZel, RVdoubleZel, Zeldovich, ZelSimple)
//
// ====================================================================================

#ifndef OUTPUT_MPI_H
#define OUTPUT_MPI_H

#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include "../precision.h"

// ====================================================================================
// AXIS CONVENTION DOCUMENTATION
// ====================================================================================
//
// NEW NOTATION (used throughout this code):
//   axis0 = Z coordinate (slab index, perpendicular to slabs)
//   axis1 = Y coordinate 
//   axis2 = X coordinate
//
// LEGACY ZELDOVICH MAPPING:
//   Legacy variable "z" -> axis0 (Z)
//   Legacy variable "y" -> axis1 (Y)
//   Legacy variable "x" -> axis2 (X)
//   Legacy out.displ[0] -> axis0 displacement (Z-displacement)
//   Legacy out.displ[1] -> axis1 displacement (Y-displacement)
//   Legacy out.displ[2] -> axis2 displacement (X-displacement)
//
// SLAB ORGANIZATION:
//   - Write axis0-distributed slabs (Z-slabs)
//   - Each file contains: all axis2 (X) x all axis1 (Y) for one axis0 (Z)
//   - File naming: rank_{rank}/z{z}_slab_N{N}.bin or ic_{z}.bin (after reassembly)
//
// MPI DATA DISTRIBUTION (after 3D FFT in MPI code):
//   - Data is X-distributed (axis2-distributed)
//   - Each rank owns: axis2_range x all_axis1 x all_axis0
//   - For each Z-slab: each rank writes its axis2_range for all axis1
//
// ITERATION ORDER (for single axis0 value):
//   for (axis1 = 0; axis1 < N; axis1++)     // All Y values
//     for (axis2 = 0; axis2 < N; axis2++)   // All X values
//       Particle at (axis2, axis1, axis0)
//       Write in order: out.i=axis0, out.j=axis1, out.k=axis2
//       Displacement order: out.displ[0]=axis0, [1]=axis1, [2]=axis2
//
// ====================================================================================

// ====================================================================================
// PARTICLE STRUCT DEFINITIONS (copied from Zeldovich-PLT output.h)
// ====================================================================================

// RVdoubleZel: positions (i,j,k) + displacements + velocities (double precision)
typedef struct {
    int32_t i, j, k;        // Grid indices (axis0, axis1, axis2)
    double displ[3];        // Displacements [axis0, axis1, axis2]
    double vel[3];          // Velocities [axis0, axis1, axis2]
} RVdoubleZelParticle;

// RVZel: positions (i,j,k) + displacements + velocities (float precision)
typedef struct {
    int32_t i, j, k;        // Grid indices (axis0, axis1, axis2)
    float displ[3];         // Displacements [axis0, axis1, axis2]
    float vel[3];           // Velocities [axis0, axis1, axis2]
} RVZelParticle;

// Zeldovich: positions (i,j,k) + displacements only (float precision)
typedef struct {
    int32_t i, j, k;        // Grid indices (axis0, axis1, axis2)
    float displ[3];         // Displacements [axis0, axis1, axis2]
} ZelParticle;

// ZelSimple: displacements only, no positions (float precision)
typedef struct {
    float displ[3];         // Displacements [axis0, axis1, axis2]
} ZelSimpleParticle;

// ====================================================================================
// OUTPUT FORMAT ENUM
// ====================================================================================

typedef enum {
    OUTPUT_RVDOUBLEZEL,     // RVdoubleZelParticle (double precision)
    OUTPUT_RVZEL,           // RVZelParticle (float precision)
    OUTPUT_ZEL,             // ZelParticle (positions + displacements)
    OUTPUT_ZEL_SIMPLE       // ZelSimpleParticle (displacements only)
} OutputFormat;

// ====================================================================================
// MAIN OUTPUT FUNCTION
// ====================================================================================

/**
 * Write a single Z-slab (axis0-distributed slab) for this MPI rank's X-range
 * 
 * This function writes particles for:
 *   - One axis0 (Z) value
 *   - All axis1 (Y) values [0, N)
 *   - This rank's axis2 (X) range [axis2_start, axis2_start + axis2_count)
 * 
 * Output file: {output_dir}/rank_{rank}/z{axis0}_slab_N{N}.bin
 * 
 * Data layout:
 *   Input arrays are [axis1][axis2_local] (row-major)
 *   Output iteration: for axis1, for axis2_local
 *   Output struct: RVZelParticle with (i=axis0, j=axis1, k=axis2_global)
 * 
 * @param axis0_index       Z value (slab index) [0, N)
 * @param density_data      Density field: [axis1][axis2_local] or NULL if not writing density
 * @param displacement_axis0 Z-displacement: [axis1][axis2_local]
 * @param displacement_axis1 Y-displacement: [axis1][axis2_local]
 * @param displacement_axis2 X-displacement: [axis1][axis2_local]
 * @param velocity_axis0    Z-velocity: [axis1][axis2_local] or NULL if no velocities
 * @param velocity_axis1    Y-velocity: [axis1][axis2_local] or NULL if no velocities
 * @param velocity_axis2    X-velocity: [axis1][axis2_local] or NULL if no velocities
 * @param N                 Grid size (total particles per dimension)
 * @param axis2_start       Starting X index for this rank
 * @param axis2_count       Number of X values this rank owns
 * @param output_dir        Output directory path
 * @param format            Particle format (OUTPUT_RVZEL, etc.)
 * @param rank              MPI rank number
 * @param write_density     Whether to write separate density file
 * @return                  0 on success, -1 on error
 */
int WriteParticlesSlabMPI(
    int axis0_index,
    const real_t *density_data,
    const real_t *displacement_axis0,
    const real_t *displacement_axis1,
    const real_t *displacement_axis2,
    const real_t *velocity_axis0,
    const real_t *velocity_axis1,
    const real_t *velocity_axis2,
    int N,
    int axis2_start,
    int axis2_count,
    const char *output_dir,
    OutputFormat format,
    int rank,
    int write_density
);

// ====================================================================================
// UTILITY FUNCTIONS
// ====================================================================================

/**
 * Parse output format string to OutputFormat enum
 * 
 * @param format_str  Format string ("RVdoubleZel", "RVZel", "Zeldovich", "ZelSimple")
 * @return            OutputFormat enum value, or -1 if invalid
 */
int ParseOutputFormat(const char *format_str);

/**
 * Get size of particle struct for given format
 * 
 * @param format  OutputFormat enum value
 * @return        Size in bytes of particle struct
 */
size_t GetParticleSize(OutputFormat format);

/**
 * Create output directory structure for this rank
 * Creates: {output_dir}/rank_{rank}/
 * 
 * @param output_dir  Base output directory
 * @param rank        MPI rank number
 * @return            0 on success, -1 on error
 */
int CreateOutputDirForRank(const char *output_dir, int rank);

/**
 * Write metadata file with axis convention and grid information
 * Creates: {output_dir}/ic_metadata.json
 * 
 * @param output_dir      Base output directory
 * @param N               Grid size
 * @param num_ranks       Total number of MPI ranks
 * @param format          Particle format
 * @param write_density   Whether density files are written
 * @return                0 on success, -1 on error
 */
int WriteMetadata(
    const char *output_dir,
    int N,
    int num_ranks,
    OutputFormat format,
    int write_density
);

#endif // OUTPUT_MPI_H

