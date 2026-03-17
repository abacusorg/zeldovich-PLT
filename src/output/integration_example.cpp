// ====================================================================================
// EXAMPLE INTEGRATION: Using output_MPI module in main MPI code
// ====================================================================================
//
// This file demonstrates how to integrate the output_MPI module into
// the main MPI IC generation code after the 3D FFT is complete.
//
// ====================================================================================

#include "output/output_MPI.h"
#include <string.h>  // for memcpy
#include <stdlib.h>
#include <mpi.h>

/**
 * Example integration function
 * 
 * This function is called after 3D FFT is complete and all data is in real space.
 * Data layout: [axis0][axis1][axis2_local] (X-distributed)
 * 
 * @param fft_output       Array of 7 3D arrays: [D, F, G, H, Fv, Gv, Hv]
 *                         Each array is [axis0][axis1][axis2_local]
 * @param N                Grid size
 * @param axis2_start      Starting X index for this rank
 * @param axis2_count      Number of X values this rank owns
 * @param output_dir       Output directory path
 * @param format_str       Particle format string ("RVZel", "RVdoubleZel", etc.)
 * @param write_density    Whether to write density files
 * @param has_velocities   Whether velocity arrays are present (PLT mode)
 * @param rank             MPI rank
 * @param num_ranks        Total number of MPI ranks
 * @param comm             MPI communicator
 */
void WriteOutputAfterFFT(
    real_t **fft_output,
    int N,
    int axis2_start,
    int axis2_count,
    const char *output_dir,
    const char *format_str,
    int write_density,
    int has_velocities,
    int rank,
    int num_ranks,
    MPI_Comm comm
) {
    // Parse output format
    OutputFormat format = (OutputFormat)ParseOutputFormat(format_str);
    if (format < 0) {
        if (rank == 0) {
            fprintf(stderr, "Error: Unknown output format '%s'\n", format_str);
        }
        return;
    }
    
    // Create output directory for this rank
    if (CreateOutputDirForRank(output_dir, rank) != 0) {
        fprintf(stderr, "[Rank %d] Failed to create output directory\n", rank);
        return;
    }
    
    // Write metadata (rank 0 only)
    if (rank == 0) {
        WriteMetadata(output_dir, N, num_ranks, format, write_density);
    }
    
    // Allocate buffers for 2D Z-slices [axis1][axis2_local]
    int64_t slice_size = (int64_t)N * axis2_count;
    
    real_t *density_slice = NULL;
    real_t *displ_axis0_slice = malloc(slice_size * sizeof(real_t));
    real_t *displ_axis1_slice = malloc(slice_size * sizeof(real_t));
    real_t *displ_axis2_slice = malloc(slice_size * sizeof(real_t));
    real_t *vel_axis0_slice = NULL;
    real_t *vel_axis1_slice = NULL;
    real_t *vel_axis2_slice = NULL;
    
    if (write_density) {
        density_slice = malloc(slice_size * sizeof(real_t));
    }
    
    if (has_velocities) {
        vel_axis0_slice = malloc(slice_size * sizeof(real_t));
        vel_axis1_slice = malloc(slice_size * sizeof(real_t));
        vel_axis2_slice = malloc(slice_size * sizeof(real_t));
    }
    
    // Check allocation
    if (!displ_axis0_slice || !displ_axis1_slice || !displ_axis2_slice) {
        fprintf(stderr, "[Rank %d] Failed to allocate slice buffers\n", rank);
        goto cleanup;
    }
    
    if (write_density && !density_slice) {
        fprintf(stderr, "[Rank %d] Failed to allocate density slice buffer\n", rank);
        goto cleanup;
    }
    
    if (has_velocities && (!vel_axis0_slice || !vel_axis1_slice || !vel_axis2_slice)) {
        fprintf(stderr, "[Rank %d] Failed to allocate velocity slice buffers\n", rank);
        goto cleanup;
    }
    
    // Progress reporting
    if (rank == 0) {
        printf("\n[Output] Writing Z-slabs to disk...\n");
        printf("[Output] Format: %s\n", format_str);
        printf("[Output] Grid size: %d\n", N);
        printf("[Output] Writing density: %s\n", write_density ? "yes" : "no");
        printf("[Output] Has velocities: %s\n", has_velocities ? "yes" : "no");
    }
    
    // Loop over all Z values (axis0)
    for (int axis0 = 0; axis0 < N; axis0++) {
        // Progress reporting (every 100 slabs)
        if (rank == 0 && axis0 % 100 == 0) {
            printf("[Output] Writing Z-slab %d / %d (%.1f%%)\n", 
                   axis0, N, 100.0 * axis0 / N);
        }
        
        // Data is already in correct format [array][y][x] for this Z-slab
        // fft_output array indices: [0]=D, [1]=F, [2]=G, [3]=H, [4]=Fv, [5]=Gv, [6]=Hv
        // Note: F, G, H correspond to axis0, axis1, axis2 displacements respectively
        // No extraction needed - data is already in [axis1][axis2_local] format
        
        if (write_density && fft_output[0] != NULL) {
            memcpy(density_slice, fft_output[0], slice_size * sizeof(real_t));
        }
        
        // Copy displacements: F->axis0(Z), G->axis1(Y), H->axis2(X)
        memcpy(displ_axis0_slice, fft_output[1], slice_size * sizeof(real_t));  // F -> Z
        memcpy(displ_axis1_slice, fft_output[2], slice_size * sizeof(real_t));  // G -> Y
        memcpy(displ_axis2_slice, fft_output[3], slice_size * sizeof(real_t));  // H -> X
        
        // Copy velocities if present: Fv->axis0(Z), Gv->axis1(Y), Hv->axis2(X)
        if (has_velocities) {
            memcpy(vel_axis0_slice, fft_output[4], slice_size * sizeof(real_t));  // Fv -> Z
            memcpy(vel_axis1_slice, fft_output[5], slice_size * sizeof(real_t));  // Gv -> Y
            memcpy(vel_axis2_slice, fft_output[6], slice_size * sizeof(real_t));  // Hv -> X
        }
        
        // Write Z-slab for this rank's X-range
        int status = WriteParticlesSlabMPI(
            axis0,                      // Z index
            density_slice,              // Density (or NULL)
            displ_axis0_slice,          // Z-displacement
            displ_axis1_slice,          // Y-displacement
            displ_axis2_slice,          // X-displacement
            vel_axis0_slice,            // Z-velocity (or NULL)
            vel_axis1_slice,            // Y-velocity (or NULL)
            vel_axis2_slice,            // X-velocity (or NULL)
            N,
            axis2_start,                // Starting X index for this rank
            axis2_count,                // Number of X values this rank owns
            output_dir,
            format,
            rank,
            write_density
        );
        
        if (status != 0) {
            fprintf(stderr, "[Rank %d] Error writing Z-slab %d\n", rank, axis0);
            goto cleanup;
        }
    }
    
    if (rank == 0) {
        printf("[Output] All Z-slabs written successfully!\n");
        printf("[Output] Files location: %s/rank_*/z*_slab_N%d.bin\n", output_dir, N);
    }
    
cleanup:
    // Free allocated buffers
    if (density_slice) free(density_slice);
    if (displ_axis0_slice) free(displ_axis0_slice);
    if (displ_axis1_slice) free(displ_axis1_slice);
    if (displ_axis2_slice) free(displ_axis2_slice);
    if (vel_axis0_slice) free(vel_axis0_slice);
    if (vel_axis1_slice) free(vel_axis1_slice);
    if (vel_axis2_slice) free(vel_axis2_slice);
}

// ====================================================================================
// INTEGRATION NOTES
// ====================================================================================
//
// To integrate into main.cpp:
//
// 1. After 3D FFT is complete and data is in real space, call WriteOutputAfterFFT()
//
// 2. Pass the FFT output arrays in this order:
//    fft_output[0] = D (density)
//    fft_output[1] = F (axis0/Z-displacement)
//    fft_output[2] = G (axis1/Y-displacement)
//    fft_output[3] = H (axis2/X-displacement)
//    fft_output[4] = Fv (axis0/Z-velocity) [if PLT]
//    fft_output[5] = Gv (axis1/Y-velocity) [if PLT]
//    fft_output[6] = Hv (axis2/X-velocity) [if PLT]
//
// 3. Data layout for each array must be: [axis1][axis2_local] (already in correct format)
//    where axis2_local is the X-distributed local range
//
// 4. Example call in main.cpp:
//    
//    real_t *fft_arrays[7] = {
//        global_output_D,
//        global_output_F,
//        global_output_G,
//        global_output_H,
//        global_output_Fv,  // or NULL if not PLT
//        global_output_Gv,  // or NULL if not PLT
//        global_output_Hv   // or NULL if not PLT
//    };
//    
//    WriteOutputAfterFFT(
//        fft_arrays,
//        N,
//        x_start,           // axis2_start
//        local_x_count,     // axis2_count
//        output_dir,
//        "RVZel",           // or parameter from config
//        qdensity,          // write_density flag
//        qPLT,              // has_velocities flag
//        rank,
//        num_ranks,
//        MPI_COMM_WORLD
//    );
//
// ====================================================================================


