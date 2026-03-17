// ====================================================================================
// MPI OUTPUT MODULE - IMPLEMENTATION
// ====================================================================================

#include "output_MPI.h"
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

// ====================================================================================
// UTILITY FUNCTIONS
// ====================================================================================

int ParseOutputFormat(const char *format_str) {
    if (strcmp(format_str, "RVdoubleZel") == 0) {
        return OUTPUT_RVDOUBLEZEL;
    } else if (strcmp(format_str, "RVZel") == 0) {
        return OUTPUT_RVZEL;
    } else if (strcmp(format_str, "Zeldovich") == 0) {
        return OUTPUT_ZEL;
    } else if (strcmp(format_str, "ZelSimple") == 0) {
        return OUTPUT_ZEL_SIMPLE;
    } else {
        return -1;
    }
}

size_t GetParticleSize(OutputFormat format) {
    switch (format) {
        case OUTPUT_RVDOUBLEZEL:
            return sizeof(RVdoubleZelParticle);
        case OUTPUT_RVZEL:
            return sizeof(RVZelParticle);
        case OUTPUT_ZEL:
            return sizeof(ZelParticle);
        case OUTPUT_ZEL_SIMPLE:
            return sizeof(ZelSimpleParticle);
        default:
            return 0;
    }
}

int CreateOutputDirForRank(const char *output_dir, int rank) {
    char rank_dir[1024];
    snprintf(rank_dir, sizeof(rank_dir), "%s/rank_%d", output_dir, rank);
    
    // Create directory (mode 0755)
    if (mkdir(rank_dir, 0755) != 0 && errno != EEXIST) {
        fprintf(stderr, "[Rank %d] Error creating directory %s: %s\n", 
                rank, rank_dir, strerror(errno));
        return -1;
    }
    
    return 0;
}

int WriteMetadata(
    const char *output_dir,
    int N,
    int num_ranks,
    OutputFormat format,
    int write_density
) {
    char metadata_path[1024];
    snprintf(metadata_path, sizeof(metadata_path), "%s/ic_metadata.json", output_dir);
    
    FILE *fp = fopen(metadata_path, "w");
    if (!fp) {
        fprintf(stderr, "Error creating metadata file %s: %s\n", 
                metadata_path, strerror(errno));
        return -1;
    }
    
    const char *format_name;
    switch (format) {
        case OUTPUT_RVDOUBLEZEL: format_name = "RVdoubleZel"; break;
        case OUTPUT_RVZEL:       format_name = "RVZel"; break;
        case OUTPUT_ZEL:         format_name = "Zeldovich"; break;
        case OUTPUT_ZEL_SIMPLE:  format_name = "ZelSimple"; break;
        default:                 format_name = "Unknown"; break;
    }
    
    fprintf(fp, "{\n");
    fprintf(fp, "  \"format_version\": \"MPI_v1_axis012\",\n");
    fprintf(fp, "  \"axis_convention\": {\n");
    fprintf(fp, "    \"axis0\": \"Z (legacy z, perpendicular to slabs)\",\n");
    fprintf(fp, "    \"axis1\": \"Y (legacy y)\",\n");
    fprintf(fp, "    \"axis2\": \"X (legacy x)\"\n");
    fprintf(fp, "  },\n");
    fprintf(fp, "  \"slab_organization\": \"axis0-distributed\",\n");
    fprintf(fp, "  \"legacy_compatible\": true,\n");
    fprintf(fp, "  \"grid_size\": %d,\n", N);
    fprintf(fp, "  \"num_slabs\": %d,\n", N);
    fprintf(fp, "  \"particle_format\": \"%s\",\n", format_name);
    fprintf(fp, "  \"mpi_ranks\": %d,\n", num_ranks);
    fprintf(fp, "  \"displacement_order\": [\"axis0\", \"axis1\", \"axis2\"],\n");
    fprintf(fp, "  \"write_density\": %s,\n", write_density ? "true" : "false");
    fprintf(fp, "  \"per_rank_output\": true,\n");
    fprintf(fp, "  \"file_pattern\": \"rank_{rank}/z{z}_slab_N{N}.bin\"\n");
    fprintf(fp, "}\n");
    
    fclose(fp);
    return 0;
}

// ====================================================================================
// MAIN OUTPUT FUNCTION
// ====================================================================================

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
) {
    // Validate inputs
    if (!displacement_axis0 || !displacement_axis1 || !displacement_axis2) {
        fprintf(stderr, "[Rank %d] Error: displacement arrays cannot be NULL\n", rank);
        return -1;
    }
    
    if (axis0_index < 0 || axis0_index >= N) {
        fprintf(stderr, "[Rank %d] Error: invalid axis0_index %d (N=%d)\n", 
                rank, axis0_index, N);
        return -1;
    }
    
    // Check if velocities are provided (needed for velocity-including formats)
    int has_velocities = (velocity_axis0 != NULL && 
                         velocity_axis1 != NULL && 
                         velocity_axis2 != NULL);
    
    if ((format == OUTPUT_RVDOUBLEZEL || format == OUTPUT_RVZEL) && !has_velocities) {
        fprintf(stderr, "[Rank %d] Warning: format requires velocities but none provided\n", rank);
        // Continue with zero velocities
    }
    
    // Build output filename
    char filename[1024];
    snprintf(filename, sizeof(filename), "%s/rank_%d/z%d_slab_N%d.bin",
             output_dir, rank, axis0_index, N);
    
    // Open file for writing
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "[Rank %d] Error opening file %s: %s\n", 
                rank, filename, strerror(errno));
        return -1;
    }
    
    // Calculate total particles for this rank's contribution to this Z-slab
    int64_t total_particles = (int64_t)N * axis2_count;
    
    // Allocate output buffer based on format
    void *output_buffer = NULL;
    size_t particle_size = GetParticleSize(format);
    
    output_buffer = malloc(total_particles * particle_size);
    if (!output_buffer) {
        fprintf(stderr, "[Rank %d] Error allocating output buffer (%ld particles)\n",
                rank, total_particles);
        fclose(fp);
        return -1;
    }
    
    // Fill buffer with particles
    // Iteration order: for axis1 (Y), for axis2_local (X in local range)
    int64_t particle_idx = 0;
    
    for (int axis1 = 0; axis1 < N; axis1++) {
        for (int axis2_local = 0; axis2_local < axis2_count; axis2_local++) {
            int axis2_global = axis2_start + axis2_local;  // Global X coordinate
            int data_idx = axis1 * axis2_count + axis2_local;  // Index into [axis1][axis2_local] arrays
            
            // Pack particle based on format
            switch (format) {
                case OUTPUT_RVDOUBLEZEL: {
                    RVdoubleZelParticle *particles = (RVdoubleZelParticle *)output_buffer;
                    RVdoubleZelParticle *p = &particles[particle_idx];
                    
                    // Grid indices (legacy order: i=z, j=y, k=x)
                    p->i = axis0_index;      // axis0 (Z)
                    p->j = axis1;            // axis1 (Y)
                    p->k = axis2_global;     // axis2 (X)
                    
                    // Displacements [0]=axis0(Z), [1]=axis1(Y), [2]=axis2(X)
                    p->displ[0] = (double)displacement_axis0[data_idx];
                    p->displ[1] = (double)displacement_axis1[data_idx];
                    p->displ[2] = (double)displacement_axis2[data_idx];
                    
                    // Velocities (if available)
                    if (has_velocities) {
                        p->vel[0] = (double)velocity_axis0[data_idx];
                        p->vel[1] = (double)velocity_axis1[data_idx];
                        p->vel[2] = (double)velocity_axis2[data_idx];
                    } else {
                        p->vel[0] = 0.0;
                        p->vel[1] = 0.0;
                        p->vel[2] = 0.0;
                    }
                    break;
                }
                
                case OUTPUT_RVZEL: {
                    RVZelParticle *particles = (RVZelParticle *)output_buffer;
                    RVZelParticle *p = &particles[particle_idx];
                    
                    p->i = axis0_index;
                    p->j = axis1;
                    p->k = axis2_global;
                    
                    p->displ[0] = (float)displacement_axis0[data_idx];
                    p->displ[1] = (float)displacement_axis1[data_idx];
                    p->displ[2] = (float)displacement_axis2[data_idx];
                    
                    if (has_velocities) {
                        p->vel[0] = (float)velocity_axis0[data_idx];
                        p->vel[1] = (float)velocity_axis1[data_idx];
                        p->vel[2] = (float)velocity_axis2[data_idx];
                    } else {
                        p->vel[0] = 0.0f;
                        p->vel[1] = 0.0f;
                        p->vel[2] = 0.0f;
                    }
                    break;
                }
                
                case OUTPUT_ZEL: {
                    ZelParticle *particles = (ZelParticle *)output_buffer;
                    ZelParticle *p = &particles[particle_idx];
                    
                    p->i = axis0_index;
                    p->j = axis1;
                    p->k = axis2_global;
                    
                    p->displ[0] = (float)displacement_axis0[data_idx];
                    p->displ[1] = (float)displacement_axis1[data_idx];
                    p->displ[2] = (float)displacement_axis2[data_idx];
                    break;
                }
                
                case OUTPUT_ZEL_SIMPLE: {
                    ZelSimpleParticle *particles = (ZelSimpleParticle *)output_buffer;
                    ZelSimpleParticle *p = &particles[particle_idx];
                    
                    p->displ[0] = (float)displacement_axis0[data_idx];
                    p->displ[1] = (float)displacement_axis1[data_idx];
                    p->displ[2] = (float)displacement_axis2[data_idx];
                    break;
                }
                
                default:
                    fprintf(stderr, "[Rank %d] Error: unknown output format %d\n", rank, format);
                    free(output_buffer);
                    fclose(fp);
                    return -1;
            }
            
            particle_idx++;
        }
    }
    
    // Write buffer to file
    size_t items_written = fwrite(output_buffer, particle_size, total_particles, fp);
    if (items_written != total_particles) {
        fprintf(stderr, "[Rank %d] Error writing particles: wrote %zu of %ld\n",
                rank, items_written, total_particles);
        free(output_buffer);
        fclose(fp);
        return -1;
    }
    
    // Clean up
    free(output_buffer);
    fclose(fp);
    
    // Optionally write density file
    if (write_density && density_data != NULL) {
        char density_filename[1024];
        snprintf(density_filename, sizeof(density_filename), 
                "%s/rank_%d/z%d_density_N%d.bin",
                output_dir, rank, axis0_index, N);
        
        FILE *dens_fp = fopen(density_filename, "wb");
        if (!dens_fp) {
            fprintf(stderr, "[Rank %d] Warning: could not open density file %s\n",
                    rank, density_filename);
            return 0;  // Don't fail if density write fails
        }
        
        // Write density as float array [axis1][axis2_local]
        float *density_float = (float *)malloc(total_particles * sizeof(float));
        if (density_float) {
            for (int64_t i = 0; i < total_particles; i++) {
                density_float[i] = (float)density_data[i];
            }
            fwrite(density_float, sizeof(float), total_particles, dens_fp);
            free(density_float);
        }
        
        fclose(dens_fp);
    }
    
    return 0;
}

