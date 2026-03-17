// ====================================================================================
// OUTPUT VERIFICATION UTILITIES - IMPLEMENTATION
// ====================================================================================

#include "output_verification.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>

// Helper: Check if value is finite
static int is_finite(double val) {
    return !isnan(val) && !isinf(val);
}

int CompareParticleFiles(
    const char *file1,
    const char *file2,
    OutputFormat format,
    int verbose
) {
    FILE *fp1 = fopen(file1, "rb");
    FILE *fp2 = fopen(file2, "rb");
    
    if (!fp1 || !fp2) {
        if (verbose) {
            fprintf(stderr, "Error opening files for comparison\n");
            if (!fp1) fprintf(stderr, "  Cannot open: %s\n", file1);
            if (!fp2) fprintf(stderr, "  Cannot open: %s\n", file2);
        }
        if (fp1) fclose(fp1);
        if (fp2) fclose(fp2);
        return -1;
    }
    
    // Get file sizes
    fseek(fp1, 0, SEEK_END);
    long size1 = ftell(fp1);
    fseek(fp1, 0, SEEK_SET);
    
    fseek(fp2, 0, SEEK_END);
    long size2 = ftell(fp2);
    fseek(fp2, 0, SEEK_SET);
    
    if (size1 != size2) {
        if (verbose) {
            fprintf(stderr, "File sizes differ: %ld vs %ld bytes\n", size1, size2);
        }
        fclose(fp1);
        fclose(fp2);
        return 1;
    }
    
    // Compare byte-by-byte
    size_t particle_size = GetParticleSize(format);
    size_t num_particles = size1 / particle_size;
    
    int num_differences = 0;
    
    for (size_t i = 0; i < num_particles; i++) {
        char buf1[256], buf2[256];
        
        fread(buf1, particle_size, 1, fp1);
        fread(buf2, particle_size, 1, fp2);
        
        if (memcmp(buf1, buf2, particle_size) != 0) {
            num_differences++;
            if (verbose && num_differences <= 10) {
                fprintf(stderr, "Difference at particle %zu\n", i);
            }
        }
    }
    
    fclose(fp1);
    fclose(fp2);
    
    if (verbose && num_differences > 0) {
        fprintf(stderr, "Total differences: %d / %zu particles\n", 
                num_differences, num_particles);
    }
    
    return num_differences > 0 ? 1 : 0;
}

int VerifyParticleFile(
    const char *filename,
    OutputFormat format,
    int N,
    int verbose
) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        if (verbose) {
            fprintf(stderr, "Error opening file: %s\n", filename);
        }
        return -1;
    }
    
    size_t particle_size = GetParticleSize(format);
    
    // Get number of particles
    fseek(fp, 0, SEEK_END);
    long file_size = ftell(fp);
    size_t num_particles = file_size / particle_size;
    fseek(fp, 0, SEEK_SET);
    
    if (verbose) {
        printf("Verifying file: %s\n", filename);
        printf("  Particles: %zu\n", num_particles);
        printf("  Format: %d\n", format);
    }
    
    int num_errors = 0;
    
    // Read and verify each particle
    for (size_t p = 0; p < num_particles; p++) {
        // Read particle based on format
        int i = -1, j = -1, k = -1;
        float displ[3] = {0, 0, 0};
        float vel[3] = {0, 0, 0};
        
        switch (format) {
            case OUTPUT_RVDOUBLEZEL: {
                RVdoubleZelParticle particle;
                fread(&particle, sizeof(particle), 1, fp);
                i = particle.i;
                j = particle.j;
                k = particle.k;
                displ[0] = particle.displ[0];
                displ[1] = particle.displ[1];
                displ[2] = particle.displ[2];
                vel[0] = particle.vel[0];
                vel[1] = particle.vel[1];
                vel[2] = particle.vel[2];
                break;
            }
            case OUTPUT_RVZEL: {
                RVZelParticle particle;
                fread(&particle, sizeof(particle), 1, fp);
                i = particle.i;
                j = particle.j;
                k = particle.k;
                displ[0] = particle.displ[0];
                displ[1] = particle.displ[1];
                displ[2] = particle.displ[2];
                vel[0] = particle.vel[0];
                vel[1] = particle.vel[1];
                vel[2] = particle.vel[2];
                break;
            }
            case OUTPUT_ZEL: {
                ZelParticle particle;
                fread(&particle, sizeof(particle), 1, fp);
                i = particle.i;
                j = particle.j;
                k = particle.k;
                displ[0] = particle.displ[0];
                displ[1] = particle.displ[1];
                displ[2] = particle.displ[2];
                break;
            }
            case OUTPUT_ZEL_SIMPLE: {
                ZelSimpleParticle particle;
                fread(&particle, sizeof(particle), 1, fp);
                displ[0] = particle.displ[0];
                displ[1] = particle.displ[1];
                displ[2] = particle.displ[2];
                break;
            }
        }
        
        // Verify indices (if present)
        if (format != OUTPUT_ZEL_SIMPLE) {
            if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N) {
                num_errors++;
                if (verbose && num_errors <= 10) {
                    fprintf(stderr, "  Particle %zu: invalid indices (%d, %d, %d)\n",
                            p, i, j, k);
                }
            }
        }
        
        // Verify displacements are finite
        for (int d = 0; d < 3; d++) {
            if (!is_finite(displ[d])) {
                num_errors++;
                if (verbose && num_errors <= 10) {
                    fprintf(stderr, "  Particle %zu: non-finite displacement[%d] = %f\n",
                            p, d, displ[d]);
                }
            }
        }
        
        // Verify velocities are finite (if present)
        if (format == OUTPUT_RVDOUBLEZEL || format == OUTPUT_RVZEL) {
            for (int d = 0; d < 3; d++) {
                if (!is_finite(vel[d])) {
                    num_errors++;
                    if (verbose && num_errors <= 10) {
                        fprintf(stderr, "  Particle %zu: non-finite velocity[%d] = %f\n",
                                p, d, vel[d]);
                    }
                }
            }
        }
    }
    
    fclose(fp);
    
    if (verbose) {
        if (num_errors == 0) {
            printf("  Verification: PASS (no errors)\n");
        } else {
            printf("  Verification: FAIL (%d errors)\n", num_errors);
        }
    }
    
    return num_errors;
}

int VerifyAxisConvention(
    const char *filename,
    OutputFormat format,
    int N,
    int verbose
) {
    if (format == OUTPUT_ZEL_SIMPLE) {
        if (verbose) {
            printf("Skipping axis convention check for ZelSimple format (no indices)\n");
        }
        return 0;
    }
    
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        if (verbose) {
            fprintf(stderr, "Error opening file: %s\n", filename);
        }
        return -1;
    }
    
    // Read first particle to check axis convention
    int i = -1, j = -1, k = -1;
    
    switch (format) {
        case OUTPUT_RVDOUBLEZEL: {
            RVdoubleZelParticle particle;
            fread(&particle, sizeof(particle), 1, fp);
            i = particle.i;
            j = particle.j;
            k = particle.k;
            break;
        }
        case OUTPUT_RVZEL: {
            RVZelParticle particle;
            fread(&particle, sizeof(particle), 1, fp);
            i = particle.i;
            j = particle.j;
            k = particle.k;
            break;
        }
        case OUTPUT_ZEL: {
            ZelParticle particle;
            fread(&particle, sizeof(particle), 1, fp);
            i = particle.i;
            j = particle.j;
            k = particle.k;
            break;
        }
        default:
            break;
    }
    
    fclose(fp);
    
    if (verbose) {
        printf("Axis convention check:\n");
        printf("  First particle: i=%d (axis0/Z), j=%d (axis1/Y), k=%d (axis2/X)\n", i, j, k);
        printf("  Expected: i=axis0 (Z-slab index), j in [0,%d), k in [0,%d)\n", N, N);
    }
    
    // Basic sanity check
    int num_errors = 0;
    if (i < 0 || i >= N) {
        if (verbose) fprintf(stderr, "  ERROR: i (axis0) out of range\n");
        num_errors++;
    }
    if (j < 0 || j >= N) {
        if (verbose) fprintf(stderr, "  ERROR: j (axis1) out of range\n");
        num_errors++;
    }
    if (k < 0 || k >= N) {
        if (verbose) fprintf(stderr, "  ERROR: k (axis2) out of range\n");
        num_errors++;
    }
    
    return num_errors;
}

int ComputeParticleStats(
    const char *filename,
    OutputFormat format,
    double *stats
) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        return -1;
    }
    
    size_t particle_size = GetParticleSize(format);
    fseek(fp, 0, SEEK_END);
    size_t num_particles = ftell(fp) / particle_size;
    fseek(fp, 0, SEEK_SET);
    
    double displ_sum[3] = {0, 0, 0};
    double displ_sum2[3] = {0, 0, 0};
    double displ_min[3] = {1e30, 1e30, 1e30};
    double displ_max[3] = {-1e30, -1e30, -1e30};
    
    double vel_sum[3] = {0, 0, 0};
    double vel_sum2[3] = {0, 0, 0};
    double vel_min[3] = {1e30, 1e30, 1e30};
    double vel_max[3] = {-1e30, -1e30, -1e30};
    
    int has_velocities = (format == OUTPUT_RVDOUBLEZEL || format == OUTPUT_RVZEL);
    
    for (size_t p = 0; p < num_particles; p++) {
        float displ[3] = {0, 0, 0};
        float vel[3] = {0, 0, 0};
        
        switch (format) {
            case OUTPUT_RVDOUBLEZEL: {
                RVdoubleZelParticle particle;
                fread(&particle, sizeof(particle), 1, fp);
                for (int d = 0; d < 3; d++) {
                    displ[d] = particle.displ[d];
                    vel[d] = particle.vel[d];
                }
                break;
            }
            case OUTPUT_RVZEL: {
                RVZelParticle particle;
                fread(&particle, sizeof(particle), 1, fp);
                for (int d = 0; d < 3; d++) {
                    displ[d] = particle.displ[d];
                    vel[d] = particle.vel[d];
                }
                break;
            }
            case OUTPUT_ZEL: {
                ZelParticle particle;
                fread(&particle, sizeof(particle), 1, fp);
                for (int d = 0; d < 3; d++) {
                    displ[d] = particle.displ[d];
                }
                break;
            }
            case OUTPUT_ZEL_SIMPLE: {
                ZelSimpleParticle particle;
                fread(&particle, sizeof(particle), 1, fp);
                for (int d = 0; d < 3; d++) {
                    displ[d] = particle.displ[d];
                }
                break;
            }
        }
        
        for (int d = 0; d < 3; d++) {
            displ_sum[d] += displ[d];
            displ_sum2[d] += displ[d] * displ[d];
            if (displ[d] < displ_min[d]) displ_min[d] = displ[d];
            if (displ[d] > displ_max[d]) displ_max[d] = displ[d];
            
            if (has_velocities) {
                vel_sum[d] += vel[d];
                vel_sum2[d] += vel[d] * vel[d];
                if (vel[d] < vel_min[d]) vel_min[d] = vel[d];
                if (vel[d] > vel_max[d]) vel_max[d] = vel[d];
            }
        }
    }
    
    fclose(fp);
    
    // Compute stats (average over 3 dimensions)
    double displ_mean = 0, displ_std = 0;
    for (int d = 0; d < 3; d++) {
        displ_mean += displ_sum[d] / num_particles / 3.0;
        double var = displ_sum2[d] / num_particles - pow(displ_sum[d] / num_particles, 2);
        displ_std += sqrt(var) / 3.0;
    }
    
    stats[0] = displ_min[0];  // Use first component for min/max
    stats[1] = displ_max[0];
    stats[2] = displ_mean;
    stats[3] = displ_std;
    
    if (has_velocities) {
        double vel_mean = 0, vel_std = 0;
        for (int d = 0; d < 3; d++) {
            vel_mean += vel_sum[d] / num_particles / 3.0;
            double var = vel_sum2[d] / num_particles - pow(vel_sum[d] / num_particles, 2);
            vel_std += sqrt(var) / 3.0;
        }
        stats[4] = vel_min[0];
        stats[5] = vel_max[0];
        stats[6] = vel_mean;
        stats[7] = vel_std;
    }
    
    return 0;
}

void PrintOutputSummary(
    const char *output_dir,
    int rank,
    int N,
    OutputFormat format
) {
    char rank_dir[1024];
    snprintf(rank_dir, sizeof(rank_dir), "%s/rank_%d", output_dir, rank);
    
    printf("\n[Rank %d] Output Summary:\n", rank);
    printf("Directory: %s\n", rank_dir);
    
    // Count files
    int num_files = 0;
    long total_size = 0;
    
    for (int z = 0; z < N; z++) {
        char filename[1024];
        snprintf(filename, sizeof(filename), "%s/z%d_slab_N%d.bin", rank_dir, z, N);
        
        struct stat st;
        if (stat(filename, &st) == 0) {
            num_files++;
            total_size += st.st_size;
        }
    }
    
    printf("Files written: %d / %d Z-slabs\n", num_files, N);
    printf("Total size: %.2f MB\n", total_size / 1024.0 / 1024.0);
    printf("Average size per slab: %.2f KB\n", total_size / num_files / 1024.0);
    
    size_t particle_size = GetParticleSize(format);
    printf("Particle size: %zu bytes\n", particle_size);
}


