// ====================================================================================
// WRITE PARTICLES FROM REASSEMBLED MPI Z-SLABS
// ====================================================================================
// Reassembles per-rank i-slabs written by main.cpp and calls WriteParticlesSlab_new
// to write complete particle data (matching serial Zeldovich output).
//
// 1. Reads command-line arguments
// 2. Loads or creates sim parameters
// 3. Initializes output buffers
// 4. For each i-slab:
//    - Reassembles full i-slab from all ranks using ReassembleISlabFromRanks
//    - Extracts pointers to the 4 arrays
//    - Calls WriteParticlesSlab_new to write complete particle data
//
// Example:
//   ./write_particles_from_reassembled_mpi . 16 4 param_N16.par 0 16
// ====================================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <vector>
#include <string>
#include <filesystem>
#include <iostream>
#include <algorithm>

// Avoid MPI C++ binding conflicts
#define MPICH_SKIP_MPICXX

// Include zeldovich-PLT headers
#include <output.h>
#include <parameters.h>
#include <omp.h>

// Include output_new.h for WriteParticlesSlab_new
// Undefine MAX_PPD to avoid conflict
#ifdef MAX_PPD
#undef MAX_PPD
#endif
#include "output/output_new.h"

// Include grid decomposition utilities
extern "C" {
#include "utils/decomposition.h"
}

namespace fs = std::filesystem;

// Precision selection: match the precision used in main.cpp
// Use BinComplx to avoid conflict with zeldovich.h's Complx (which is always double)
#ifdef USE_DOUBLE_PRECISION
    using BinComplx = std::complex<double>;
    constexpr size_t EXPECTED_BYTES_PER_COMPLEX = 16;
#else
    using BinComplx = std::complex<float>;
    constexpr size_t EXPECTED_BYTES_PER_COMPLEX = 8;
#endif

// ====================================================================================
// REASSEMBLE I-SLAB FROM RANKS
// ====================================================================================
// Load per-rank i-slabs written by main.cpp:
//   rank_%d/i%d_slab_N%d.bin
// and assemble them into a full i-slab with all X and all Y:
//   layout: [narray][Y][X] stored as:
//       idx = array_idx * N * N + j * N + k
// This uses the grid decomposition to determine each rank's X-range.
// ====================================================================================

std::vector<BinComplx> ReassembleISlabFromRanks(
    const std::string &output_dir,
    int i,          // i index (Z slab)
    int N,
    int narray,
    int num_ranks
) {
    std::vector<BinComplx> full_slab;
    full_slab.resize((size_t)narray * (size_t)N * (size_t)N);

    // Zero-initialize
    std::fill(full_slab.begin(), full_slab.end(), BinComplx(0.0, 0.0));

    // Determine 2D grid factors for domain decomposition
    int grid_x = 0;
    int grid_z = 0;
    calculate_grid_factors(num_ranks, &grid_x, &grid_z);

    for (int rank = 0; rank < num_ranks; rank++) {
        // Build filename for this rank and slab
        std::string filename = output_dir + "/rank_" + std::to_string(rank) + 
                              "/i" + std::to_string(i) + "_slab_N" + std::to_string(N) + ".bin";

        if (!fs::exists(filename)) {
            continue;  // This rank does not own this i-slab
        }

        // Get X-range for this rank from grid decomposition
        GridBounds bounds = get_grid_bounds(rank, N, num_ranks);
        int x_start = bounds.x_start;
        int x_end   = bounds.x_end;
        int x_count_expected = x_end - x_start;

        // Open file and read contents
        FILE *fp = fopen(filename.c_str(), "rb");
        if (!fp) {
            fprintf(stderr,
                    "Error: cannot open %s for reading (rank=%d, errno=%d).\n",
                    filename.c_str(), rank, errno);
            continue;
        }

        // Determine file size
        if (fseek(fp, 0, SEEK_END) != 0) {
            fprintf(stderr, "Error: fseek failed for %s (rank=%d).\n",
                    filename.c_str(), rank);
            fclose(fp);
            continue;
        }
        long file_size = ftell(fp);
        if (file_size < 0) {
            fprintf(stderr, "Error: ftell failed for %s (rank=%d).\n",
                    filename.c_str(), rank);
            fclose(fp);
            continue;
        }
        rewind(fp);

        size_t elem_size = sizeof(BinComplx);
        if (file_size % (long)elem_size != 0) {
            fprintf(stderr,
                    "Warning: file %s size %ld not multiple of complex size %zu.\n",
                    filename.c_str(), file_size, elem_size);
        }

        // Infer x_count from file size: narray * N * x_count * sizeof(BinComplx)
        int64_t num_elems = file_size / (long)elem_size;
        int64_t denom = (int64_t)narray * (int64_t)N;
        if (denom == 0) {
            fprintf(stderr, "Error: invalid narray or N in ReassembleISlabFromRanks.\n");
            fclose(fp);
            continue;
        }
        int x_count = (int)(num_elems / denom);

        if (x_count != x_count_expected) {
            fprintf(stderr,
                    "Warning: %s x_count=%d differs from expected %d (rank=%d).\n",
                    filename.c_str(), x_count, x_count_expected, rank);
        }

        std::vector<BinComplx> local_slab;
        local_slab.resize((size_t)narray * (size_t)N * (size_t)x_count);

        size_t read_count = fread(local_slab.data(), elem_size,
                                  (size_t)narray * (size_t)N * (size_t)x_count, fp);
        fclose(fp);

        if (read_count != (size_t)narray * (size_t)N * (size_t)x_count) {
            fprintf(stderr,
                    "Warning: short read from %s (rank=%d, read %zu, expected %zu).\n",
                    filename.c_str(), rank, read_count,
                    (size_t)narray * (size_t)N * (size_t)x_count);
        }

        // Map from [array][Y][k_local] -> [array][Y][k_global]
        for (int array_idx = 0; array_idx < narray; array_idx++) {
            for (int j = 0; j < N; j++) {
                for (int k_local = 0; k_local < x_count; k_local++) {
                    int k_global = x_start + k_local;
                    if (k_global < 0 || k_global >= N) {
                        continue;
                    }

                    size_t src_idx =
                        (size_t)array_idx * (size_t)N * (size_t)x_count +
                        (size_t)j * (size_t)x_count +
                        (size_t)k_local;

                    size_t dst_idx =
                        (size_t)array_idx * (size_t)N * (size_t)N +
                        (size_t)j * (size_t)N +
                        (size_t)k_global;

                    full_slab[dst_idx] = local_slab[src_idx];
                }
            }
        }
    }

    return full_slab;
}

// ====================================================================================
// MAIN
// ====================================================================================

int main(int argc, char* argv[]) {
    // Force single OpenMP thread to avoid segfault in zeldovich-PLT initialization
    omp_set_num_threads(1);
    
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <output_dir> <N> <num_ranks> [param_file] [i_start] [i_end]\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "  output_dir: Directory containing rank_*/ subdirectories\n");
        fprintf(stderr, "  N: Resolution\n");
        fprintf(stderr, "  num_ranks: Number of MPI ranks\n");
        fprintf(stderr, "  param_file: Optional parameter file (default: create minimal)\n");
        fprintf(stderr, "  i_start: First i-slab to process (default: 0)\n");
        fprintf(stderr, "  i_end: Last i-slab to process (default: N)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Example:\n");
        fprintf(stderr, "  %s . 16 4 param_N16.par 0 16\n", argv[0]);
        return 1;
    }
    
    std::string output_dir = argv[1];
    int N = std::stoi(argv[2]);
    int num_ranks = std::stoi(argv[3]);
    
    std::string param_file = "";
    int i_start = 0;
    int i_end = N;
    
    if (argc >= 5) {
        param_file = argv[4];
    }
    if (argc >= 6) {
        i_start = std::stoi(argv[5]);
    }
    if (argc >= 7) {
        i_end = std::stoi(argv[6]);
    }
    
    printf("================================================================================\n");
    printf("WRITE PARTICLES FROM REASSEMBLED MPI Z-SLABS\n");
    printf("================================================================================\n");
    printf("Output directory: %s\n", output_dir.c_str());
    printf("Grid size N: %d\n", N);
    printf("Number of MPI ranks: %d\n", num_ranks);
    printf("i-slab range: [%d, %d)\n", i_start, i_end);
    printf("Parameter file: %s\n", param_file.empty() ? "(create minimal)" : param_file.c_str());
    printf("================================================================================\n\n");
    
    // Create Parameters object
    Parameters* param = nullptr;
    if (!param_file.empty() && fs::exists(param_file)) {
        try {
            param = new Parameters(fs::path(param_file));
            printf("Loaded parameters from: %s\n", param_file.c_str());
            printf("  ppd: %ld\n", param->ppd);
            printf("  boxsize: %f\n", param->boxsize);
            printf("  separation: %f\n", param->separation);
        } catch (const std::exception& e) {
            fprintf(stderr, "ERROR: Failed to load parameter file: %s\n", e.what());
            fprintf(stderr, "Creating minimal parameters...\n");
            param = nullptr;
        }
    }
    
    // Determine narray from parameter file (must match how main.cpp sets it)
    // This is critical for correctly reading .bin files!
    int narray;
    if (param) {
        if (param->qdensity == 2) {
            narray = 1;  // Density only
        } else {
            // Normal mode: Set narray based on PLT
            narray = param->qPLT ? 4 : 2;
        }
    } else {
        // Default to 2 if no param file (PLT disabled by default)
        narray = 2;
    }
    printf("Number of arrays: %d (determined from parameter file: qPLT=%d, qdensity=%d)\n", 
           narray, param ? param->qPLT : 0, param ? param->qdensity : 0);
    
    // Create minimal Parameters if needed
    if (!param) {
        std::string tmp_param_file = "/tmp/write_particles_reassembled_mpi_params.par";
        FILE* tmp_fp = fopen(tmp_param_file.c_str(), "w");
        if (tmp_fp) {
            fprintf(tmp_fp, "BoxSize = 1000.0\n");
            fprintf(tmp_fp, "ZD_Version = 2\n");
            fprintf(tmp_fp, "ZD_PPD = %d\n", N);
            fprintf(tmp_fp, "ZD_Pk_scale = 1.0\n");
            fprintf(tmp_fp, "NP = %lld\n", (long long)(int64_t)N * (int64_t)N * (int64_t)N);
            fprintf(tmp_fp, "ZD_NumBlock = 2\n");
            fprintf(tmp_fp, "CPD = %d\n", N);
            fprintf(tmp_fp, "ZD_Pk_norm = 1.0\n");
            fprintf(tmp_fp, "ZD_Pk_sigma = 1.0\n");
            fprintf(tmp_fp, "ZD_Pk_smooth = 0.0\n");
            fprintf(tmp_fp, "ZD_Pk_powerlaw_index = -2.0\n");
            fprintf(tmp_fp, "InitialConditionsDirectory = \"./output_complete\"\n");
            fprintf(tmp_fp, "InitialRedshift = 0.0\n");
            fprintf(tmp_fp, "ICFormat = \"RVZel\"\n");
            fprintf(tmp_fp, "ZD_Seed = 12345\n");
            fprintf(tmp_fp, "ZD_qPLT = 0\n");
            fprintf(tmp_fp, "ZD_qdensity = 1\n");
            fprintf(tmp_fp, "ZD_qascii = 0\n");
            fclose(tmp_fp);
            
            try {
                param = new Parameters(fs::path(tmp_param_file));
                printf("Created minimal parameters (ppd=%ld, boxsize=%f)\n", param->ppd, param->boxsize);
            } catch (const std::exception& e) {
                fprintf(stderr, "ERROR: Failed to create minimal parameters: %s\n", e.what());
                return 1;
            }
        } else {
            fprintf(stderr, "ERROR: Cannot create temporary parameter file\n");
            return 1;
        }
    }
    
    // Initialize output buffers (required by WriteParticlesSlab_new)
    printf("Initializing output...\n");
    SetupOutputDir(*param);
    double buffer_size = InitOutputBuffers(*param);
    printf("  Output directory: %s\n", param->output_dir.c_str());
    printf("  Buffer size: %.3f GiB\n", buffer_size);
    
    // WriteParticlesSlab_new writes to internal buffers, FILE* can be NULL
    FILE* output_fp = NULL;
    
    printf("\nProcessing i-slabs...\n");
    
    // Determine grid factors
    int grid_x, grid_z;
    calculate_grid_factors(num_ranks, &grid_x, &grid_z);
    printf("Grid decomposition: %d x %d = %d ranks\n", grid_x, grid_z, num_ranks);
    
    // Process each i-slab
    int processed = 0;
    for (int i = i_start; i < i_end; i++) {
        // Reassemble full i-slab from all ranks (as BinComplx - matches .bin file precision)
        std::vector<BinComplx> full_slab_bin = ReassembleISlabFromRanks(
            output_dir, i, N, narray, num_ranks
        );
        
        if (full_slab_bin.empty()) {
            printf("  i=%d: Skipping (reassembly failed or no data)\n", i);
            continue;
        }
        
        // Verify we have the right amount of data
        if (full_slab_bin.size() != (size_t)(narray * N * N)) {
            fprintf(stderr, "ERROR: i=%d wrong size: %zu != %d\n", i, full_slab_bin.size(), narray * N * N);
            continue;
        }
        
        // Convert from BinComplx (float or double, matches .bin file) to Complx (double, required by WriteParticlesSlab_new)
        // WriteParticlesSlab_new always expects Complx = std::complex<double>
        std::vector<Complx> full_slab;
        full_slab.reserve(full_slab_bin.size());
        for (const auto& val : full_slab_bin) {
            full_slab.push_back(Complx(static_cast<double>(val.real()), static_cast<double>(val.imag())));
        }
        
        // Create 2D slab pointers for WriteParticlesSlab_new
        // Layout: [narray][Y][X]
        // Mapping:
        //   displ[X] = imag(slab1) = imag(Array=0)
        //   displ[Y] = real(slab2) = real(Array=1)
        //   displ[Z] = imag(slab2) = imag(Array=1)
        Complx* slab1 = &full_slab[0 * N * N];  // Array 0: D + i*F (X displacement in imag part)
        Complx* slab2 = &full_slab[1 * N * N];  // Array 1: G + i*H (Y displacement in real part, Z displacement in imag part)
        Complx* slab3 = &full_slab[2 * N * N];  // Array 2: Z-velocity (if PLT enabled)
        Complx* slab4 = &full_slab[3 * N * N];  // Array 3: X-velocity + Y-velocity (if PLT enabled)
        
        // Call WriteParticlesSlab_new
        try {
            WriteParticlesSlab_new(output_fp, i, slab1, slab2, slab3, slab4, *param);
            processed++;
            printf("  i=%d: OK\n", i);
        } catch (const std::exception& e) {
            fprintf(stderr, "ERROR: WriteParticlesSlab_new failed for i=%d: %s\n", i, e.what());
        } catch (...) {
            fprintf(stderr, "ERROR: Unknown exception in WriteParticlesSlab_new for i=%d\n", i);
        }
    }
    
    // Note: output_fp is NULL since WriteParticlesSlab_new uses internal buffers
    // Don't call fclose(NULL) - it causes undefined behavior
    
    printf("\n================================================================================\n");
    printf("COMPLETE: Processed %d i-slabs\n", processed);
    printf("Output files written to: %s\n", param->output_dir.c_str());
    printf("  Files: ic_* (complete particle data, matches serial Zeldovich)\n");
    if (param->qdensity) {
        printf("  Density file: %s\n", param->density_filename.c_str());
    }
    printf("================================================================================\n");
    
    // Cleanup
    TeardownOutput();
    delete param;
    
    return 0;
}
