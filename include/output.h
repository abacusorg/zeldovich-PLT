#pragma once

#include <sys/stat.h>
#include <dirent.h>
#include <libgen.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>     /* PATH_MAX */
#include <errno.h>
#include <mutex>

#include "zeldovich.h"
#include "parameters.h"
#include "block_array.h"

extern double density_variance;
extern double max_disp[3];

class ZelParticle {
public:
    unsigned short i,j,k;
    double displ[3];
};

class ZelSimpleParticle {
public:
    float displ[3];
};

class RVZelParticle {
public:
    unsigned short i,j,k;
    float displ[3];
    float vel[3];
};

class RVdoubleZelParticle {
public:
    unsigned short i,j,k;
    double displ[3];
    double vel[3];
};

enum OutputType {OUTPUT_ZEL, OUTPUT_RVZEL, OUTPUT_RVDOUBLEZEL, OUTPUT_ZEL_SIMPLE};

void WriteParticlesSlab(FILE *output, 
    int z, Complx *slab1, Complx *slab2, Complx *slab3, Complx *slab4,
    BlockArray& array, Parameters& param);

void SetupOutputDir(Parameters &param);
int CleanDirectory(const char *path);
int CreateDirectories(const char *path);

// Returns GiB size of allocated buffer
double InitOutputBuffers(Parameters &param);
void TeardownOutput();
