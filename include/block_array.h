#pragma once

#include <complex>
#include <mutex>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <filesystem>

#include "STimer.h"
#include "zeldovich.h"

namespace fs = std::filesystem;

#ifdef DIRECTIO
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "file.h"
#endif

// We use a set of X-Z arrays of Complx numbers (ordered by A and Y).
// ppd is 64-bit, so this expression should be safe from overflow
#define BLK_AYZX(_slab, _a, _y, _z, _x)                                 \
    _slab[(int64_t) (_x) + ppd * ((_z) + ppd * ((_a) + narray * (_y)))]
// We use a set of X-Y arrays of Complx numbers (ordered by A and Z).
#define BLK_AZYX(_slab, _a, _z, _y, _x)                                 \
    _slab[(int64_t) (_x) + ppd * ((_y) + ppd * ((_a) + narray * (_z)))]

class BlockArray {
    // double complex data[zblock=0..NB-1][yblock=0..NB-1]
    //     [arr=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
    unsigned long long size;
    std::mutex array_mutex;

#ifndef DISK
    Complx *arr;  // But this array may never be allocated!
#endif

public:
    // These are some stats to accumulate
    STimer wtimer;  // write timer
    STimer rtimer;  // read timer
    int64_t bytes_written = 0;
    int64_t bytes_read    = 0;
    int files_written     = 0;

    int numblock, block;
    int64_t
       ppd;  // Making ppd a 64-bit integer avoids trouble with ppd^3 int32_t overflow
    int64_t ppdhalf;  // ppd/2
    int64_t narray;   // Make int64 to avoid overflow in later calculations
    fs::path TMPDIR;
    int ramdisk;
    int quickdelete;
    int part;

    BlockArray(
       int _ppd,
       int _numblock,
       int _narray,
       const fs::path &_dir,
       int _ramdisk,
       int _quickdelete,
       int _part
    );
    ~BlockArray();

private:
#ifdef DISK
#ifdef DIRECTIO_BROKEN  // BROKEN CODE, I think.
    // These routines are for DIRECTIO
    off_t fileoffset;
    int diskbuffer;

    void bopen(int yblock, int zblock, const std::string &mode);
    void bclose();
    void bwrite(Complx *buffer, size_t num);
    void bread(Complx *buffer, size_t num);
#else
    // These routines are for reading blocks on and off disk without DIRECTIO
private:
    FILE *bopen(int yblock, int zblock, const std::string &mode);
    void bclose(FILE *fp);
    void bwrite(FILE *fp, Complx *buffer, size_t num);
    void bread(FILE *fp, Complx *buffer, size_t num);
#endif

    void bremove(int yblock, int zblock);

public:
    void StoreBlock(int yblock, int zblock, Complx *slab);
    void LoadBlock(int yblock, int zblock, Complx *slab);
    void StoreBlockForward(int yblock, int zblock, Complx *slab);
    void LoadBlockForward(int yblock, int zblock, Complx *slab);

#else  // not DISK
    // These routines are for reading in and out of a big array in memory

public:
    void StoreBlock(int yblock, int zblock, Complx *slab);
    void LoadBlock(int yblock, int zblock, Complx *slab);
    void StoreBlockForward(int yblock, int zblock, Complx *slab);
    void LoadBlockForward(int yblock, int zblock, Complx *slab);

private:
    Complx *bopen(int yblock, int zblock, const std::string &mode);
    void bclose(Complx *&IOptr);
    void bwrite(Complx *&IOptr, Complx *buffer, size_t num);
    void bread(Complx *&IOptr, Complx *buffer, size_t num);
#endif
};
