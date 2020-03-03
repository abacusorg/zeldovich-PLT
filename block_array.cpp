#include "STimer.cc"

// We use a set of X-Z arrays of Complx numbers (ordered by A and Y).
// ppd is 64-bit, so this expression should be safe from overflow
#define BLK_AYZX(_slab,_a,_y,_z,_x) _slab[(int64_t)(_x)+ppd*((_z)+ppd*((_a)+narray*(_y)))]
// We use a set of X-Y arrays of Complx numbers (ordered by A and Z).
#define BLK_AZYX(_slab,_a,_z,_y,_x) _slab[(int64_t)(_x)+ppd*((_y)+ppd*((_a)+narray*(_z)))]

class BlockArray {
    // double complex data[zblock=0..NB-1][yblock=0..NB-1]
    //     [arr=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
    unsigned long long size;
    
    STimer wtimer;  // write timer
    STimer rtimer;  // read timer
    int64_t bytes_written = 0;
    int64_t bytes_read = 0;

#ifndef DISK
    Complx *arr;   // But this array may never be allocated!
#endif

public:
    int numblock, block;
    int64_t ppd;  // Making ppd a 64-bit integer avoids trouble with ppd^3 int32_t overflow
    int64_t narray;  // Make int64 to avoid overflow in later calculations
    char TMPDIR[1024];
    int ramdisk;
#ifdef DISK
    Complx *StoreBlock_tmp;  // One extra block for fast transposes
#endif

    BlockArray(int _ppd, int _numblock, int _narray, char *_dir, int _ramdisk) {
        ppd = (int64_t) _ppd;
        numblock = _numblock;
        block = ppd/numblock;
        narray = _narray;
        strcpy(TMPDIR,_dir);
        assert(ppd%2==0);    // PPD must be even, due to incomplete Nyquist code
        assert(numblock%2==0);    // Number of blocks must be even
        assert(ppd==numblock*block);   // We'd like the blocks to divide evenly
        size = ppd*ppd*ppd*narray;
        ramdisk = _ramdisk;  // just to silence the compiler; might be overriden below

#ifdef DISK
        StoreBlock_tmp = new Complx[(int64_t)narray*block*block*ppd];
#ifdef DIRECTIO
        fileoffset = 0;
        diskbuffer = 1024*512;  // Magic number pulled from io_dio.cpp
#endif
#else
        arr = new Complx[size];
#endif

        // We parallelize over planes within a block.
        // If there are fewer planes than threads; we're underutilizing the CPU!
        // But it probably only matters for large problem sizes.
        int nthread = omp_get_max_threads();
        fprintf(stderr, "Using %d OpenMP threads\n", nthread);
        if(block < nthread && ppd >= 512){
            fprintf(stderr, R"(
*** Note: the number of particles per block (%d) is fewer than the number of threads (%d),
    so the CPU will be under-utilized.  You may wish to decrease ZD_NumBlock (%d) if memory
    constraints allow.

)", block, nthread, numblock);
        }
    }
    ~BlockArray() {
#ifdef DISK
        delete[] StoreBlock_tmp;

        // Clean up the "zeldovich.*.*" files
        for(int yblock = 0; yblock < numblock; yblock++){
            for(int zblock = 0; zblock < numblock; zblock++){
                int n = snprintf(filename, 1024, "%s/zeldovich.%1d.%1d",TMPDIR,yblock,zblock);
                assert(n < 1024);
                remove(filename);
            }
        }

        fprintf(stderr, "Block IO took %.2f sec to write %.2f GB ==> %.1f MB/sec\n",
            wtimer.Elapsed(), bytes_written/1e9, bytes_written/1e6/wtimer.Elapsed());
        fprintf(stderr, "Block IO took %.2f sec to read %.2f GB ==> %.1f MB/sec\n",
            rtimer.Elapsed(), bytes_read/1e9, bytes_read/1e6/rtimer.Elapsed());

#else
        delete []arr;
#endif
    }

#ifdef DISK
#ifdef DIRECTIO
        // These routines are for DIRECTIO
private:
    char filename[1024];
    off_t fileoffset;
    int diskbuffer;

    void bopen(int yblock, int zblock, const char *mode) {
        // DirectIO actually opens and closes the files on demand, so we don't need to open the file here
        assert(yblock>=0&&yblock<numblock);
        assert(zblock>=0&&zblock<numblock);
        sprintf(filename,"%s/zeldovich.%1d.%1d",TMPDIR,yblock,zblock);
        FILE * outfile = fopen(filename,"a");
        assert(outfile != NULL);
        fclose(outfile);

        fileoffset = 0;  // We're opening a new file, so reset the file offset
        return;
    }
    void bclose() { return; }
    void bwrite(Complx *buffer, size_t num) {
        // Write num Complx numbers to the buffer, increment the pointer
        WriteDirect WD(ramdisk, diskbuffer);
        int sizebytes = num*sizeof(Complx);
        WD.BlockingAppend(filename, (char*)buffer, sizebytes);
    }
    void bread(Complx *buffer, size_t num) {
        // Read num Complx numbers into the buffer, increment the pointer
        ReadDirect RD(ramdisk, diskbuffer);
        size_t sizebytes = num*sizeof(Complx);
        RD.BlockingRead( filename, (char*)buffer, sizebytes, fileoffset);
        fileoffset += sizebytes;
    }
#else
    // These routines are for reading blocks on and off disk without DIRECTIO
private: 
    FILE *fp;
    char filename[1024];

public:
    void StoreBlock(int yblock, int zblock, Complx *slab) {
        // We must be sure to store the block sequentially.
        // data[zblock=0..NB-1][yblock=0..NB-1]
        //     [array=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
        wtimer.Start();

        int64_t i = 0;

        unsigned int a;
        int yres,zres,z;
        for (a=0;a<narray;a++) 
            for (zres=0;zres<block;zres++) {
                z = zres+block*zblock;
                for (yres=0;yres<block;yres++) {
                    // Copy the whole X skewer
                    for(int x = 0; x < ppd; x++)
                        StoreBlock_tmp[i++] = BLK_AYZX(slab,a,yres,z,x);
                }
            }

        assert(i == (int64_t)narray*block*block*ppd);

        bopen(yblock,zblock,"w");
        bwrite(StoreBlock_tmp, narray*block*block*ppd);
        bclose();

        wtimer.Stop();

        int64_t totsize = narray*block*block*ppd*sizeof(Complx);
        bytes_written += totsize;
        //fprintf(stderr, "StoreBlock took %.3g sec to write %.3g MB ==> %.3g MB/sec\n",
        //    wtimer.Elapsed(), totsize/1e6, totsize/1e6/wtimer.Elapsed());
    }

    void LoadBlock(int yblock, int zblock, Complx *slab) {
        rtimer.Start();

        // We must be sure to access the block sequentially.
        // data[zblock=0..NB-1][yblock=0..NB-1]
        //     [array=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
        // Can't openMP an I/O loop.
        unsigned int a;
        int yres,y,zres,yshift;
        
        // TODO: could probably just read in two halves
        bopen(yblock,zblock,"r");
        bread(StoreBlock_tmp, ppd*narray*block*block);
        bclose();

        int64_t i = 0;
        for (a=0;a<narray;a++)
            for (zres=0;zres<block;zres++) 
                for (yres=0;yres<block;yres++) {
                    y = yres+block*yblock;
                    // Copy the whole X skewer.  However, we want to
                    // shift the y frequencies in the reflected half
                    // by one.
                    // FLAW: Assumes ppd is even.
                    if (y>=ppd/2) yshift=y+1; else yshift=y;
                    if (yshift==ppd) yshift=ppd/2;
                    // Put it somewhere; this is about to be overwritten
                    for(int x = 0; x < ppd; x++)
                        BLK_AZYX(slab,a,zres,yshift,x) = StoreBlock_tmp[i++];
                }
        assert(i == ppd*narray*block*block);

        rtimer.Stop();
        int64_t totsize = narray*block*block*ppd*sizeof(Complx);
        bytes_read += totsize;
        //fprintf(stderr, "LoadBlock took %.3g sec to read %.3g MB ==> %.3g MB/sec\n",
        //        rtimer.Elapsed(), totsize/1e6, totsize/1e6/rtimer.Elapsed());

        return;
    }

private:

    void bopen(int yblock, int zblock, const char *mode) {
        // Set up for reading or writing this block
        assert(yblock>=0&&yblock<numblock);
        assert(zblock>=0&&zblock<numblock);
        sprintf(filename,"%s/zeldovich.%1d.%1d",TMPDIR,yblock,zblock);

        fp = fopen(filename,mode);
        if(fp == NULL) fprintf(stderr,"bad filename: %s\n",filename);
        assert(fp!=NULL);
        
        return;
    }
    void bclose() { fclose(fp); return; }
    void bwrite(Complx *buffer, size_t num) {
        // Write num Complx numbers to the buffer, increment the pointer
        fwrite(buffer,sizeof(Complx),num,fp);
    }
    void bread(Complx *buffer, size_t num) {
        // Read num Complx numbers into the buffer, increment the pointer
        size_t nread = fread(buffer,sizeof(Complx),num,fp);
        assert(nread == num);
    }
#endif
#else
    // These routines are for reading in and out of a big array in memory

private: 
    Complx *IOptr;
public:
    void StoreBlock(int yblock, int zblock, Complx *slab) {
        wtimer.Start();

        // We must be sure to store the block sequentially.
        // data[zblock=0..NB-1][yblock=0..NB-1]
        //     [array=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
        // Can't openMP an I/O loop.
        unsigned int a;
        int yres,zres,z;
        bopen(yblock,zblock,"w");
        for (a=0;a<narray;a++) 
        for (zres=0;zres<block;zres++) 
        for (yres=0;yres<block;yres++) {
            z = zres+block*zblock;
            //y = yres+block*yblock;  // the slab is in the y coordinate, so we never need the absolute y index
            // Copy the whole X skewer
            bwrite(&(BLK_AYZX(slab,a,yres,z,0)),ppd);
        }
        bclose();

        wtimer.Stop();
        bytes_written += narray*block*block*ppd*sizeof(Complx);
    }

    void LoadBlock(int yblock, int zblock, Complx *slab) {
        rtimer.Start();

        // We must be sure to access the block sequentially.
        // data[zblock=0..NB-1][yblock=0..NB-1]
        //     [array=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
        // Can't openMP an I/O loop.
        unsigned int a;
        int yres,y,zres,yshift;
        bopen(yblock,zblock,"r");
        for (a=0;a<narray;a++)
        for (zres=0;zres<block;zres++) 
        for (yres=0;yres<block;yres++) {
            //z = zres+block*zblock;  // slabs are in z; never need the absolute z coord
            y = yres+block*yblock;
            // Copy the whole X skewer.  However, we want to
            // shift the y frequencies in the reflected half
            // by one.
            // FLAW: Assumes ppd is even.
            if (y>=ppd/2) yshift=y+1; else yshift=y;
            if (yshift==ppd) yshift=ppd/2;
            // Put it somewhere; this is about to be overwritten
            bread(&(BLK_AZYX(slab,a,zres,yshift,0)),ppd);
        }
        bclose();

        rtimer.Stop();
        int64_t totsize = narray*block*block*ppd*sizeof(Complx);
        bytes_read += totsize;
        return;
    }

private: 

    void bopen(int yblock, int zblock, const char *mode) {
        // Set up for reading or writing this block
        assert(yblock>=0&&yblock<numblock);
        assert(zblock>=0&&zblock<numblock);
        IOptr = arr+((int64_t)zblock*numblock+yblock)*(block*block*ppd*narray);
        return;
    }
    void bclose() { IOptr = NULL; return; }
    void bwrite(Complx *buffer, size_t num) {
        // Write num Complx numbers to the buffer, increment the pointer
        memcpy(IOptr,buffer,sizeof(Complx)*num); IOptr+=num;
    }
    void bread(Complx *buffer, size_t num) {
        // Read num Complx numbers into the buffer, increment the pointer
        memcpy(buffer,IOptr,sizeof(Complx)*num); IOptr+=num;
    }
#endif
};