#include <mutex>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>
#include <omp.h>

#include "STimer.h"
#include "block_array.h"

#ifdef DIRECTIO
#include "file.cpp"
#include "iolib.cpp"
#endif

BlockArray::BlockArray(int _ppd, int _numblock, int _narray, char *_dir, int _ramdisk, int _quickdelete, int _part) {
    ppd = (int64_t) _ppd;
    ppdhalf = ppd/2;
    numblock = _numblock;
    block = ppd/numblock;
    narray = _narray;
    strcpy(TMPDIR,_dir);
    assert(ppd%2==0);    // PPD must be even, due to incomplete Nyquist code
    assert(numblock%2==0);    // Number of blocks must be even
    assert(ppd==numblock*block);   // We'd like the blocks to divide evenly
    size = ppd*ppd*ppd*narray;
    ramdisk = _ramdisk;  // just to silence the compiler; might be overriden below
    quickdelete = _quickdelete;   // If set, then delete block_array files immediately upon reading, to save the space for the outputs.
    part = _part;  // part 1 (Z), part 2 (XY), or part -1 (both)

#ifdef DISK
#ifdef DIRECTIO
    fileoffset = 0;
    diskbuffer = 1024*512;  // Magic number pulled from io_dio.cpp
#endif

    if(part != 2){
        // Make directories for the zeldovich blocks
        int ret = mkdir(TMPDIR, 0775);
        int reason = errno;
        if(ret != 0 && reason != EEXIST){
            fprintf(stderr, "mkdir(\"%s\") failed for reason %s\n", TMPDIR, strerror(reason));
            exit(1);
        }
        for(int yblock = 0; yblock < numblock; yblock++){
            char blockdir[1100];
            sprintf(blockdir, "%s/zeldovich.%1d", TMPDIR, yblock);
            int ret = mkdir(blockdir, 0775);
            int reason = errno;
            if(ret != 0 && reason != EEXIST){
                fprintf(stderr, "mkdir(\"%s\") failed for reason %s\n", blockdir, strerror(reason));
                exit(1);
            }
        }
    }
#else
    //arr = new Complx[size];  // careful about hidden single-threaded zero-init!
    
    int ret = posix_memalign((void **) &arr, 4096, sizeof(Complx)*size);
    assert(ret == 0);
    #pragma omp parallel for schedule(static)
    for(uint64_t i = 0; i < size; i++){
        arr[i] = Complx(0.,0.);
    }
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

BlockArray::~BlockArray() {
#ifdef DISK
    if(part != 1){
        // Clean up the "zeldovich.*.*" files
        for(int yblock = 0; yblock < numblock; yblock++){
            if(!quickdelete){
                for(int zblock = 0; zblock < numblock; zblock++){
                    bremove(yblock, zblock);
                }
            }
            // Remove the block directory
            char blockdir[1100];
            sprintf(blockdir, "%s/zeldovich.%1d", TMPDIR, yblock);
            int ret = rmdir(blockdir);
            int reason = errno;
            if(ret != 0){
                fprintf(stderr, "rmdir(\"%s\") failed for reason %s\n", blockdir, strerror(reason));
                exit(1);
            }
        }
        // Remove the TMPDIR if it's empty
        rmdir(TMPDIR);
    }

    double serial_time = wtimer.Elapsed()/omp_get_max_threads();
    fprintf(stderr, "Block IO took %.2f sec to write %.2f GB ==> %.1f MB/sec\n",
        serial_time, bytes_written/1e9, bytes_written/1e6/serial_time);
    serial_time = rtimer.Elapsed()/omp_get_max_threads();
    fprintf(stderr, "Block IO took %.2f sec to read %.2f GB ==> %.1f MB/sec\n",
        serial_time, bytes_read/1e9, bytes_read/1e6/serial_time);
#else
    free(arr);
#endif
}

#ifdef DISK
#ifdef DIRECTIO_BROKEN
void BlockArray::bopen(int yblock, int zblock, const char *mode) {
    // DirectIO actually opens and closes the files on demand, so we don't need to open the file here
    assert(yblock>=0&&yblock<numblock);
    assert(zblock>=0&&zblock<numblock);
    char filename[1024];
    sprintf(filename,"%s/zeldovich.%1d/zeldovich.%1d.%1d",TMPDIR,yblock,yblock,zblock);
    FILE * outfile = fopen(filename,"a");
    assert(outfile != NULL);
    fclose(outfile);

    fileoffset = 0;  // We're opening a new file, so reset the file offset
    return;
}

void BlockArray::bclose() {
    return;
}

void BlockArray::bwrite(Complx *buffer, size_t num) {
    // Write num Complx numbers to the buffer, increment the pointer
    WriteDirect WD(ramdisk, diskbuffer);
    int sizebytes = num*sizeof(Complx);
    WD.BlockingAppend(filename, (char*)buffer, sizebytes);
}

void BlockArray::bread(Complx *buffer, size_t num) {
    // Read num Complx numbers into the buffer, increment the pointer
    ReadDirect RD(ramdisk, diskbuffer);
    size_t sizebytes = num*sizeof(Complx);
    RD.BlockingRead( filename, (char*)buffer, sizebytes, fileoffset);
    fileoffset += sizebytes;
}
  
#else

// These routines are for reading blocks on and off disk without DIRECTIO

FILE *BlockArray::bopen(int yblock, int zblock, const char *mode) {
    // Set up for reading or writing this block
    assert(yblock>=0&&yblock<numblock);
    assert(zblock>=0&&zblock<numblock);
    char filename[1100];
    sprintf(filename,"%s/zeldovich.%1d/zeldovich.%1d.%1d",TMPDIR,yblock,yblock,zblock);

    FILE *fp;
    fp = fopen(filename,mode);
    if(fp == NULL) fprintf(stderr,"bad filename: %s\n",filename);
    assert(fp!=NULL);
    
    return fp;
}
void BlockArray::bclose(FILE *fp) { fclose(fp); return; }
void BlockArray::bwrite(FILE *fp, Complx *buffer, size_t num) {
    // Write num Complx numbers to the buffer, increment the pointer
    fwrite(buffer,sizeof(Complx),num,fp);
}
void BlockArray::bread(FILE *fp, Complx *buffer, size_t num) {
    // Read num Complx numbers into the buffer, increment the pointer
    size_t nread = fread(buffer,sizeof(Complx),num,fp);
    assert(nread == num);
}
#endif

void BlockArray::bremove(int yblock, int zblock) {
    char filename[1100];
    sprintf(filename, "%s/zeldovich.%1d/zeldovich.%1d.%1d",TMPDIR,yblock,yblock,zblock);
    remove(filename);
}

void BlockArray::StoreBlock(int yblock, int zblock, Complx *slab) {
    // We must be sure to store the block sequentially.
    // data[zblock=0..NB-1][yblock=0..NB-1]
    //     [array=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
    STimer thiswtimer;  // write timer
    thiswtimer.Start();
    int64_t totsize = narray*block*block*ppd;
    Complx *StoreBlock_tmp;  // One extra block for fast transposes
    int ret = posix_memalign((void **)&StoreBlock_tmp, 4096, totsize*sizeof(Complx));
    assert(ret==0);

    int64_t i = 0;

    unsigned int a;
    int yres,zres,z;
    for (a=0;a<narray;a++) 
        for (zres=0;zres<block;zres++) {
            z = zres+block*zblock;
            for (yres=0;yres<block;yres++) {
                // Copy the whole X skewer
                memcpy(StoreBlock_tmp+i, &(BLK_AYZX(slab,a,yres,z,0)), ppd*sizeof(Complx));
                i+=ppd;
                // for(int x = 0; x < ppd; x++)
                    // StoreBlock_tmp[i++] = BLK_AYZX(slab,a,yres,z,x);
            }
        }

    assert(i == totsize);

    FILE *fp = bopen(yblock,zblock,"w");
    bwrite(fp, StoreBlock_tmp, totsize);
    bclose(fp);
    free(StoreBlock_tmp);

    thiswtimer.Stop();
    array_mutex.lock();
    wtimer.increment(thiswtimer.timer);
    bytes_written += totsize*sizeof(Complx);
    files_written++;
    //fprintf(stderr, "StoreBlock took %.3g sec to write %.3g MB ==> %.3g MB/sec\n",
    //    wtimer.Elapsed(), totsize/1e6, totsize/1e6/wtimer.Elapsed());
    array_mutex.unlock();
}

void BlockArray::LoadBlock(int yblock, int zblock, Complx *slab) {
    STimer thisrtimer;  // read timer
    thisrtimer.Start();

    // We must be sure to access the block sequentially.
    // data[zblock=0..NB-1][yblock=0..NB-1]
    //     [array=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
    // Can't openMP an I/O loop.
    unsigned int a;
    int yres,y,zres,yshift;
    
    // TODO: could probably just read in two halves
    int64_t totsize = narray*block*block*ppd;
    Complx *StoreBlock_tmp;  // One extra block for fast transposes
    int ret = posix_memalign((void **)&StoreBlock_tmp, 4096, totsize*sizeof(Complx));
    assert(ret==0);
    FILE *fp = bopen(yblock,zblock,"r");
    bread(fp,StoreBlock_tmp, totsize);
    bclose(fp);
    if (quickdelete) bremove(yblock, zblock);

    int64_t i = 0;
    for (a=0;a<narray;a++)
        for (zres=0;zres<block;zres++) 
            for (yres=0;yres<block;yres++) {
                y = yres+block*yblock;
                // Copy the whole X skewer.  However, we want to
                // shift the y frequencies in the reflected half
                // by one.
                // FLAW: Assumes ppd is even.
                if (y>=ppdhalf) yshift=y+1; else yshift=y;
                if (yshift==ppd) yshift=ppdhalf;
                // Put it somewhere; this is about to be overwritten
                memcpy(&(BLK_AZYX(slab,a,zres,yshift,0)), StoreBlock_tmp+i, ppd*sizeof(Complx));
                i+=ppd;
                // for(int x = 0; x < ppd; x++)
                    // BLK_AZYX(slab,a,zres,yshift,x) = StoreBlock_tmp[i++];
            }
    assert(i == ppd*narray*block*block);
    free(StoreBlock_tmp);

    thisrtimer.Stop();
    array_mutex.lock();
    rtimer.increment(thisrtimer.timer);
    bytes_read += totsize*sizeof(Complx);
    //fprintf(stderr, "LoadBlock took %.3g sec to read %.3g MB ==> %.3g MB/sec\n",
    //        rtimer.Elapsed(), totsize/1e6, totsize/1e6/rtimer.Elapsed());
    array_mutex.unlock();
    return;
}

void BlockArray::StoreBlockForward(int yblock, int zblock, Complx *slab) {
    STimer thiswtimer;  // write timer
    thiswtimer.Start();
    int64_t totsize = narray*block*block*ppd;
    Complx *StoreBlock_tmp;  // One extra block for fast transposes
    int ret = posix_memalign((void **)&StoreBlock_tmp, 4096, totsize*sizeof(Complx));
    assert(ret==0);

    int64_t i = 0;

    for (unsigned int a=0;a<narray;a++) 
        for (int yres=0;yres<block;yres++) {
            int y = yres+block*yblock;
            for (int zres=0;zres<block;zres++) {
                // Copy the whole X skewer
                memcpy(StoreBlock_tmp+i, &(BLK_AZYX(slab,a,zres,y,0)), ppd*sizeof(Complx));
                i+=ppd;
            }
        }

    assert(i == totsize);

    FILE *fp = bopen(yblock,zblock,"w");
    bwrite(fp, StoreBlock_tmp, totsize);
    bclose(fp);
    free(StoreBlock_tmp);

    thiswtimer.Stop();
    array_mutex.lock();
    wtimer.increment(thiswtimer.timer);
    bytes_written += totsize*sizeof(Complx);
    files_written++;
    //fprintf(stderr, "StoreBlock took %.3g sec to write %.3g MB ==> %.3g MB/sec\n",
    //    wtimer.Elapsed(), totsize/1e6, totsize/1e6/wtimer.Elapsed());
    array_mutex.unlock();
}

void BlockArray::LoadBlockForward(int yblock, int zblock, Complx *slab) {
    STimer thisrtimer;  // read timer
    thisrtimer.Start();
    
    int64_t totsize = narray*block*block*ppd;
    Complx *StoreBlock_tmp;  // One extra block for fast transposes
    int ret = posix_memalign((void **)&StoreBlock_tmp, 4096, totsize*sizeof(Complx));
    assert(ret==0);
    FILE *fp = bopen(yblock,zblock,"r");
    bread(fp,StoreBlock_tmp, totsize);
    bclose(fp);
    if (quickdelete) bremove(yblock, zblock);

    int64_t i = 0;
    for (unsigned int a=0;a<narray;a++)
        for (int yres=0;yres<block;yres++)
            for (int zres=0;zres<block;zres++) {
                int z = zres+block*zblock;
                memcpy(&(BLK_AYZX(slab,a,yres,z,0)), StoreBlock_tmp+i, ppd*sizeof(Complx));
                i+=ppd;
            }
    assert(i == ppd*narray*block*block);
    free(StoreBlock_tmp);

    thisrtimer.Stop();
    array_mutex.lock();
    rtimer.increment(thisrtimer.timer);
    bytes_read += totsize*sizeof(Complx);
    //fprintf(stderr, "LoadBlock took %.3g sec to read %.3g MB ==> %.3g MB/sec\n",
    //        rtimer.Elapsed(), totsize/1e6, totsize/1e6/rtimer.Elapsed());
    array_mutex.unlock();
    return;
}

#else   // not DISK
// These routines are for reading in and out of a big array in memory

void BlockArray::StoreBlock(int yblock, int zblock, Complx *slab) {
    STimer thiswtimer;
    thiswtimer.Start();

    // We must be sure to store the block sequentially.
    // data[zblock=0..NB-1][yblock=0..NB-1]
    //     [array=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
    // Can't openMP an I/O loop.
    unsigned int a;
    int yres,zres,z;
    Complx *IOptr = bopen(yblock,zblock,"w");
    for (a=0;a<narray;a++) 
    for (zres=0;zres<block;zres++) 
    for (yres=0;yres<block;yres++) {
        z = zres+block*zblock;
        //y = yres+block*yblock;  // the slab is in the y coordinate, so we never need the absolute y index
        // Copy the whole X skewer
        bwrite(IOptr, &(BLK_AYZX(slab,a,yres,z,0)),ppd);
    }
    bclose(IOptr);

    thiswtimer.Stop();
    array_mutex.lock();
    wtimer.increment(thiswtimer.timer);
    bytes_written += narray*block*block*ppd*sizeof(Complx);
    array_mutex.unlock();
}

void BlockArray::StoreBlockForward(int yblock, int zblock, Complx *slab) {
    STimer thiswtimer;
    thiswtimer.Start();

    unsigned int a;
    int yres,zres,y;
    Complx *IOptr = bopen(yblock,zblock,"w");
    for (a=0;a<narray;a++) 
    for (yres=0;yres<block;yres++)
    for (zres=0;zres<block;zres++)  {
        y = yres+block*yblock;
        // Copy the whole X skewer
        bwrite(IOptr, &(BLK_AZYX(slab,a,zres,y,0)),ppd);
    }
    bclose(IOptr);

    thiswtimer.Stop();
    array_mutex.lock();
    wtimer.increment(thiswtimer.timer);
    bytes_written += narray*block*block*ppd*sizeof(Complx);
    array_mutex.unlock();
}

void BlockArray::LoadBlockForward(int yblock, int zblock, Complx *slab) {
    STimer thisrtimer;
    thisrtimer.Start();

    unsigned int a;
    int yres,zres,z;
    Complx *IOptr = bopen(yblock,zblock,"r");
    for (a=0;a<narray;a++)
    for (yres=0;yres<block;yres++)
    for (zres=0;zres<block;zres++) {
        z = zres+block*zblock;
        // TODO: need shift for forward?
        // if (y>=ppdhalf) yshift=y+1; else yshift=y;
        // if (yshift==ppd) yshift=ppdhalf;
        bread(IOptr, &(BLK_AYZX(slab,a,yres,z,0)),ppd);
    }
    bclose(IOptr);

    thisrtimer.Stop();
    int64_t totsize = narray*block*block*ppd*sizeof(Complx);
    array_mutex.lock();
    rtimer.increment(thisrtimer.timer);
    bytes_read += totsize;
    array_mutex.unlock();
    return;
}

void BlockArray::LoadBlock(int yblock, int zblock, Complx *slab) {
    STimer thisrtimer;
    thisrtimer.Start();

    // We must be sure to access the block sequentially.
    // data[zblock=0..NB-1][yblock=0..NB-1]
    //     [array=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
    // Can't openMP an I/O loop.
    unsigned int a;
    int yres,y,zres,yshift;
    Complx *IOptr = bopen(yblock,zblock,"r");
    for (a=0;a<narray;a++)
    for (zres=0;zres<block;zres++) 
    for (yres=0;yres<block;yres++) {
        //z = zres+block*zblock;  // slabs are in z; never need the absolute z coord
        y = yres+block*yblock;
        // Copy the whole X skewer.  However, we want to
        // shift the y frequencies in the reflected half
        // by one.
        // FLAW: Assumes ppd is even.
        if (y>=ppdhalf) yshift=y+1; else yshift=y;
        if (yshift==ppd) yshift=ppdhalf;
        // Put it somewhere; this is about to be overwritten
        bread(IOptr, &(BLK_AZYX(slab,a,zres,yshift,0)),ppd);
    }
    bclose(IOptr);

    thisrtimer.Stop();
    int64_t totsize = narray*block*block*ppd*sizeof(Complx);
    array_mutex.lock();
    rtimer.increment(thisrtimer.timer);
    bytes_read += totsize;
    array_mutex.unlock();
    return;
}

Complx *BlockArray::bopen(int yblock, int zblock, const char *mode) {
    // Set up for reading or writing this block
    assert(yblock>=0&&yblock<numblock);
    assert(zblock>=0&&zblock<numblock);
    Complx *IOptr = arr+((int64_t)zblock*numblock+yblock)*(block*block*ppd*narray);
    return IOptr;
}

void BlockArray::bclose(Complx * &IOptr) { IOptr = NULL; return; }

void BlockArray::bwrite(Complx * &IOptr, Complx *buffer, size_t num) {
    // Write num Complx numbers to the buffer, increment the pointer
    memcpy(IOptr,buffer,sizeof(Complx)*num); IOptr+=num;
}

void BlockArray::bread(Complx * &IOptr, Complx *buffer, size_t num) {
    // Read num Complx numbers into the buffer, increment the pointer
    memcpy(buffer,IOptr,sizeof(Complx)*num); IOptr+=num;
}
#endif
