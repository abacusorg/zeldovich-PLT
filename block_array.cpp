class BlockArray {
    // double complex data[zblock=0..NB-1][yblock=0..NB-1]
    //     [arr=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
    unsigned long long size;
    Complx *arr;   // But this array may never be allocated!

public:
    int ppd, numblock, block;
    long long unsigned int narray;  // Make uint64 to avoid overflow in later calculations
    char TMPDIR[1024];
    int ramdisk;
    BlockArray(int _ppd, int _numblock, int _narray, char *_dir, int ramdisk) {
        ppd = _ppd;
        numblock = _numblock;
        block = ppd/numblock;
        narray = _narray;
        strcpy(TMPDIR,_dir);
        arr = NULL;
        assert(ppd%2==0);    // PPD must be even, due to incomplete Nyquist code
        assert(numblock%2==0);    // Number of blocks must be even
        assert(ppd==numblock*block);   // We'd like the blocks to divide evenly
        size = 1llu*ppd*ppd*ppd*narray;
#ifndef DISK
        arr = new Complx[size];
#elif defined DIRECTIO
        fileoffset = 0;
        diskbuffer = 1024*512;  // Magic number pulled from io_dio.cpp
        ramdisk = 0;
#endif
    }
    ~BlockArray() {
#ifdef DISK
        // Clean up the "zeldovich.*.*" files
        for(int yblock = 0; yblock < numblock; yblock++){
            for(int zblock = 0; zblock < numblock; zblock++){
                sprintf(filename,"%s/zeldovich.%1d.%1d",TMPDIR,yblock,zblock);
                remove(filename);
            }
        }
#else
        delete []arr;
#endif
    }

    //    Complx *point(int a, int z, int y, int x) {
    //    // Return a pointer to this part of the buffer
    //        int zblock = z/block;
    //      int yblock = y/block;
    //      return arr+x+ppd*(y-block*yblock+block*(z-block*zblock+block*(a+narray*(yblock+numblock*zblock))));
    //    }

#ifdef DISK
#ifdef DIRECTIO
        // These routines are for DIRECTIO
private:
    char filename[1024];
    off_t fileoffset;
    int diskbuffer;
public:
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
    // These routines are for reading blocks on and off disk
private: 
    FILE *fp;
    char filename[1024];
public:
    void bopen(int yblock, int zblock, const char *mode) {
        // Set up for reading or writing this block
        assert(yblock>=0&&yblock<numblock);
        assert(zblock>=0&&zblock<numblock);
        sprintf(filename,"%s/zeldovich.%1d.%1d",TMPDIR,yblock,zblock);

        fp = fopen(filename,mode);
        if(fp == NULL) printf("bad filename: %s\n",filename);
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
    void bopen(int yblock, int zblock, const char *mode) {
        // Set up for reading or writing this block
        assert(yblock>=0&&yblock<numblock);
        assert(zblock>=0&&zblock<numblock);
        IOptr = arr+(zblock*numblock+yblock)*(block*block*ppd*narray);
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