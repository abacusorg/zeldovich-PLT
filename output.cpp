// #include "counts_in_cell.h"
// CountCell *cic;

#include <sys/stat.h>
#include <dirent.h>
#include <libgen.h>
#include <unistd.h>

#include <string.h>
#include <limits.h>     /* PATH_MAX */
#include <errno.h>

#define YX(_slab,_y,_x) _slab[(_x)+array.ppd*(_y)]

#define WRAP(_x) if (_x<0.0) _x += param.boxsize; \
    if (_x>=param.boxsize) _x -= param.boxsize;

double density_variance;
void *output_tmp;  // output buffer
STimer outtimer;
size_t output_bytes_written = 0;

class ZelParticle {
public:
    unsigned short i,j,k;
    double displ[3];
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

enum OutputType {OUTPUT_ZEL, OUTPUT_RVZEL, OUTPUT_RVDOUBLEZEL};
OutputType param_icformat;
size_t sizeof_outputtype;

void WriteParticlesSlab(FILE *output, FILE *densoutput, 
int z, Complx *slab1, Complx *slab2, Complx *slab3, Complx *slab4,
BlockArray& array, Parameters& param) {
    outtimer.Start();

    // Write out one slab of particles
    int x,y;
    double pos[4], vel[3];
    double norm, densitynorm, vnorm;
    // We also need to fix the normalizations, which come from many places:
    // 1) We used ik/k^2 to apply the velocities.
    //    We did correctly divide by param.fundamental when doing this.
    // 2) The inverse FFT carries a pre-factor of 1/N
    // 3) Our power spectrum convention requires a prefactor of sqrt(Volume)
    //    We fixed these in the power spectrum.
    // No care has been paid to the normalization of the velocity;
    // this is simply the displacement field.  However, this does happen
    // to be v/H, which is correct for the redshift-space displacement in EdS
    norm = 1.0; densitynorm = 1.0;
    
    // We may also need to apply an f_growth factor to the velocities
    // from the addition of a smooth, non-clustering background component
    // (indicated by f_cluster != 1)
    // If we used PLT, then we took care of this while computing
    // the PLT growth rate effects.  If not, take care of it here.
    if(param.qPLT){
        vnorm = 1.0;
    } else {
        vnorm = (sqrt(1. + 24*param.f_cluster) - 1)*.25;
    }
    //norm = 1e-2;

    int64_t i = 0;
    for (y=0;y<array.ppd;y++)
    for (x=0;x<array.ppd;x++) {
        // The displacements are in YX(slab,y,x) and the
        // base positions are in z,y,x
        //       pos[0] = x*param.separation+imag(YX(slab1,y,x))*norm;
        //       pos[1] = y*param.separation+real(YX(slab2,y,x))*norm;
        //       pos[2] = z*param.separation+imag(YX(slab2,y,x))*norm;
        pos[0] = imag(YX(slab1,y,x))*norm;
        pos[1] = real(YX(slab2,y,x))*norm;
        pos[2] = imag(YX(slab2,y,x))*norm;
        pos[3] = real(YX(slab1,y,x))*densitynorm;
        if(param.qPLT){
            vel[0] = imag(YX(slab3,y,x))*vnorm;
            vel[1] = real(YX(slab4,y,x))*vnorm;
            vel[2] = imag(YX(slab4,y,x))*vnorm;
        } else {
            vel[0] = imag(YX(slab1,y,x))*vnorm;
            vel[1] = real(YX(slab2,y,x))*vnorm;
            vel[2] = imag(YX(slab2,y,x))*vnorm;
            //vel[0] = 0; vel[1] = 0;vel[2] = 0;
        }
        //            WRAP(pos[0]);
        //            WRAP(pos[1]);
        //            WRAP(pos[2]);
        if (param.qascii) {
            fprintf(output,"%d %d %d %f %f %f %f %f %f %f\n",x,y,z,pos[0],pos[1],pos[2],pos[3], vel[0], vel[1], vel[2]);
        } else {
            switch(param_icformat){
                case OUTPUT_RVDOUBLEZEL: {
                    RVdoubleZelParticle out;
                    out.i = z; out.j =y; out.k = x;
                    out.displ[0] = pos[2]; out.displ[1] = pos[1]; out.displ[2] = pos[0];
                    out.vel[0] = vel[2]; out.vel[1] = vel[1]; out.vel[2] = vel[0];
                    ((RVdoubleZelParticle *) output_tmp)[i++] = out;
                    break;
                }

                case OUTPUT_RVZEL: {
                    RVZelParticle out;
                    out.i = z; out.j =y; out.k = x;
                    out.displ[0] = pos[2]; out.displ[1] = pos[1]; out.displ[2] = pos[0];
                    out.vel[0] = vel[2]; out.vel[1] = vel[1]; out.vel[2] = vel[0];
                    ((RVZelParticle *) output_tmp)[i++] = out;
                    break;
                }

                case OUTPUT_ZEL: {
                    ZelParticle out;
                    out.i = z; out.j =y; out.k = x;
                    out.displ[0] = pos[2]; out.displ[1] = pos[1]; out.displ[2] = pos[0];
                    ((ZelParticle *) output_tmp)[i++] = out;
                    break;
                }

                default:
                    fprintf(stderr, "Error: unknown ICFormat \"%s\". Aborting.\n", param.ICFormat);
                    exit(1);
            }
            //           fwrite(vel,sizeof(float),3,output);
            //           int *id;
            //           id = new int;
            //           *id = x+(y*param.ppd)+(z*param.ppd*param.ppd);
            //           if (*id <0) std::cout <<"Bad id: "<<*id<<"\n";
            //           fwrite(id,sizeof(int),1,output);
            //           fwrite(id,sizeof(int),1,output);
            //           if (densoutput!=NULL) fwrite(pos+3,sizeof(float),1,densoutput);
            //           delete id;
        }
        density_variance += pos[3]*pos[3];
        
        // Track the global max displacement
        for(int i = 0; i < 3; i++){
            max_disp[i] = pos[i] > max_disp[i] ? pos[i] : max_disp[i];
        }

        // cic->add_cic(param.boxsize,pos);
    }

    assert(i == param.ppd*param.ppd);

    char fn[1080];
    sprintf(fn, "%s/ic_%ld",param.output_dir,z*param.cpd/param.ppd);
    //fprintf(stderr,"z: %d goes goes to ic_%d /n", z, z*param.cpd/param.ppd);
    output = fopen(fn,"ab");
    fwrite(output_tmp, sizeof_outputtype, param.ppd*param.ppd, output);
    fclose(output);


    outtimer.Stop();
    int64_t totsize = array.ppd*array.ppd*sizeof_outputtype;
    output_bytes_written += totsize;
    
    return;
}

int CleanDirectory(const char *path);
int CreateDirectories(const char *path);

// Returns GiB size of allocated buffer
double InitOutput(Parameters &param){
    CleanDirectory(param.output_dir);
    CreateDirectories(param.output_dir);

    if(strcmp(param.ICFormat, "RVdoubleZel") == 0){
        param_icformat = OUTPUT_RVDOUBLEZEL;
        output_tmp = new RVdoubleZelParticle[param.ppd*param.ppd];
        sizeof_outputtype = sizeof(RVdoubleZelParticle);
    } else if (strcmp(param.ICFormat, "RVZel") == 0){
        param_icformat = OUTPUT_RVZEL;
        output_tmp = new RVZelParticle[param.ppd*param.ppd];
        sizeof_outputtype = sizeof(RVZelParticle);
    } else if (strcmp(param.ICFormat, "Zeldovich") == 0){
        param_icformat = OUTPUT_ZEL;
        output_tmp = new ZelParticle[param.ppd*param.ppd];
        sizeof_outputtype = sizeof(ZelParticle);
    }
    else {
        fprintf(stderr, "Error: unknown ICFormat \"%s\". Aborting.\n", param.ICFormat);
        exit(1);
    }

    return sizeof_outputtype*param.ppd*param.ppd/1024./1024./1024.;
}

void TeardownOutput(){
    switch(param_icformat){
        case OUTPUT_RVDOUBLEZEL:
            delete[] (RVdoubleZelParticle *) output_tmp;
            break;

        case OUTPUT_RVZEL:
            delete[] (RVZelParticle *) output_tmp;
            break;

        case OUTPUT_ZEL:
            delete[] (ZelParticle *) output_tmp;
            break;
    }

    printf("WriteParticlesSlab took %.3g sec to write %.3g MB ==> %.3g MB/sec\n",
            outtimer.Elapsed(), output_bytes_written/1e6, output_bytes_written/1e6/outtimer.Elapsed());
}

// Remove files named ic_* and zeldovich.* from the current directory
int CleanDirectory(const char *path){
    DIR *d = opendir(path);
    size_t path_len = strlen(path);
    int r = -1;

    if (d){
        struct dirent *p;

        r = 0;

        while (!r && (p=readdir(d))){
            int r2 = -1;
            char *buf;
            size_t len;

            /* Only process our target file names */
            if (strncmp(p->d_name, "ic_", 3) && strncmp(p->d_name, "zeldovich.", 10)){
                continue;
            }

            len = path_len + strlen(p->d_name) + 2; 
            buf = (char *) malloc(len);
            assert(buf != NULL);

            struct stat statbuf;

            snprintf(buf, len, "%s/%s", path, p->d_name);

            if (!lstat(buf, &statbuf)){  // stat or lstat?
                if (!S_ISDIR(statbuf.st_mode)){
                    r2 = unlink(buf);
                }
            }

            free(buf);

            r = r2;
        }

        closedir(d);
    }

    return r;
}

// A recursive mkdir function.
// This is semantically similar to Python's os.makedirs()
// from https://gist.github.com/JonathonReinhart/8c0d90191c38af2dcadb102c4e202950
int CreateDirectories(const char *path){
    const size_t len = strlen(path);
    char _path[PATH_MAX];
    char *p; 

    errno = 0;

    /* Copy string so it is mutable */
    if (len > sizeof(_path)-1) {
        errno = ENAMETOOLONG;
        return -1; 
    }   
    strcpy(_path, path);

    /* Iterate the string */
    for (p = _path + 1; *p; p++) {
        if (*p == '/') {
            /* Temporarily truncate */
            *p = '\0';

            if (mkdir(_path, S_IRWXU) != 0) {
                if (errno != EEXIST)
                    return -1; 
            }

            *p = '/';
        }
    }   

    if (mkdir(_path, S_IRWXU) != 0) {
        if (errno != EEXIST)
            return -1;
    }   

    return 0;
}
