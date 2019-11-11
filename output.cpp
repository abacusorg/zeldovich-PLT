// #include "counts_in_cell.h"
// CountCell *cic;

double density_variance;
#define YX(_slab,_y,_x) _slab[(_x)+array.ppd*(_y)]

#define WRAP(_x) if (_x<0.0) _x += param.boxsize; \
    if (_x>=param.boxsize) _x -= param.boxsize;
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

void WriteParticlesSlab(FILE *output, FILE *densoutput, 
int z, Complx *slab1, Complx *slab2, Complx *slab3, Complx *slab4,
BlockArray& array, Parameters& param) {
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
    norm = 1.0; densitynorm = 1.0; vnorm = 1.0;
    //norm = 1e-2;

    char fn[1080];
    sprintf(fn, "%s/ic_%d",param.output_dir,z*param.cpd/param.ppd);
    //fprintf(stderr,"z: %d goes goes to ic_%d /n", z, z*param.cpd/param.ppd);
    output = fopen(fn,"ab");

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
            if(strcmp(param.ICFormat, "RVdoubleZel") == 0){
                RVdoubleZelParticle out;
                out.i = z; out.j =y; out.k = x;
                out.displ[0] = pos[2]; out.displ[1] = pos[1]; out.displ[2] = pos[0];
                out.vel[0] = vel[2]; out.vel[1] = vel[1]; out.vel[2] = vel[0];
                fwrite(&out,sizeof(out),1,output);
            } else if (strcmp(param.ICFormat, "RVZel") == 0){
                RVZelParticle out;
                out.i = z; out.j =y; out.k = x;
                out.displ[0] = pos[2]; out.displ[1] = pos[1]; out.displ[2] = pos[0];
                out.vel[0] = vel[2]; out.vel[1] = vel[1]; out.vel[2] = vel[0];
                fwrite(&out,sizeof(out),1,output);
            } else if (strcmp(param.ICFormat, "Zeldovich") == 0){
                ZelParticle out;
                out.i = z; out.j =y; out.k = x;
                out.displ[0] = pos[2]; out.displ[1] = pos[1]; out.displ[2] = pos[0];
                fwrite(&out,sizeof(ZelParticle),1,output);
            }
            else {
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
    fclose(output);
    return;
}