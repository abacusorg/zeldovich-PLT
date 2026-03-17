// ====================================================================================
// PLT EIGENMODE (load + interpolate)
// ====================================================================================

#include "plt_eigenmodes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

// Same as zeldovich-PLT's global variables
static double *eig_vecs = NULL;
static int64_t eig_vecs_ppd = 0;

// ====================================================================================
// EIGENMODE FILE LOADING
// ====================================================================================

int plt_load_eigenmodes(const char *filename) {
    if (!filename || filename[0] == '\0') {
        fprintf(stderr, "[Error] PLT eigenmode filename is empty\n");
        return -1;
    }
    
    FILE *eigf = fopen(filename, "rb");
    if (!eigf) {
        fprintf(stderr, "[Error] Could not open eigenmode file \"%s\"\n", filename);
        return -1;
    }
    
    // Get file size
    fseek(eigf, 0, SEEK_END);
    long size = ftell(eigf);
    fseek(eigf, 0, SEEK_SET);
    
    // Read ppd (32-bit int)
    int32_t eig_vecs_ppd_32;
    if (fread(&eig_vecs_ppd_32, sizeof(eig_vecs_ppd_32), 1, eigf) != 1) {
        fprintf(stderr, "[Error] Could not read ppd from eigenmode file \"%s\"\n", filename);
        fclose(eigf);
        return -1;
    }
    eig_vecs_ppd = (int64_t)eig_vecs_ppd_32;
    
    // Calculate expected size
    int64_t halfppd = eig_vecs_ppd / 2 + 1;
    size_t nelem = (size_t)eig_vecs_ppd * eig_vecs_ppd * halfppd * 4;
    size_t nbytes = nelem * sizeof(double);
    size_t expected_size = nbytes + sizeof(eig_vecs_ppd_32);
    
    if ((size_t)size != expected_size) {
        fprintf(stderr, "[Error] Eigenmode file \"%s\" of size %ld did not match expected size %zu "
                "from eig_vecs_ppd %lld\n", filename, size, expected_size, (long long)eig_vecs_ppd);
        fclose(eigf);
        return -1;
    }
    
    // Allocate memory (aligned to 4096 bytes like zeldovich-PLT)
    int ret = posix_memalign((void **)&eig_vecs, 4096, nbytes);
    if (ret != 0) {
        fprintf(stderr, "[Error] Could not allocate memory for eigenmodes (posix_memalign failed)\n");
        fclose(eigf);
        return -1;
    }
    
    // Read eigenmode data
    if (fread(eig_vecs, sizeof(double), nelem, eigf) != nelem) {
        fprintf(stderr, "[Error] Could not read eigenmode data from file \"%s\"\n", filename);
        free(eig_vecs);
        eig_vecs = NULL;
        eig_vecs_ppd = 0;
        fclose(eigf);
        return -1;
    }
    
    fclose(eigf);
    
    fprintf(stderr, "[PLT] Loaded eigenmodes from \"%s\" (ppd=%lld, %zu elements)\n",
            filename, (long long)eig_vecs_ppd, nelem);
    
    return 0;
}

void plt_free_eigenmodes(void) {
    if (eig_vecs) {
        free(eig_vecs);
        eig_vecs = NULL;
        eig_vecs_ppd = 0;
    }
}

int64_t plt_get_eigenmode_ppd(void) {
    return eig_vecs_ppd;
}

// ====================================================================================
// EIGENMODE INTERPOLATION
// ====================================================================================

// Accessing eigenmode data (see zeldovich-PLT's EIGMODE macro)
#define EIGMODE(_kx, _ky, _kz, _i) \
    (eig_vecs[(int64_t)(_kx) * eig_vecs_ppd * halfppd * 4 + (_ky) * halfppd * 4 + (_kz) * 4 + (_i)])

int plt_get_eigenmode(int ikx, int iky, int ikz, int64_t ppd, eigenmode *e) {
    if (!eig_vecs || eig_vecs_ppd == 0) {
        fprintf(stderr, "[Error] Eigenmodes not loaded\n");
        return -1;
    }
    
    if (!e) {
        fprintf(stderr, "[Error] NULL eigenmode output pointer\n");
        return -1;
    }
    
    int64_t halfppd = eig_vecs_ppd / 2 + 1;
    int64_t ppdhalf = eig_vecs_ppd / 2;
    
    // Exact match case (no interpolation needed)
    if (eig_vecs_ppd % ppd == 0) {
        int64_t scale = eig_vecs_ppd / ppd;
        int64_t ekx = ikx * scale;
        int64_t eky = iky * scale;
        int64_t ekz = ikz * scale;
        
        // Handle negative k indices (wrap to positive)
        if (ekx < 0) ekx += eig_vecs_ppd;
        if (eky < 0) eky += eig_vecs_ppd;
        if (ekz < 0) ekz += eig_vecs_ppd;
        
        // Handle z index (only positive half-space stored)
        if (ekz > ppdhalf) ekz = eig_vecs_ppd - ekz;
        
        e->vec[0] = EIGMODE(ekx, eky, ekz, 0);
        e->vec[1] = EIGMODE(ekx, eky, ekz, 1);
        e->vec[2] = EIGMODE(ekx, eky, ekz, 2);
        e->val = EIGMODE(ekx, eky, ekz, 3);
        return 0;
    }
    
    // Interpolation case (matches zeldovich-PLT's interp_eigmode)
    double fx = ((double)eig_vecs_ppd) / ppd * ikx;
    double fy = ((double)eig_vecs_ppd) / ppd * iky;
    double fz = ((double)eig_vecs_ppd) / ppd * ikz;
    
    // Handle negative k indices
    if (fx < 0) fx += eig_vecs_ppd;
    if (fy < 0) fy += eig_vecs_ppd;
    if (fz < 0) fz += eig_vecs_ppd;
    
    // Don't interpolate across Nyquist boundary
    if (fx > ppdhalf && fx < halfppd) fx = floor(fx + 1);
    if (fy > ppdhalf && fy < halfppd) fy = floor(fy + 1);
    if (fz > ppdhalf && fz < halfppd) fz = floor(fz + 1);
    
    // Get grid point indices
    int ikx_l = (int)fx;
    int ikx_h = ikx_l + 1;
    int iky_l = (int)fy;
    int iky_h = iky_l + 1;
    int ikz_l = (int)fz;
    int ikz_h = ikz_l + 1;
    
    // Handle wrap-around at boundaries
    if (ikx_h == eig_vecs_ppd) ikx_h = 0;
    if (iky_h == eig_vecs_ppd) iky_h = 0;
    if (ikz_h == eig_vecs_ppd) ikz_h = 0;
    
    // Fractional positions
    fx -= ikx_l;
    fy -= iky_l;
    fz -= ikz_l;
    
    // Trilinear interpolation coefficients
    double f[8];
    f[0] = (1 - fx) * (1 - fy) * (1 - fz);
    f[1] = (1 - fx) * (1 - fy) * (fz);
    f[2] = (1 - fx) * (fy) * (1 - fz);
    f[3] = (1 - fx) * (fy) * (fz);
    f[4] = (fx) * (1 - fy) * (1 - fz);
    f[5] = (fx) * (1 - fy) * (fz);
    f[6] = (fx) * (fy) * (1 - fz);
    f[7] = (fx) * (fy) * (fz);
    
    // Handle z index wrap-around for interpolation
    int ekz_l = ikz_l;
    int ekz_h = ikz_h;
    if (ekz_l > ppdhalf) ekz_l = eig_vecs_ppd - ekz_l;
    if (ekz_h > ppdhalf) ekz_h = eig_vecs_ppd - ekz_h;
    
    // Interpolate each component
    for (int i = 0; i < 4; i++) {
        double result = 0.0;
        result += f[0] * EIGMODE(ikx_l, iky_l, ekz_l, i);
        result += f[1] * EIGMODE(ikx_l, iky_l, ekz_h, i);
        result += f[2] * EIGMODE(ikx_l, iky_h, ekz_l, i);
        result += f[3] * EIGMODE(ikx_l, iky_h, ekz_h, i);
        result += f[4] * EIGMODE(ikx_h, iky_l, ekz_l, i);
        result += f[5] * EIGMODE(ikx_h, iky_l, ekz_h, i);
        result += f[6] * EIGMODE(ikx_h, iky_h, ekz_l, i);
        result += f[7] * EIGMODE(ikx_h, iky_h, ekz_h, i);
        
        if (i < 3) {
            e->vec[i] = result;
        } else {
            e->val = result;
        }
    }
    
    return 0;
}
