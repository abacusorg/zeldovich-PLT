#ifndef PLT_EIGENMODES_H
#define PLT_EIGENMODES_H

#include "../types.h"  // For eigenmode structure
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Load PLT eigenmodes from binary file
// Returns: 0 on success, non-zero on error
// Side effects: Sets global eig_vecs and eig_vecs_ppd
// File format:
//   - First 4 bytes: eig_vecs_ppd (32-bit int)
//   - Then: eig_vecs_ppd * eig_vecs_ppd * (eig_vecs_ppd/2 + 1) * 4 doubles
//   - Each eigenmode is 4 doubles: vec[3] + val
int plt_load_eigenmodes(const char *filename);

// Free emode data
void plt_free_eigenmodes(void);

// Get the ppd (grid size) of the loaded eigenmodes
// Returns: eig_vecs_ppd, or 0 if not loaded
int64_t plt_get_eigenmode_ppd(void);

// Get PLT eigenmode for given k-vector using interpolation
int plt_get_eigenmode(int ikx, int iky, int ikz, int64_t ppd, eigenmode *e);

#ifdef __cplusplus
}
#endif

#endif
