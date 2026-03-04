#ifndef HERMITIAN_TYPES_H
#define HERMITIAN_TYPES_H

#include "precision.h"

// ====================================================================================
// Structure for Particle Linear Theory (PLT) eigenmodes
// Used when qPLT is enabled to get correct growth rates and eigenvectors
// ====================================================================================

typedef struct {
    double vec[3];  // evector components [x, y, z]
    double val;     // eval = growth rate
} eigenmode;

// ====================================================================================
// Grid bounds structure for 2D decomposition in (X,Z) 
// Ranges: [x_start, x_end) and [z_start, z_end)
// ====================================================================================

typedef struct {
    int x_start, x_end;  // X-range [x_start, x_end)
    int z_start, z_end;  // Z-range [z_start, z_end)
} GridBounds;

// ====================================================================================
// Extended grid bounds with padded (periodic, may overlap with neighbors) regions
// Core region  : Non-overlapping, primary responsibility of this rank
// Padded region: Extended with X_PADDING on each side, may overlap with neighbors
//
// Ex. (N=1024, 3x3 decomp, X_PADDING=10):
//   Rank 0 core: X=[0, 342), Z=[0, 342)
//   Rank 0 padded: X=[-10, 352), Z=[-10, 352)  // Wraps at boundaries
// ====================================================================================

typedef struct {
    GridBounds core;              // Core [x_start, x_end) 
    GridBounds padded;            // Padded region with overlap
    int num_pencils_core;         // # pencils in core
    int num_pencils_padded;       // # pencils in padded region
} ExtendedGridBounds;

// ====================================================================================
// Y-slice idxing: Access array 'array_idx' in slice 'slice_idx' at (x, z)
// Memory order: [Slice][Array][Z][X] (X is stride-1)
// Formula: x + N * (z + N * (array_idx + narray * slice_idx))
//
// Ex. N=8, narray=4, slice_idx=2, array_idx=1, x=3, z=5:
//   Index = 3 + 8*(5 + 8*(1 + 4*2)) = 3 + 8*(5 + 8*9) = 3 + 8*77 = 619
// ====================================================================================

/* Use (int64_t)(N) to avoid overflow: N*N can exceed INT_MAX for N > ~46340 */
#define Y_SLICE(slice_idx, array_idx, x, z, N, narray) \
    local_y_slices[(int64_t)(x) + (int64_t)(N) * ((z) + (int64_t)(N) * ((array_idx) + (narray) * (slice_idx)))]

// ====================================================================================
// Pencil indexing: Access array 'array_idx' in pencil 'pencil_idx' at Y position 'y'
// Memory order: [Pencil][Array][Y] (Y stride-1 cuz 1D FFT)
// Formula: y + N * (array_idx + narray * pencil_idx)
// ====================================================================================

#define PENCIL(pencil_idx, array_idx, y, N, narray) \
    local_pencils[(int64_t)(y) + (N) * ((array_idx) + (narray) * (pencil_idx))]

// ====================================================================================
// Z-slab idxing: [Array][X][Y]
// Formula: y + N * (x_idx + x_count * array_idx)
// Note: This format has better cache locality than [X][Array][Y] for unpacking
// For final output, transpose to [Array][Y][X] format during write!
// ====================================================================================
#define ZSLAB(array_idx, x_idx, y, N, narray, x_count) \
    local_z_slab[(int64_t)(y) + (N) * ((x_idx) + (x_count) * (array_idx))]

// ====================================================================================
// PBC MACROS
// ====================================================================================
// Maps logical coordinates (can be negative or >= N) to actual array indices [0, N).
//
// Examples (N=1024, X_PADDING=10):
//   PERIODIC_X(50, 1024)   -> 50    (common case, no-op)
//   PERIODIC_X(-5, 1024)   -> 1019  (wraps from right: X=-5 -> X=N-5)
//   PERIODIC_X(-10, 1024)  -> 1014  (wraps from right: X=-10 -> X=N-10)
//   PERIODIC_X(1030, 1024) -> 6     (wraps to left: X=N+6 -> X=6)

#if USE_X_PADDING
    #define PERIODIC_X(x, N) \
        ((x) >= 0 ? \
         ((x) < (N) ? (x) : (x) - (N)) : \
         (((x) % (N) + (N)) % (N)))
    
    #define PERIODIC_Z(z, N) \
        ((z) >= 0 ? \
         ((z) < (N) ? (z) : (z) - (N)) : \
         (((z) % (N) + (N)) % (N)))
#else
    // No padding: coordinates are always in [0, N), so no wrapping needed
    #define PERIODIC_X(x, N) (x)
    #define PERIODIC_Z(z, N) (z)
#endif

// ====================================================================================
// SIZE AND COUNT FUNCTIONS
// ====================================================================================
// Inline functions for computing sizes and counts based on grid configuration.

// Get number of elements in a GridBounds region
static inline int get_grid_bounds_size(const GridBounds *bounds) {
    int x_count = bounds->x_end - bounds->x_start;
    int z_count = bounds->z_end - bounds->z_start;
    return x_count * z_count;
}

// Get X-count from GridBounds
static inline int get_x_count(const GridBounds *bounds) {
    return bounds->x_end - bounds->x_start;
}

// Get Z-count from GridBounds
static inline int get_z_count(const GridBounds *bounds) {
    return bounds->z_end - bounds->z_start;
}

// Get total memory size for Y-slices
static inline size_t get_y_slice_memory_bytes(int num_slices, int N, int narray) {
    return (size_t)num_slices * narray * N * N * BYTES_PER_COMPLEX;
}

// Get total memory size for pencils
static inline size_t get_pencil_memory_bytes(int num_pencils, int N, int narray) {
    return (size_t)num_pencils * narray * N * BYTES_PER_COMPLEX;
}

// Get total memory size for Z-slab
static inline size_t get_z_slab_memory_bytes(int x_count, int N, int narray) {
    return (size_t)x_count * narray * N * BYTES_PER_COMPLEX;
}

// ====================================================================================
// Grid debugging prints
// ====================================================================================
static inline void print_grid_bounds(const GridBounds *bounds, const char *label, int rank) {
    printf("[Rank %d] %s: X=[%d, %d), Z=[%d, %d), size=%d pencils\n",
           rank, label,
           bounds->x_start, bounds->x_end,
           bounds->z_start, bounds->z_end,
           get_grid_bounds_size(bounds));
}

static inline void print_extended_grid_bounds(const ExtendedGridBounds *ext_bounds, int rank) {
    printf("[Rank %d] Extended Grid Bounds:\n", rank);
    print_grid_bounds(&ext_bounds->core, "  Core", rank);
    print_grid_bounds(&ext_bounds->padded, "  Padded", rank);
    printf("[Rank %d]   Core pencils: %d, Padded pencils: %d\n",
           rank, ext_bounds->num_pencils_core, ext_bounds->num_pencils_padded);
}

#endif

