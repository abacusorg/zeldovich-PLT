#ifndef HERMITIAN_PRECISION_H
#define HERMITIAN_PRECISION_H

// ====================================================================================
// Default: Single precision (float)
// Double:  Compile with -DUSE_DOUBLE_PRECISION
// ====================================================================================

#include <fftw3.h>
#include <mpi.h>
#include <math.h>

#ifdef USE_DOUBLE_PRECISION

    typedef double real_t;
    typedef fftw_complex fftw_complex_t;
    typedef fftw_plan fftw_plan_t;
    
    // FFTW 
    #define FFTW_INIT_THREADS        fftw_init_threads
    #define FFTW_PLAN_WITH_NTHREADS  fftw_plan_with_nthreads
    #define FFTW_PLAN_DFT_2D        fftw_plan_dft_2d
    #define FFTW_PLAN_MANY_DFT      fftw_plan_many_dft
    #define FFTW_PLAN_DFT_1D        fftw_plan_dft_1d
    #define FFTW_EXECUTE_DFT        fftw_execute_dft
    #define FFTW_DESTROY_PLAN fftw_destroy_plan
    #define FFTW_MALLOC fftw_malloc
    #define FFTW_FREE fftw_free
    #define FFTW_IMPORT_WISDOM_FROM_FILENAME fftw_import_wisdom_from_filename
    #define FFTW_EXPORT_WISDOM_TO_FILENAME   fftw_export_wisdom_to_filename
    #define FFTW_EXPORT_WISDOM_TO_STRING     fftw_export_wisdom_to_string
    #define FFTW_IMPORT_WISDOM_FROM_STRING   fftw_import_wisdom_from_string
    #define FFTW_FORGET_WISDOM               fftw_forget_wisdom
    
    // MPI
    #define MPI_COMPLEX_TYPE MPI_C_DOUBLE_COMPLEX
    #define MPI_REAL_TYPE MPI_DOUBLE
    
    #define fabs_t fabs
    #define fmax_t fmax
    #define sqrt_t sqrt
    
    // Constants
    #define PRECISION_NAME "Double"
    #define BYTES_PER_COMPLEX 16
    
#else // Single precision (default)
    typedef float real_t;
    typedef fftwf_complex fftw_complex_t;
    typedef fftwf_plan fftw_plan_t;
    
    // FFTW 
    #define FFTW_INIT_THREADS        fftwf_init_threads
    #define FFTW_PLAN_WITH_NTHREADS  fftwf_plan_with_nthreads
    #define FFTW_PLAN_DFT_2D        fftwf_plan_dft_2d
    #define FFTW_PLAN_MANY_DFT      fftwf_plan_many_dft
    #define FFTW_PLAN_DFT_1D        fftwf_plan_dft_1d
    #define FFTW_EXECUTE_DFT        fftwf_execute_dft
    #define FFTW_DESTROY_PLAN fftwf_destroy_plan
    #define FFTW_MALLOC fftwf_malloc
    #define FFTW_FREE fftwf_free
    #define FFTW_IMPORT_WISDOM_FROM_FILENAME fftwf_import_wisdom_from_filename
    #define FFTW_EXPORT_WISDOM_TO_FILENAME   fftwf_export_wisdom_to_filename
    #define FFTW_EXPORT_WISDOM_TO_STRING     fftwf_export_wisdom_to_string
    #define FFTW_IMPORT_WISDOM_FROM_STRING   fftwf_import_wisdom_from_string
    #define FFTW_FORGET_WISDOM               fftwf_forget_wisdom
    
    // MPI
    #define MPI_COMPLEX_TYPE MPI_C_FLOAT_COMPLEX
    #define MPI_REAL_TYPE MPI_FLOAT
    
    #define fabs_t fabs
    #define fmax_t fmax
    #define sqrt_t sqrt
    
    // Constants
    #define PRECISION_NAME "Single (float)"
    #define BYTES_PER_COMPLEX 8
    
#endif

#endif

