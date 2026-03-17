# NOTE: USE_DOUBLE_PRECISION can be set in src/config.h, but CFLAGS takes precedence.
#       If using double precision, also pass -DUSE_DOUBLE_PRECISION in CFLAGS
#       so the correct FFTW libraries are linked.

# Be explicit about compiler to match b/w fftw and omp
COMPILER_FAMILY ?= gcc
FFTW_DIR ?=

ifeq ($(COMPILER_FAMILY),gcc)
    MPICXX := $(shell which mpicxx 2>/dev/null)
    MPICC  := $(shell which mpicc 2>/dev/null)
    ifneq ($(MPICXX),)
        CXX = mpicxx
        CC  = mpicc
        MPI_INCLUDES =
        MPI_LIBS =
    else
        CXX = g++
        CC  = gcc
        MPI_INCLUDES =
        MPI_LIBS =
    endif
    OPENMP_FLAGS ?= -fopenmp
else ifeq ($(COMPILER_FAMILY),intel)
    MPIICPX := $(shell which mpiicpx 2>/dev/null)
    MPIICC  := $(shell which mpiicc 2>/dev/null)
    ifneq ($(MPIICPX),)
        CXX = mpiicpx
        CC  = mpiicc
        MPI_INCLUDES =
        MPI_LIBS =
    else
        CXX = icpx
        CC  = icx
        MPI_INCLUDES =
        MPI_LIBS =
    endif
    OPENMP_FLAGS ?= -qopenmp
else
    $(error Unknown COMPILER_FAMILY='$(COMPILER_FAMILY)'. Use gcc or intel)
endif

ifeq ($(FFTW_DIR),)
    $(error FFTW_DIR is not set. Example: make COMPILER_FAMILY=gcc FFTW_DIR=/path/to/fftw)
endif

FFTW_INCLUDES = -I$(FFTW_DIR)/include
FFTW_LIB_FLAG = -L$(FFTW_DIR)/lib

# User-provided flags
CFLAGS ?=

# ============================================================================
# BUILD FLAGS
# ============================================================================

BASE_CXXFLAGS = -Wall -Wextra -Wno-cast-function-type -std=c++17 -O3 $(OPENMP_FLAGS) -DFMT_HEADER_ONLY
ALL_CXXFLAGS = $(BASE_CXXFLAGS) $(CFLAGS)

INCLUDES = -Ideps/ParseHeader/include \
           -Ideps/ParseHeader/generated \
           -Ideps/zeldovich_core/include \
           -Ideps/zeldovich_core/pcg \
           -Ideps/fmt/include \
           -Isrc -Isrc/utils -Isrc/fft -Isrc/generation -Isrc/communication -Isrc/streaming -Isrc/output \
           -Ideps -I../.. -I../../../.. \
           $(MPI_INCLUDES) $(FFTW_INCLUDES)

# FFTW library selection
ifeq ($(findstring -DUSE_DOUBLE_PRECISION,$(CFLAGS)),)
    FFTW_OMP_LIB = $(FFTW_DIR)/lib/libfftw3f_omp.so
    FFTW_BASE_LIB = $(FFTW_DIR)/lib/libfftw3f.so
    FFTW_LIBS = $(FFTW_LIB_FLAG) -lfftw3f_omp -lfftw3f
else
    FFTW_OMP_LIB = $(FFTW_DIR)/lib/libfftw3_omp.so
    FFTW_BASE_LIB = $(FFTW_DIR)/lib/libfftw3.so
    FFTW_LIBS = $(FFTW_LIB_FLAG) -lfftw3_omp -lfftw3
endif

ifeq ($(wildcard $(FFTW_OMP_LIB)),)
    $(error Missing $(FFTW_OMP_LIB). FFTW must be built with OpenMP using the same compiler family)
endif

ifeq ($(wildcard $(FFTW_BASE_LIB)),)
    $(error Missing $(FFTW_BASE_LIB))
endif

BASE_LDFLAGS = $(MPI_LIBS) $(FFTW_LIBS) -lm -lstdc++
ifeq ($(findstring -fsanitize=address,$(CFLAGS)),)
    LDFLAGS = $(BASE_LDFLAGS)
else
    LDFLAGS = $(BASE_LDFLAGS) -fsanitize=address
endif

# ============================================================================
# SOURCE FILES
# ============================================================================

SRC = src/main.cpp

UTILS_SRC = src/utils/printing.c \
            src/utils/verification.c \
            src/utils/decomposition.c \
            src/utils/batch_helpers.c \
            src/utils/rng.c \
            src/utils/power_spectrum.c \
            src/utils/plt_eigenmodes.c

PARSEHEADER_SRC = \
    deps/ParseHeader/src/HeaderStream.cc \
    deps/ParseHeader/src/ParseHeader.cc \
    deps/ParseHeader/src/phDriver.cc \
    deps/ParseHeader/src/stringutil.cc \
    deps/ParseHeader/generated/phParser.tab.cc \
    deps/ParseHeader/generated/phScanner.cc

ZELDOVICH_CORE_SRC = \
    deps/zeldovich_core/src/power_spectrum.cpp \
    deps/zeldovich_core/src/parameters.cpp \
    deps/zeldovich_core/src/STimer.cc

ZELDOVICH_WRAPPER_SRC = src/utils/zeldovich_wrapper.cpp

MODULE_SRC = src/fft/fft_setup.c \
             src/generation/hermitian_generation.c \
             src/communication/mpi_exchange.c \
             src/streaming/z_streaming.c

OUTPUT_SRC = src/output/output_new.cpp

TARGET = hermitian_3d_matrix

# ============================================================================
# BUILD RULES
# ============================================================================

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC) $(UTILS_SRC) $(MODULE_SRC) $(ZELDOVICH_WRAPPER_SRC) $(OUTPUT_SRC) \
           $(PARSEHEADER_SRC) $(ZELDOVICH_CORE_SRC) src/config.h
	$(CXX) $(ALL_CXXFLAGS) $(INCLUDES) -o $(TARGET) \
        $(SRC) $(UTILS_SRC) $(MODULE_SRC) $(ZELDOVICH_WRAPPER_SRC) $(OUTPUT_SRC) \
        $(PARSEHEADER_SRC) $(ZELDOVICH_CORE_SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o