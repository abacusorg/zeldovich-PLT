# NOTE: USE_DOUBLE_PRECISION can be set in src/config.h, but CFLAGS takes precedence.
#       The Makefile checks CFLAGS to select the correct FFTW libraries (single vs double).
#       If set USE_DOUBLE_PRECISION in config.h, should also add it to CFLAGS
#       to ensure the correct libraries are linked
# ====================================================================================

# Compiler and flags
# Use mpicxx if available (from frameworks module), otherwise fall back to g++
MPICXX := $(shell which mpicxx 2>/dev/null)
ifneq ($(MPICXX),)
    CXX = mpicxx
    CC = mpicc
    # mpicxx handles MPI includes and libraries automatically
    MPI_INCLUDES =
    MPI_LIBS =
else
    CXX = g++
    CC = gcc
    # MPI settings (OpenMPI) - fallback for systems without mpicxx
    MPI_INCLUDES = -I/usr/local/openmpi/4.1.6/gcc/include
    MPI_LIBS = -L/usr/local/openmpi/4.1.6/gcc/lib64 -lmpi -lmpi_cxx
endif

# OpenMP settings
# Prefer -qopenmp if compiler accepts it (Intel), else -fopenmp (GCC).
# More robust than parsing --version: mpicxx on Aurora may not print "icpx" in first line.
# Override: make OPENMP_FLAGS=-qopenmp  (or -fopenmp)
ifeq ($(OPENMP_FLAGS),)
    OPENMP_FLAGS := $(shell echo "int main(){return 0;}" | $(CXX) -x c++ - -c -qopenmp -o /dev/null >/dev/null 2>&1 && echo -qopenmp || echo -fopenmp)
endif
# Ensure -fopenmp/-qopenmp matches FFTW build (mismatch => threads may not spawn)

# FFTW3 settings
# CAVEAT: FFTW must be built with OpenMP support. libfftw3*_omp existence is a sign;
# ensure it was built with the same OpenMP runtime/flag as this project (e.g. Intel -qopenmp or -fopenmp)
# Use environment variables if available (from fftw module), otherwise fall back to hardcoded path
ifneq ($(C_INCLUDE_PATH),)
    # Extract FFTW path from C_INCLUDE_PATH (module system sets this)
    FFTW_INCLUDES = $(shell echo $(C_INCLUDE_PATH) | tr ':' '\n' | grep fftw | head -1 | sed 's|^|-I|')
    ifeq ($(FFTW_INCLUDES),)
        # Fallback: try to find in standard module locations
        FFTW_INCLUDES = $(shell find /opt/aurora -name "fftw3.h" 2>/dev/null | head -1 | xargs dirname | sed 's|^|-I|')
    endif
else
    # Fallback to hardcoded path
    FFTW_INCLUDES = -I/usr/local/fftw/gcc/openmpi-4.1.6/3.3.10/include
endif

# User-provided flags (override at command line)
CFLAGS ?=

# Base compilation flags
# v15.2: Updated to C++17 for zeldovich-PLT compatibility (ParseHeader requires C++17)
# Suppress -Wcast-function-type warnings from OpenMPI C++ bindings (known issue in OpenMPI 4.x)
BASE_CXXFLAGS = -Wall -Wextra -Wno-cast-function-type -std=c++17 -O3 $(OPENMP_FLAGS) -DFMT_HEADER_ONLY

# Address Sanitizer support (use with: make CFLAGS="-fsanitize=address -g -fno-omit-frame-pointer")
# Note: Address Sanitizer slows down execution significantly but detects heap corruption
# Usage: make clean && make CFLAGS="-fsanitize=address -g -fno-omit-frame-pointer"

# Combined flags
ALL_CXXFLAGS = $(BASE_CXXFLAGS) $(CFLAGS)

# Include paths (add src/ and module subdirs for headers)
# Note: zeldovich-PLT headers are in ../zeldovich-PLT/include
# v15.2: Added ParseHeader include path for zeldovich-PLT parameter file parsing
# NOTE: ParseHeader must be built with meson first to generate parser files
# Meson puts generated files in build/subprojects/ParseHeader/ (not subprojects/ParseHeader/build/)
# fmt is downloaded as a subproject by meson in subprojects/fmt-11.2.0/
# Include paths: zeldovich-PLT first to avoid conflicts with local headers
INCLUDES = -Ideps/ParseHeader/include \
           -Ideps/ParseHeader/generated \
           -Ideps/zeldovich_core/include \
           -Ideps/zeldovich_core/pcg \
           -Ideps/fmt/include \
           -Isrc -Isrc/utils -Isrc/fft -Isrc/generation -Isrc/communication -Isrc/streaming -Isrc/output \
           -Ideps -I../.. -I../../../.. \
           $(MPI_INCLUDES) $(FFTW_INCLUDES)

# FFTW library selection (default: single precision)
# Use environment variables if available (from fftw module), otherwise fall back to hardcoded path
ifneq ($(LD_LIBRARY_PATH),)
    # Extract FFTW lib path from LD_LIBRARY_PATH (module system sets this)
    FFTW_LIB_DIR = $(shell echo $(LD_LIBRARY_PATH) | tr ':' '\n' | grep fftw | head -1)
    ifeq ($(FFTW_LIB_DIR),)
        # Fallback: try to find in standard module locations
        FFTW_LIB_DIR = $(shell find /opt/aurora -name "libfftw3*.so" 2>/dev/null | head -1 | xargs dirname)
    endif
    ifneq ($(FFTW_LIB_DIR),)
        FFTW_LIB_FLAG = -L$(FFTW_LIB_DIR)
    else
        FFTW_LIB_FLAG = -L/usr/local/fftw/gcc/openmpi-4.1.6/3.3.10/lib
    endif
else
    FFTW_LIB_FLAG = -L/usr/local/fftw/gcc/openmpi-4.1.6/3.3.10/lib
endif

ifeq ($(findstring -DUSE_DOUBLE_PRECISION,$(CFLAGS)),)
    # Single precision (default): _omp for threading, base lib required (provides execute_dft etc.)
    FFTW_LIBS = $(FFTW_LIB_FLAG) -lfftw3f_omp -lfftw3f
else
    # Double precision: _omp for threading, base lib required
    FFTW_LIBS = $(FFTW_LIB_FLAG) -lfftw3_omp -lfftw3
endif

# Libraries
# v15.2: Added zeldovich-PLT libraries and RPATH for runtime library loading
# Address Sanitizer: Add -fsanitize=address to LDFLAGS if present in CFLAGS
BASE_LDFLAGS = $(MPI_LIBS) $(FFTW_LIBS) -lm -lstdc++
ifeq ($(findstring -fsanitize=address,$(CFLAGS)),)
    LDFLAGS = $(BASE_LDFLAGS)
else
    LDFLAGS = $(BASE_LDFLAGS) -fsanitize=address
endif

# Source files
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
             src/fft/fft_wisdom.c \
             src/generation/hermitian_generation.c \
             src/communication/mpi_exchange.c \
             src/streaming/z_streaming.c
OUTPUT_SRC = src/output/output_new.cpp
REASSEMBLY_SRC = src/write_particles_from_reassembled_mpi.cpp

# Output binaries
TARGET = mpi_zeldovich

# ====================================================================================
# BUILD RULES
# ====================================================================================

.PHONY: all clean

all: $(TARGET)

# Main executable (does not include reassembly tool due to main() conflict)
$(TARGET): $(SRC) $(UTILS_SRC) $(MODULE_SRC) $(ZELDOVICH_WRAPPER_SRC) $(OUTPUT_SRC) \
           $(PARSEHEADER_SRC) $(ZELDOVICH_CORE_SRC) src/config.h
	$(CXX) $(ALL_CXXFLAGS) $(INCLUDES) -o $(TARGET) \
        $(SRC) $(UTILS_SRC) $(MODULE_SRC) $(ZELDOVICH_WRAPPER_SRC) $(OUTPUT_SRC) \
        $(PARSEHEADER_SRC) $(ZELDOVICH_CORE_SRC) $(LDFLAGS)
clean:
	rm -f $(TARGET) *.o
