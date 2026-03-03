# It uses zeldovich-PLT files - need to link properly
# It also builds the reassembly tool to write the particles to a file.
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
# Use -qopenmp for Intel compiler (mpicxx on Aurora), -fopenmp for GCC
ifeq ($(findstring icpx,$(shell $(CXX) --version 2>&1 | head -1)),)
    OPENMP_FLAGS = -fopenmp
else
    OPENMP_FLAGS = -qopenmp
endif

# FFTW3 settings
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
BASE_CXXFLAGS = -Wall -Wextra -Wno-cast-function-type -std=c++17 -O3 $(OPENMP_FLAGS)

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
INCLUDES = -I../zeldovich-PLT/include \
           -I../zeldovich-PLT/subprojects/ParseHeader/include \
           -I../zeldovich-PLT/build/subprojects/ParseHeader \
           -I../zeldovich-PLT/subprojects/fmt-11.2.0/include \
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
    # Single precision (default)
    FFTW_LIBS = $(FFTW_LIB_FLAG) -lfftw3f -lfftw3f_omp
else
    # Double precision
    FFTW_LIBS = $(FFTW_LIB_FLAG) -lfftw3 -lfftw3_omp
endif

# zeldovich-PLT libraries (v15.2: Phase 7 - linking)
# Note: zeldovich-PLT may have been built with GSL support
# If linking fails with GSL errors, rebuild zeldovich-PLT on this system
ZELDOVICH_LIB_DIRS = -L../zeldovich-PLT/build \
                     -L../zeldovich-PLT/build/subprojects/ParseHeader \
                     -L../zeldovich-PLT/build/subprojects/fmt-11.2.0
# Try to find GSL libraries (zeldovich-PLT may depend on them)
GSL_LIB_DIR = $(shell find /usr/lib* /opt/aurora -name "libgsl.so*" 2>/dev/null | head -1 | xargs dirname 2>/dev/null)
ifneq ($(GSL_LIB_DIR),)
    GSL_LIBS = -L$(GSL_LIB_DIR) -lgsl -lgslcblas
else
    # GSL not found - may need to rebuild zeldovich-PLT without GSL or install GSL
    GSL_LIBS =
endif
ZELDOVICH_LIBS = -lzeldovich -lparseheader -lfmt $(GSL_LIBS)

# Libraries
# v15.2: Added zeldovich-PLT libraries and RPATH for runtime library loading
# Address Sanitizer: Add -fsanitize=address to LDFLAGS if present in CFLAGS
BASE_LDFLAGS = $(MPI_LIBS) $(FFTW_LIBS) $(ZELDOVICH_LIB_DIRS) $(ZELDOVICH_LIBS) -lm -lstdc++ \
               -Wl,-rpath,$(shell cd ../zeldovich-PLT/build && pwd):$(shell cd ../zeldovich-PLT/build/subprojects/ParseHeader && pwd):$(shell cd ../zeldovich-PLT/build/subprojects/fmt-11.2.0 && pwd)
ifeq ($(findstring -fsanitize=address,$(CFLAGS)),)
    LDFLAGS = $(BASE_LDFLAGS)
else
    LDFLAGS = $(BASE_LDFLAGS) -fsanitize=address
endif

# External dependencies
STIMER_CC = src/STimer.cc

# Source files
SRC = src/main.cpp
UTILS_SRC = src/utils/printing.c \
            src/utils/verification.c \
            src/utils/decomposition.c \
            src/utils/batch_helpers.c \
            src/utils/rng.c \
            src/utils/power_spectrum.c \
            src/utils/plt_eigenmodes.c

# C++ wrapper for zeldovich-PLT (Option B: C wrappers)
# v15.2: Enabled for direct zeldovich-PLT integration
# TODO: Replace with library linking when Option A is implemented
ZELDOVICH_WRAPPER_SRC = src/utils/zeldovich_wrapper.cpp
MODULE_SRC = src/fft/fft_setup.c \
             src/generation/hermitian_generation.c \
             src/communication/mpi_exchange.c \
             src/streaming/z_streaming.c
OUTPUT_SRC = src/output/output_new.cpp
REASSEMBLY_SRC = src/write_particles_from_reassembled_mpi.cpp

# Output binaries
TARGET = hermitian_3d_matrix
REASSEMBLY_TARGET = write_particles_from_reassembled_mpi

# ====================================================================================
# BUILD RULES
# ====================================================================================

.PHONY: all clean reassembly

all: $(TARGET)

# Main executable (does not include reassembly tool due to main() conflict)
$(TARGET): $(SRC) $(UTILS_SRC) $(MODULE_SRC) $(STIMER_CC) $(ZELDOVICH_WRAPPER_SRC) $(OUTPUT_SRC) src/config.h
	$(CXX) $(ALL_CXXFLAGS) $(INCLUDES) -o $(TARGET) $(SRC) $(UTILS_SRC) $(MODULE_SRC) $(STIMER_CC) $(ZELDOVICH_WRAPPER_SRC) $(OUTPUT_SRC) $(LDFLAGS)

# Reassembly tool (separate executable)
reassembly: $(REASSEMBLY_TARGET)

# Reassembly uses zeldovich-PLT's STimer (via <output.h> and -lzeldovich); do not link our STimer.cc to avoid redefinition
$(REASSEMBLY_TARGET): $(REASSEMBLY_SRC) $(UTILS_SRC) $(ZELDOVICH_WRAPPER_SRC) $(OUTPUT_SRC) src/config.h
	$(CXX) $(ALL_CXXFLAGS) $(INCLUDES) -o $(REASSEMBLY_TARGET) $(REASSEMBLY_SRC) $(UTILS_SRC) $(ZELDOVICH_WRAPPER_SRC) $(OUTPUT_SRC) $(LDFLAGS)
	@echo ""
	@echo "Build successful!"
	@echo "Binary: $(TARGET)"
	@echo "Configuration: see src/config.h"
	@echo ""

clean:
	rm -f $(TARGET) $(REASSEMBLY_TARGET) *.o

# ====================================================================================
# HELP
# ====================================================================================

help:
	@echo "Hermitian 3D Matrix MPI - Build System"
	@echo ""
	@echo "Targets:"
	@echo "  make              - Build main executable (hermitian_3d_matrix)"
	@echo "  make reassembly   - Build reassembly tool (write_particles_from_reassembled_mpi)"
	@echo "  make clean        - Remove build artifacts"
	@echo "  make help         - Show this help message"
	@echo ""
	@echo "Configuration:"
	@echo "  Edit src/config.h or use CFLAGS to override"
	@echo ""
	@echo "Examples:"
	@echo "  make CFLAGS=\"-DUSE_DOUBLE_PRECISION\""
	@echo "  make CFLAGS=\"-DPRODUCTION_MODE\""
	@echo "  make CFLAGS=\"-DUSE_X_PADDING=0 -DDEBUG_PRINTS=0\""
	@echo ""

