# Load Abacus's make settings if available
ROOT_DIR := ..
-include $(ROOT_DIR)/common.mk

CXX ?= g++
# Set DISK if you want to run the big BlockArray explicitly out of core.
# Set -DDIRECTIO and -I../Convolution if you want to use lib_dio
CXXFLAGS ?= -O3 -fopenmp -march=native 
CXXFLAGS += -Wall
PARSEHEADER_CPPFLAGS ?= -I ParseHeader
PARSEHEADER_LIBS ?= -L ParseHeader -lparseheader

GSL_LIBS ?= $(shell gsl-config --libs)
GSL_CPPFLAGS ?= $(shell gsl-config --cflags)

CPPFLAGS += -DDISK $(TCMALLOC_CPPFLAGS) $(THREAD_CPPFLAGS) $(PARSEHEADER_CPPFLAGS) $(GSL_CPPFLAGS)

LIBS += $(TCMALLOC_LIBS) $(PARSEHEADER_LIBS) -lfftw3 $(GSL_LIBS) -lstdc++

all: zeldovich 

zeldovich: zeldovich.o | ParseHeader
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ $(LIBS)

zeldovich.o: zeldovich.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MMD -c $<

-include zeldovich.d

rng_test: rng_test.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $^ -o $@ $(LIBS)

run_rng_test: rng_test
	@./rng_test | cmp rng_test.out - && (echo 'Passed RNG test.') || (echo 'Error: your platform did not produce the expected RNG values, and may thus generate IC files with unexpected phases.' ; exit 1)

.PHONY: all clean distclean run_rng_test default
clean:
	$(MAKE) -C ParseHeader $@
	$(RM) *.o *.gch *~

distclean: clean
	$(MAKE) -C ParseHeader $@
	$(RM) zeldovich rng_test *.d

ifndef HAVE_COMMON_MK
ParseHeader: ParseHeader/libparseheader.a
ParseHeader/libparseheader.a:
	$(MAKE) -C ParseHeader libparseheader.a
endif