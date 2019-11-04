# Load Abacus's make settings if available
ROOT_DIR := ..
-include $(ROOT_DIR)/common.mk

CXX ?= g++
# Set DISK if you want to run the big BlockArray explicitly out of core.
# Set -DDIRECTIO and -I../Convolution if you want to use lib_dio
CXXFLAGS ?= -O3 -fopenmp -march=native -std=c++11
CXXFLAGS += -Wall
PARSEHEADER_CPPFLAGS ?= -I ParseHeader
PARSEHEADER_LIBS ?= -L ParseHeader -lparseheader

GSL_LIBS ?= $(shell gsl-config --libs)
GSL_CPPFLAGS ?= $(shell gsl-config --cflags)

CPPFLAGS +=$(TCMALLOC_CPPFLAGS) $(THREAD_CPPFLAGS) $(PARSEHEADER_CPPFLAGS) $(GSL_CPPFLAGS) -DDISK

LIBS += $(TCMALLOC_LIBS) $(PARSEHEADER_LIBS) -lfftw3 $(GSL_LIBS) -lstdc++

all: zeldovich 

zeldovich: zeldovich.o | ParseHeader
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ $(LIBS)

zeldovich.o: zeldovich.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MMD -c $<

-include zeldovich.d

.PHONY: all clean distclean default
clean:
	$(MAKE) -C ParseHeader $@
	$(RM) *.o *.gch *~

distclean: clean
	$(MAKE) -C ParseHeader $@
	$(RM) zeldovich *.d

ifndef HAVE_COMMON_MK
ParseHeader: ParseHeader/libparseheader.a
ParseHeader/libparseheader.a:
	$(MAKE) -C ParseHeader libparseheader.a
endif
