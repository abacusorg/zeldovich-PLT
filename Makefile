CXX := $(shell which icc >/dev/null && echo icc || echo g++)
# Set DISK if you want to run the big BlockArray explicitly out of core.
# Set -DDIRECTIO and -I../Convolution if you want to use lib_dio
CXXFLAGS = -O3 -fopenmp -march=native -Wall -DDISK
INCL = -IParseHeader
LIBS = -LParseHeader -lparseheader -lfftw3 -lgsl -lgslcblas -lstdc++

all: zeldovich run_rng_test

zeldovich: zeldovich.o 
	make -C ParseHeader
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

zeldovich.o: zeldovich.cpp
	$(CXX) $(CXXFLAGS) $(INCL) -c $^

rng_test: rng_test.c
	$(CXX) $(CXXFLAGS) $(INCL) $^ -o $@ $(LIBS)

run_rng_test: rng_test
	@./rng_test | cmp rng_test.out - && (echo 'Passed RNG test.') || (echo 'Error: your platform did not produce the expected RNG values, and may thus generate IC files with unexpected phases.' ; exit 1)

default: zeldovich

.PHONY: all clean distclean run_rng_test default
clean:
	make -C ParseHeader $@
	$(RM) *.o *.gch *~
distclean:
	make -C ParseHeader $@
	$(RM) *.o *.gch zeldovich rng_test *~
