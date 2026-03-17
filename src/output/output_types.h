#pragma once

#include <cstdint>

class ZelParticle {
public:
    unsigned short i, j, k;
    double displ[3];
};

class ZelSimpleParticle {
public:
    float displ[3];
};

class RVZelParticle {
public:
    unsigned short i, j, k;
    float displ[3];
    float vel[3];
};

class RVdoubleZelParticle {
public:
    unsigned short i, j, k;
    double displ[3];
    double vel[3];
};

enum OutputType {
    OUTPUT_ZEL,
    OUTPUT_RVZEL,
    OUTPUT_RVDOUBLEZEL,
    OUTPUT_ZEL_SIMPLE
};