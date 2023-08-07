#include <string>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include "ParseHeader.hh"
#include "phDriver.hh"

ParseHeader::ParseHeader(void) {
    phdriver = new phDriver(1,1);
    phdriver->Debug = false;
}

ParseHeader::~ParseHeader(void) {
    delete phdriver;
}

template <typename T>
void ParseHeader::installscalar(std::string name, T var, bool must_define) {
    std::cerr << "ERROR: installscalar: " << name << " is not a recognized type\n";
    exit(1);
}

void ParseHeader::installscalar(std::string name, int &var, bool must_define) {
    phdriver->InstallSym(name.c_str(), &var, INTEGER, 1, 0, must_define);
}

void ParseHeader::installscalar(std::string name, float &var, bool must_define) {
    phdriver->InstallSym(name.c_str(), &var, FLOAT, 1, 0, must_define);
}

void ParseHeader::installscalar(std::string name, double &var, bool must_define) {
    phdriver->InstallSym(name.c_str(), &var, DOUBLE, 1, 0, must_define);
}

void ParseHeader::installscalar(std::string name, long long int &var, bool must_define) {
    phdriver->InstallSym(name.c_str(), &var, LONG, 1, 0, must_define);
}

void ParseHeader::installscalar(std::string name, char *var, bool must_define) {
    phdriver->InstallSym(name.c_str(), var, STRING, 1, 0, must_define);
}

template <typename T>
void ParseHeader::installvector(std::string name, T *var, int len, int stride, bool must_define) {
    std::cerr << "ERROR: installvector: " << name << " is not a recognized type\n";
    exit(1);
}

void ParseHeader::installvector(std::string name, int *var, int len, int stride, bool must_define) {
    phdriver->InstallSym(name.c_str(), var, INTEGER, len, stride, must_define);
}

void ParseHeader::installvector(std::string name, float *var, int len, int stride, bool must_define) {
    phdriver->InstallSym(name.c_str(), var, FLOAT, len, stride, must_define);
}

void ParseHeader::installvector(std::string name, double *var, int len, int stride, bool must_define) {
    phdriver->InstallSym(name.c_str(), var, DOUBLE, len, stride, must_define);
}

void ParseHeader::installvector(std::string name, long long int *var, int len, int stride, bool must_define) {
    phdriver->InstallSym(name.c_str(), var, LONG, len, stride, must_define);
}

template <>
void ParseHeader::installvector(std::string name, char **var, int len, int stride, bool must_define) {
    phdriver->InstallSym(name.c_str(), var[0], STRING, len, stride, must_define);
}

void ParseHeader::ReadHeader(HeaderStream &in) {
    in.ReadHeader();
    assert(in.buffer[in.bufferlength-1]==0x0 && in.buffer[in.bufferlength-2]==0x0);
    phdriver->trace_parsing = false;
    phdriver->trace_scanning = false;
    int result = phdriver->parse(in.buffer, in.bufferlength, in.name, 1, 0, false);
    if(result!=0) {
        std::cerr << "HS::parseit: there were errors parsing \"" << in.name << "\" ...exiting.\n";
        exit(1);
    }
    phdriver->ResetParser();
}
