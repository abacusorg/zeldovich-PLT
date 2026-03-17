#include <string>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include "ParseHeader.hh"
#include "detail/phDriver.hh"

#include <fmt/base.h>
#include <fmt/std.h>
#include <fmt/ostream.h>

ParseHeader::ParseHeader(void) {
    phdriver = new phDriver(1,1);
    phdriver->Debug = false;
}

ParseHeader::~ParseHeader(void) {
    delete phdriver;
}

void ParseHeader::ReadHeader(HeaderStream &in) {
    in.ReadHeader();
    assert(in.buffer[in.bufferlength-1]==0x0 && in.buffer[in.bufferlength-2]==0x0);
    phdriver->trace_parsing = false;
    phdriver->trace_scanning = false;
    int result = phdriver->parse(in.buffer, in.bufferlength, in.name, 1, 0, false);
    if(result!=0) {
        fmt::print(std::cerr, "HS::parseit: there were errors parsing \"{}\" ...exiting.\n", in.name);
        exit(1);
    }
    resize_vectors();
    phdriver->ResetParser();
}

void ParseHeader::resize_vectors(void) {
    // We over-allocated the user's vectors so that ParseHeader could read into the
    // underlying buffers, C-style. Now we need to resize the vectors to the used length.
    for(auto &v : vectors) {
        SYMENT *sym = phdriver->lookup(v.first.c_str(), 1);
        if(sym) {
            std::visit([n=sym->nvals](auto&& arg){ arg->resize(n); }, v.second);
        } else {
            fmt::print(std::cerr, "ParseHeader::resize_vectors: symbol \"{}\" not found.\n", v.first);
            exit(1);
        }
    }
}
