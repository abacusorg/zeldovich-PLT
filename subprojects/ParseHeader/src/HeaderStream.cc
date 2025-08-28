#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "ParseHeader.hh"

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/std.h>

HeaderStream::HeaderStream(const fs::path &fn) {
    name = fn;
    buffer = (char *) NULL;
    bufferlength = 0;
    fp = (FILE *) NULL;
}

HeaderStream::~HeaderStream(void) {
    if(buffer != (char *) NULL) delete[] buffer;
}

void HeaderStream::OpenForRead(void) {
    if(name.empty()) {
        fmt::print(std::cerr, "HeaderStream::OpenForRead: filename is empty\n");
        exit(1);
    }
    if(fp != (FILE *)NULL) {
        fmt::print(std::cerr, "HeaderStream::OpenForRead: file is already open\n");
        exit(1);
    }
    fp=fopen(name.c_str(),"rb"); 
    if(fp==(FILE *) NULL) {
        fmt::print(std::cerr, "HeaderStream::OpenForRead:  cannot open filename \"{}\"\n", name);
        exit(1);
    }
}

void HeaderStream::Close(void) {
    if(fp == (FILE *)NULL) {
        fmt::print(std::cerr, "HeaderStream::Close: file is already closed\n");
        exit(1);
    }
    fclose(fp);
    fp = (FILE *) NULL;
}

// returns length of header, including the ^B^B at end, and leave the
//     file open, at the end of the header.
void HeaderStream::SkipHeader(void) {
    OpenForRead();
    bufferlength = HeaderStream::SkipHeaderFromFP(fp);
    return;
}

// returns length of header, including the ^B^B at end, and leave the
//     file open, at the end of the header.
size_t HeaderStream::SkipHeaderFromFP(FILE *fp) {
    size_t bufferlength;
    char buf[2];
    if(fread(&(buf[1]), 1, 1, fp)<=0){
        bufferlength = 2;  // pseudo ^B^B
        return bufferlength;
    }
    size_t len = 1;
    do {
        buf[0] = buf[1];
        if(fread(&(buf[1]), 1, 1, fp)<=0){
            len += 2; // pseudo ^B^B
            break;
        }
        len++;
    } while(!(buf[0]==0x2 && buf[1]=='\n'));
    
    bufferlength = len;
    return bufferlength;
}

// get length of the header, including the ^B^B at the end
void HeaderStream::GetHeaderLength(void) {
    SkipHeader();
    Close();
}

void HeaderStream::ReadHeader(void) {
    GetHeaderLength();  // this guarantees there is a valid header
    buffer = new char[bufferlength];
    OpenForRead();
    size_t nread = fread(&(buffer[0]), 1, bufferlength-2, fp);
    assert(nread == bufferlength-2);
    // replace the end-of-header token with two nulls as required by the parser
    buffer[bufferlength-2] = 0x0;
    buffer[bufferlength-1] = 0x0;
}

FILE* OpenForWrite(const fs::path &fn, bool overwrite) {
    if(fn.empty()) {
        fmt::print(std::cerr, "OpenForWrite: filename is empty\n");
        exit(1);
    }

    FILE *fp=fopen(fn.c_str(),"rb"); 
    if(fp!=(FILE *) NULL) {
        if(!overwrite) {
            fmt::print(std::cerr, "OpenForWrite:  file \"{}\" exists and overwrite == false\n", fn);
            exit(1);
        }
        else
            fclose(fp);
    }

    fp=fopen(fn.c_str(),"wb"); 
    if(fp==(FILE *) NULL) {
        fmt::print(std::cerr, "OpenForWrite:  cannot open filename \"{}\"\n", fn);
        exit(1);
    }
    return fp;
}

void OpenStreamForWrite(std::ofstream& outfile, const fs::path &fn, bool overwrite) {
    if(fs::exists(fn)) {
        if(!overwrite) {
            fmt::print(std::cerr, "OpenForWrite:  file \"{}\" exists and overwrite == false\n", fn);
            abort();
        }
        else
            outfile.close();
    }

    outfile.open(fn, std::fstream::app);
    if(!outfile.is_open()) {
        fmt::print(std::cerr, "OpenForWrite:  cannot open filename \"{}\"\n", fn);
        exit(1);
    }
}

void WriteHStream(FILE *fp, const std::string &m) {
    fmt::print(fp, m);
}

void WriteHStream(FILE *fp, const std::string &m, const std::string &pre) {
    fmt::print(fp, "{}{}", pre, m);
}

void WriteHStream(FILE *fp, HeaderStream &in) {
    fmt::print(fp, "{}", fmt::string_view(in.buffer, in.bufferlength-2));
}

void WriteHStream(FILE *fp, HeaderStream &in, const std::string &pre) {
    std::stringstream ss(std::stringstream::in | std::stringstream::out);
    std::string line;
    ss << in.buffer;
    while(getline(ss,line)) {
        std::string s = fmt::format("{}{}\n", pre, line);
        fmt::print(std::cout, "{}", s);  // TODO: why are we echoing to stdout?
        fmt::print(fp, "{}", s);
    }
}

void FinalizeHeader(FILE *fout) {
    if(fout == (FILE *) NULL) {
        fmt::print(std::cerr, "FinalizeHeader: file pointer is NULL\n");
        exit(1);
    }
    char tag[2] = {0x02, '\n'};
    fmt::print(fout, "{}", fmt::string_view(tag,2));
}
