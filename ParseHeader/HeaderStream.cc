#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "ParseHeader.hh"

HeaderStream::HeaderStream(std::string fn) {
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
        std::cerr << "HeaderStream::OpenForRead: filename is empty\n";
        exit(1);
    }
    if(fp != (FILE *)NULL) {
        std::cerr << "HeaderStream::OpenForRead: file is already open\n";
        exit(1);
    }
    fp=fopen(name.c_str(),"rb"); 
    if(fp==(FILE *) NULL) {
        std::cerr << "HeaderStream::OpenForRead:  cannot open filename \"" + name + "\"\n";
        exit(1);
    }
}

void HeaderStream::Close(void) {
    if(fp == (FILE *)NULL) {
        std::cerr << "HeaderStream::Close: file is already closed\n";
        exit(1);
    }
    fclose(fp);
    fp = (FILE *) NULL;
}

// returns length of header, including the ^B^B at end, and leave the
//     file open, at the end of the header.
void HeaderStream::SkipHeader(void) {
    OpenForRead();
    char buf[2];
    int len = 1;
    fread(&(buf[1]),1,1,fp);
    do {
        buf[0] = buf[1];
        if(fread(&(buf[1]), 1, 1, fp)<=0) goto error;
        len++;
    } while(!(buf[0]==0x2 && buf[1]=='\n'));
    bufferlength = len;
    return;

 error:
    std::cerr << "HeaderStream::SkipHeader: error in skipping header (no ^B\\n at end?)\n";
    fclose(fp);
    exit(1);
}

// get length of the header, including the ^B^B at the end
void HeaderStream::GetHeaderLength(void) {
    SkipHeader();
    Close();
}

void HeaderStream::ReadHeader(void) {
    GetHeaderLength();
    buffer = new char[bufferlength];
    OpenForRead();
    int len = 0;
    do {
        fread(&(buffer[len]), 1, 1, fp);
        len++;
    } while(!(buffer[len-2]==0x2 && buffer[len-1]=='\n'));
    assert(len == bufferlength);
    // replace the end-of-header token with two nulls as required by the parser
    buffer[len-2] = 0x0;
    buffer[len-1] = 0x0;
}

FILE* OpenForWrite(std::string fn, bool overwrite) {
    if(fn.empty()) {
        std::cerr << "OpenForWrite: filename is empty\n";
        exit(1);
    }

    FILE *fp=fopen(fn.c_str(),"rb"); 
    if(fp!=(FILE *) NULL) {
        if(!overwrite) {
            std::cerr << "OpenForWrite:  file \"" + fn +
                "\" exists and overwrite == false\n";
            exit(1);
        }
        else
            fclose(fp);
    }

    fp=fopen(fn.c_str(),"wb"); 
    if(fp==(FILE *) NULL) {
        std::cerr << "OpenForWrite:  cannot open filename \"" + fn + "\"\n";
        exit(1);
    }
    return fp;
}

void OpenStreamForWrite(std::ofstream& outfile, std::string fn, bool overwrite) {
    if(fn.empty()) {
        std::cerr << "OpenForWrite: filename is empty\n";
        exit(1);
    }

    outfile.open(fn.c_str());
    if(outfile.is_open()) {
        if(!overwrite) {
            std::cerr << "OpenForWrite:  file \"" + fn +
                "\" exists and overwrite == false\n";
            abort();
        }
        else
            outfile.close();
    }

    outfile.open(fn.c_str(), std::fstream::app);
    if(!outfile.is_open()) {
        std::cerr << "OpenForWrite:  cannot open filename \"" + fn + "\"\n";
        exit(1);
    }
}

void WriteHStream(FILE *fp, std::string m) {
    fwrite(m.c_str(), sizeof(char), m.length(), fp);
}

void WriteHStream(FILE *fp, std::string m, std::string pre) {
    std::string line(pre + m);
    fwrite(line.c_str(), sizeof(char), line.length(), fp);
}

void WriteHStream(FILE *fp, HeaderStream &in) {
    fwrite(in.buffer, sizeof(char), in.bufferlength-2, fp); // write all but the two terminal nulls
}

void WriteHStream(FILE *fp, HeaderStream &in, std::string pre) {
    std::stringstream ss(std::stringstream::in | std::stringstream::out);
    std::string line;
    ss << in.buffer;
    while(getline(ss,line)!=0) {
        std::cout << pre + line + "\n";
        line = pre + line + "\n";
        fwrite(line.c_str(), sizeof(char), line.length(), fp);
    }
}

void FinalizeHeader(FILE *fout) {
    if(fout == (FILE *) NULL) {
        std::cerr << "FinalizeHeader: file pointer is NULL\n";
        exit(1);
    }
    char tag[2] = {0x02, '\n'};
    fwrite(tag, sizeof(char), 2, fout);
}

