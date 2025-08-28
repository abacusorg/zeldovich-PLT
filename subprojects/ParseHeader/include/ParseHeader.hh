#ifndef __PARSEHEADER_HH__
#define __PARSEHEADER_HH__

#include <stdio.h>
#include <string>
#include <filesystem>
#include <vector>
#include <unordered_map>

#include "detail/phDriver.hh"

namespace fs = std::filesystem;

#define MUST_DEFINE true
#define DONT_CARE false

class HeaderStream {
public:
    HeaderStream(const fs::path &fn);
    virtual ~HeaderStream(void);

    void OpenForRead(void);
    void Close(void);
    void SkipHeader(void);
    static size_t SkipHeaderFromFP(FILE *);
    void ReadHeader(void);

    void FinalizeHeader(void);

    fs::path name;
    char *buffer;
    size_t bufferlength;
    FILE *fp;

private:
    void GetHeaderLength(void);
};

void OpenStreamForWrite(std::ofstream& stream, const fs::path &fn, bool overwrite);
FILE* OpenForWrite(const fs::path &fn, bool overwrite);
void WriteHStream(FILE *fp, const std::string &m);
void WriteHStream(FILE *fp, const std::string &m, const std::string &pre);
void WriteHStream(FILE *fp, HeaderStream &in);
void WriteHStream(FILE *fp, HeaderStream &in, const std::string &pre);
void FinalizeHeader(FILE *fout);

class phDriver;

class ParseHeader {
public:
    ParseHeader(void);
    ~ParseHeader();

    // register variables with the parser
    template <typename T>
    void register_vars(T &param);

    // install a scalar
    template <typename T>
    void installscalar(const std::string &name, T& var, bool must_define);

    // Install a vector
    template <typename T>
    void installvector(const std::string &name, std::vector<T> &var, bool must_define, size_t maxlen = 1024);

    void ReadHeader(HeaderStream &in);

private:
    phDriver *phdriver;

    std::unordered_map<
        std::string, 
        std::variant<
            std::vector<int> *,
            std::vector<double> *,
            std::vector<std::string> *,
            std::vector<fs::path> *
        >
    > vectors;

    void resize_vectors(void);
};

template<typename T>
inline void ParseHeader::installscalar(const std::string &name, T &var, bool must_define) {
    phdriver->InstallSym(name, &var, 1, 0, must_define);
}

template <typename T>
inline void ParseHeader::installvector(const std::string &name, std::vector<T> &var, bool must_define, size_t maxlen) {
    var.resize(maxlen);
    phdriver->InstallSym(name, var.data(), maxlen, 1, must_define);
    vectors[name] = &var;  // will resize to the final length after parsing
}

#endif // __PARSEHEADER_HH__
