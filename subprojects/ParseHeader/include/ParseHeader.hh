#ifndef __PARSEHEADER_HH__
#define __PARSEHEADER_HH__

#include <stdio.h>
#include <string>

#define MUST_DEFINE true
#define DONT_CARE false

class HeaderStream {
public:
    HeaderStream(std::string fn);
    virtual ~HeaderStream(void);

    void OpenForRead(void);
    void Close(void);
    void SkipHeader(void);
    void ReadHeader(void);

    void FinalizeHeader(void);

    std::string name;
    char *buffer;
    size_t bufferlength;
    FILE *fp;

private:
    void GetHeaderLength(void);
};

void OpenStreamForWrite(std::ofstream& stream, std::string fn, bool overwrite);
FILE* OpenForWrite(std::string fn, bool overwrite);
void WriteHStream(FILE *fp, std::string m);
void WriteHStream(FILE *fp, std::string m, std::string pre);
void WriteHStream(FILE *fp, HeaderStream &in);
void WriteHStream(FILE *fp, HeaderStream &in, std::string pre);
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
    void installscalar(std::string name, T var, bool must_define);

    void installscalar(std::string name, int &var, bool must_define);

    void installscalar(std::string name, float &var, bool must_define);

    void installscalar(std::string name, double &var, bool must_define);

    void installscalar(std::string name, long long int &var, bool must_define);

    void installscalar(std::string name, char *var, bool must_define);

// Install a vector
    template <typename T>
    void installvector(std::string name, T *var, int len, int stride, bool must_define);

    void installvector(std::string name, int *var, int len, int stride, bool must_define);

    void installvector(std::string name, float *var, int len, int stride, bool must_define);

    void installvector(std::string name, double *var, int len, int stride, bool must_define);

    void installvector(std::string name, long long int *var, int len, int stride, bool must_define);

    void installvector(std::string name, std::string *var, int len, int stride, bool must_define);

    void ReadHeader(HeaderStream &in);

private:
    phDriver *phdriver;
};
#endif // __PARSEHEADER_HH__
