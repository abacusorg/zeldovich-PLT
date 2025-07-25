#ifndef __STRINGUTIL_H__
#define __STRINGUTIL_H__
#include <cstdlib>
#include <sstream>
#include <string>

// format anything which knows how to format itself into a std::string
template <typename T>
std::string ToString(T x) {
    std::stringstream ss(std::stringstream::in | std::stringstream::out);
    ss << x;
    return ss.str();
}

// lowercase a char*
void strnlwr(char *cstr, int len);

// Lowercase a std::string
void strlwr(std::string &str);

// get the root of a filename
std::string fileroot(std::string s);

// glob expand a filename
void get_absolute_pathname(std::string &name);

// following for doing printf-style formats
std::ostream &do_setformat(std::ostream &os, const char *fmt);

typedef struct fmtT {
    const char *format;
} fmtT;

std::basic_ostream<char> &operator<<(std::basic_ostream<char> &os, const fmtT &sep);

fmtT fmt(const char *t);

#endif  // __STRINGUTIL_H__
