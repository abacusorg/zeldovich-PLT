#include <iomanip>
#include <iostream>
#include <ios>
#include <wordexp.h>
#include "detail/stringutil.hh"

namespace stringutil {

// lowercase a char *
void strnlwr(char *cstr, int len) {
    for(int i=0; i<len; i++)
        if(cstr[i]>='A' && cstr[i]<='Z') cstr[i] = cstr[i] + 'a' - 'A';
}
// Lowercase a std::string
void strlwr(std::string &str) {
    for(std::string::iterator i=str.begin(); i<str.end(); i++) 
        if(*i>='A' && *i<='Z') *i = *i + 'a' - 'A';
}

std::ostream& do_setformat(std::ostream& os, const char *fmt) {
    int i = 0;
    while (fmt[i] != 0) {
        if (fmt[i] != '%') {  // pass through character strings not beginning with %
            os << fmt[i]; i++;
        }
        else {
            i++;
            if (fmt[i] == '%') {   // double % is escaped
                os << fmt[i]; i++;
            }
            else {
                bool ok = true;
                int istart = i;
                bool more = true;
                int width = 0;
                int precision = 6;
                std::ios::fmtflags flags = static_cast<std::ios::fmtflags>(0);
                char fill = ' ';
                bool alternate = false;
                while (more) {
                    switch (fmt[i]) {  
                    case '+':
                        flags |= std::ios::showpos;
                        break;
                    case '-':
                        flags |= std::ios::left;
                        break;
                    case '0':
                        flags |= std::ios::internal;
                        fill = '0';
                        break;
                    case '#':
                        alternate = true;
                        break;
                    case ' ':
                        flags |= std::ios::right;
                        break;
                    default:
                        more = false;
                        break;
                    }
                    if (more) i++;
                }
                if (isdigit(fmt[i])) {
                    width = atoi(fmt+i); 
                    do i++; while (isdigit(fmt[i]));
                }
                if (fmt[i] == '.') {
                    i++;
                    precision = atoi(fmt+i); 
                    while (isdigit(fmt[i])) i++;
                }
                switch (fmt[i]) {
                case 'd':
                    flags |= std::ios::dec | std::ios::right;
                    break;
                case 'x':
                    flags |= std::ios::hex;
                    if (alternate) flags |= std::ios::showbase;
                    break;
                case 'X':
                    flags |= std::ios::hex | std::ios::uppercase;
                    if (alternate) flags |= std::ios::showbase;
                    break;
                case 'o':
                    flags |= std::ios::hex;
                    if (alternate) flags |= std::ios::showbase;
                    break;
                case 'f':
                    flags |= std::ios::fixed;
                    if (alternate) flags |= std::ios::showpoint;
                    break;
                case 'e':
                    flags |= std::ios::scientific;
                    if (alternate) flags |= std::ios::showpoint;
                    break;
                case 'E':
                    flags |= std::ios::scientific | std::ios::uppercase;
                    if (alternate) flags |= std::ios::showpoint;
                    break;
                case 'g':
                    if (alternate) flags |= std::ios::showpoint;
                    break;
                case 'G':
                    flags |= std::ios::uppercase;
                    if (alternate) flags |= std::ios::showpoint;
                    break;
                default:
                    ok = false;
                    break;
                }
                i++;
                if (fmt[i] != 0) ok = false;
                if (ok) {
                    os.unsetf(static_cast<std::ios::fmtflags>(0xFFFFFFFF));
                    os.setf(flags);
                    os.width(width);
                    os.precision(precision);
                    os.fill(fill);
                }
                else i = istart;
            }
        }
    }
    return os;
}


std::basic_ostream<char>& operator<<(std::basic_ostream<char>& os, const fmtT& sep) {
    return do_setformat(os, sep.format);
}

fmtT fmt(const char *t) {
    fmtT sep;
    sep.format = t;
    return sep;
}

} // namespace stringutil


#ifdef TEST
int main(void) {

    std::string text("Hello There!");
    stringutil::strlwr(text);
    std::cout << text << "\n";

    double x = 17.45;
    double y = 43.55;

    std::cout << fmt("%7.1e") << x << y <<"\n";

    std::string name("~/file");
    get_absolute_pathname(name);
    std::cout << "The name was: " << name << "\n";


}
#endif
