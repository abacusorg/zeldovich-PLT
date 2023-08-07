# ParseHeader

A library for parsing plain text "key = value" headers. Used with Abacus inputs and outputs.

Original author: Phil Pinto

## Compiling
This is a meson project; to build, run:
```
$ meson setup build
$ meson compile -C build
```

## Building against ParseHeader
To use this library, include the `ParseHeader.h` header file in your source code, and link against the `libparseheader` library.

To use from another meson project, put this project in the `subprojects` directory and include it as a subproject in your `meson.build`:
```
libparseheader_dep = subproject('ParseHeader').get_variable('libparseheader_dep')
executable('myexe', 'myexe.c', dependencies: [libparseheader_dep])
```

## Usage
Directions for using the header classes:

```c++
#include "ParseHeader.hh"
```

Declare a typedef for each of the input files to be parsed containing
the variables one wishes to set. For example,

```c++
   typedef struct {
        double foo, goo;
        int zip[100];
        double *blah;
        char svec[10][1024];
   } codeparam;
```

Write a version of 

```c++
   template <>
   void ParseHeader::register_vars(codeparam &P) {
      installscalar("foo", P.foo, MUST_DEFINE);
      installscalar("goo", P.goo, MUST_DEFINE);
      installvector("blah", P.blah, 1024, 1, MUST_DEFINE); // can store a maximum of 1024 elements
      installvector("svec", P.svec, 10, 1024, 0);          // must be allocated contiguously   
   }
```

Then execute:

```c++
   codeparam cp;
   ParseHeader PH;
   PH.register_vars(cp);
   HeaderStream cpstream("myfile.in");
   PH.ReadHeader(cpstream);
```

## Future
The interpretation of ParseHeader headers is parser-dependent, which has led to some problems, especially between Python and C++. It would probably be healthier long-term to switch to a standardized format like TOML where the interpretation (tokenization, type inference, etc) is unambiguous.  It's also somewhat labor-intensive to define a new header field right now.  Also, there's no way to know (e.g. magic number) if a file has a ParseHeader header or not.

Also, the library can be hard to build on new platforms because of the bison/flex dependencies.  Furthermore, the code emits compiler warnings and triggers some sanitizer errors, and it's not clear if these are real bugs or not.
