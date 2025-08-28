#ifndef __PHDRIVER_HH__
# define __PHDRIVER_HH__

#include <filesystem>
#include <string>
#include <map>
#include "detail/stringutil.hh"
#include "phParser.tab.hh"

#include <fmt/base.h>
#include <fmt/format.h>
#include <fmt/std.h>

namespace fs = std::filesystem;

template <> struct fmt::formatter<yy::location> : ostream_formatter {};

// Helper template to trigger static_assert for unsupported types
template<typename T>
struct always_false : std::false_type {};

// variable types for "type", below
enum class PHType { BAD, INTEGER, SIZE_T, FLOAT, DOUBLE, STRING, LOGICAL, LONG, PATH };

// types for values
enum class ValType { INTEGER, DOUBLE, STRING, LOGICAL };

// Tell Flex the lexer's prototype ...
# define YY_DECL                                   \
    yy::phParser::token_type                       \
    yylex (yy::phParser::semantic_type* yylval,    \
           yy::phParser::location_type* yylloc,    \
           phDriver& driver)
// ... and declare it for the parser's sake.
YY_DECL;

typedef long long int LLint;

// value stack -- these are the kinds of values which are
// assigned when the parser assembles a val_stack
typedef struct val_stack {
    ValType type;
    struct {
        unsigned u;
        LLint l;
        float f;
        double d;
        char *s;
    } value;
} VAL_STACK;

// this is the value assigned to a symbol table entry
typedef union {
    unsigned *u;
    int *i;
    LLint *l;
    size_t *z;
    float *f;
    double *d;
    std::string *s;
    fs::path *p;
} PVAL;

// parser symbol table entry
typedef struct {
    char *name;
    PHType type;
    PVAL pval;
    PVAL pvalorig;
    int stride;  // TODO: should stride always be 1 if we're using container types?
    int dim;
    int nvals;
    int *init;
    char* mapped;
} SYMENT;

typedef struct filestack {
  FILE *file_pointer;
  char *file_name;
  int line_number;
} FILE_STACK;

#define NSYMS 1000
#define MAX_STK 1000

// Conducting the whole scanning and parsing of Calc++.
class phDriver {
public:
    phDriver (int _Warnings, int _No_Undefined);
    virtual ~phDriver ();
     
    // external interface
    void ResetParser(void);

    template<typename T>
    void InstallSym(const std::string &name, T *ptr, int dim, int stride, bool init);
    int parse (const fs::path fname, bool Warn, int NoUndefined, bool stop);
    int parse (const char* buffer, int len, const fs::path fromfilename, 
                     bool Warn, int NoUndefined, bool stop);

    // Handling the scanner.
    void scan_begin(const fs::path &fn);
    void scan_begin(const char *buffer, int len);
    void scan_end ();

    // Error handling.
    void error (const yy::location& l, const std::string& m);
    void error (const std::string& m);

    void ERROR(std::string m, const yy::location& l);
    void PHWARNING(std::string m, const yy::location& l);
    void DEBUGOUT(std::string m, const yy::location& l);
    void MESSAGE(std::string m, const yy::location& l);

    // other
    void checkinit(void);
    void InitSymtab(void);    
    void DumpAllValues(void);

    // grammer utilies
    SYMENT *lookup(const char *id, int check);
    void TypeValError(PHType type, ValType val, const yy::location& l);
    void stuffit(SYMENT **spec, VAL_STACK *list, int len, int varzone, const yy::location& l);
        
    // grammer actions
    int InstallMapVar(char *val, const yy::location& yylloc);
    int InstallVCounter(char *val, const yy::location& yylloc);
    int ProcessEOS(const yy::location& yylloc);
    int ProcessVec_Start(const yy::location& yylloc);
    void ProcessValue_List(const yy::location& yylloc);
    int ProcessEquals(char *val, const yy::location& yylloc);
    void ReSetID_List(char *val);
    void SetID_List(char *val);
    void AddToValStack(LLint c);
    void AddToValStack(double c);
    void AddToValStack(bool c);
    void AddToValStack(char *c);
    void ReSetValStack(LLint c);
    void ReSetValStack(double c);
    void ReSetValStack(bool c);
    void ReSetValStack(char *c);

    // variables
    int fatal;                        // nonzero if a fatal error encountered
    int *nzones_prs;                  // number of zones in zone lists
    int *nzones_init;
    int nzprs_old;                    // number of zones in zone lists
    SYMENT *zonecounter_psym_old;
    int check_flag;
    int pp_len;
    bool Warnings;                     // flag governing whether we issue warnings
    int no_undefined;                 // flag governing whether undefined variables are ignored
    int init_symtab;                  // flag to initialize symtab
    bool StopOnError;                  // if true, don't keep parsing after fatal error

    SYMENT symtab[NSYMS];             // the symbol table 
    SYMENT *psym, *mapbase;
    char *id_stack[MAX_STK];
    VAL_STACK val_stack[MAX_STK];
    SYMENT *zonespec[MAX_STK];

    std::string filename;
    int trace_parsing;
    int trace_scanning;
    bool Debug;

    // pointer to the current parser
    yy::phParser *parser;

    int *yynerrs;
    int nwarn;

};

// install a variable in the symbol table and initialize its entry
template<typename T>
void phDriver::InstallSym(const std::string &name, T *ptr, int dim, int stride, bool init) {
    SYMENT *sym;

    if(init_symtab==1) InitSymtab();

    sym = lookup(name.c_str(),1);
    if(sym) {
        fmt::print(std::cerr, "symbol \"{}\" is already defined.\n", sym->name);
        exit(1);
    }
  
    sym = lookup(name.c_str(),0);
    if constexpr (std::is_same_v<T, int>) {
        if(Debug) fmt::print(std::cerr, "installing an int: \"{}\"\n", name);
        sym->pval.i = ptr;
        sym->pvalorig.i = ptr;
        sym->type = PHType::INTEGER;
    } else if constexpr (std::is_same_v<T, size_t>) {
        if(Debug) fmt::print(std::cerr, "installing a size_t: \"{}\"\n", name);
        sym->pval.z = ptr;
        sym->pvalorig.z = ptr;
        sym->type = PHType::SIZE_T;
    } else if constexpr (std::is_same_v<T, LLint>) {
        if(Debug) fmt::print(std::cerr, "installing a long long: \"{}\"\n", name);
        sym->pval.l = ptr;
        sym->pvalorig.l = ptr;
        sym->type = PHType::LONG;
    } else if constexpr (std::is_same_v<T, float>) {
        if(Debug) fmt::print(std::cerr, "installing a float: \"{}\"\n", name);
        sym->pval.f = (float *)ptr;
        sym->pvalorig.f = ptr;
        sym->type = PHType::FLOAT;
    } else if constexpr (std::is_same_v<T, double>) {
        if(Debug) fmt::print(std::cerr, "installing a double: \"{}\"\n", name);
        sym->pval.d = ptr;
        sym->pvalorig.d = ptr;
        sym->type = PHType::DOUBLE;
    } else if constexpr (std::is_same_v<T, std::string>) {
        if(Debug) fmt::print(std::cerr, "installing a string: \"{}\"\n", name);
        sym->pval.s = ptr;
        sym->pvalorig.s = ptr;
        sym->type = PHType::STRING;
    } else if constexpr (std::is_same_v<T, fs::path>) {
        if(Debug) fmt::print(std::cerr, "installing a path: \"{}\"\n", name);
        sym->pval.p = ptr;
        sym->pvalorig.p = ptr;
        sym->type = PHType::PATH;
    } else {
        static_assert(always_false<T>::value, "bad type\n");
    }

    sym->stride = stride;
    sym->dim = dim;
    if(init)                          // a false argument means don't care 
        *sym->init = -1;              // must be defined 
    else
        *sym->init = 1;               // don't care 
  
    sym->mapped = (char *)NULL;       // not mapped yet 
  
}


#endif //  __PHDRIVER_HH__
