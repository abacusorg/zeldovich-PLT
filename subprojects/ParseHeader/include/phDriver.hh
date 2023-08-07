#ifndef __PHDRIVER_HH__
# define __PHDRIVER_HH__

# include <string>
# include <map>
# include "stringutil.hh"
#include "phParser.tab.hh"

// variable types for "type", below
enum { BAD, INTEGER, FLOAT, DOUBLE, STRING, LOGICAL, LONG };

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
    int type;
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
    float *f;
    double *d;
    char *s;
} PVAL;

// parser symbol table entry
typedef struct {
    char *name;
    int type;
    PVAL pval;
    PVAL pvalorig;
    int stride;
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

// types for values
#define VALINTEGER 1
#define VALFLOAT 2
#define VALDOUBLE 3
#define VALSTRING 4
#define VALLOGICAL 5

// Conducting the whole scanning and parsing of Calc++.
class phDriver {
public:
    phDriver (int _Warnings, int _No_Undefined);
    virtual ~phDriver ();
     
    // external interface
    void ResetParser(void);
    void InstallSym(const char *name, void *ptr, int type, int dim, int stride, bool init);
    int parse (const std::string &fname, bool Warn, int NoUndefined, bool stop);
    int parse (const char* buffer, int len, std::string fromfilename, 
                     bool Warn, int NoUndefined, bool stop);

    // Handling the scanner.
    void scan_begin(std::string fn);
    void scan_begin(const char *buffer, int len);
    void scan_end ();

    // Error handling.
    void error (const yy::location& l, const std::string& m);
    void error (const std::string& m);

    void ERROR(std::string m, const yy::location& l);
    void WARNING(std::string m, const yy::location& l);
    void DEBUGOUT(std::string m, const yy::location& l);
    void MESSAGE(std::string m, const yy::location& l);

    // other
    void checkinit(void);
    void InitSymtab(void);    
    void DumpAllValues(void);

    // grammer utilies
    SYMENT *lookup(const char *id, int check);
    void TypeValError(int type, int val, const yy::location& l);
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

#endif //  __PHDRIVER_HH__
