#include <string.h>
#include <sstream>
#include <cmath>
//#include <cfloat>
#include <limits>
#include "phScanner.hh"
#include "detail/phDriver.hh"
#include "detail/stringutil.hh"

#include <fmt/base.h>
#include <fmt/format.h>
#include <fmt/std.h>

using std::string;
using std::endl;

phDriver::phDriver (int _Warnings, int _No_Undefined) 
    : Warnings(_Warnings),
      no_undefined(_No_Undefined),
      trace_parsing (false),
      trace_scanning (false),
      Debug(false) {
    InitSymtab();
    ResetParser();
}
     
phDriver::~phDriver () { }

int phDriver::parse(const fs::path fn, bool Warn, int NoUndefined, bool stop) {
    Warnings = Warn;
    no_undefined = NoUndefined;
    filename = fn;
    StopOnError = stop;

    scan_begin(fn);
    yy::phParser parser (*this);
    parser.set_debug_level (trace_parsing);
    int res = parser.parse ();
    if(*yynerrs>0) res = 1;
    if(nwarn>0) error("there were " + stringutil::ToString(nwarn) + " warnings\n");
    scan_end ();

    checkinit();

    return res;
}

int phDriver::parse(const char* buffer, int len, const fs::path fromfilename, 
                    bool Warn, int NoUndefined, bool stop) {
    Warnings = Warn;
    no_undefined = NoUndefined;
    filename = fromfilename;
    StopOnError = stop;

    scan_begin(buffer, len);
    yy::phParser parser (*this);
    parser.set_debug_level (trace_parsing);
    int res = parser.parse ();
    if(*yynerrs>0) res = 1;
    if(nwarn>0) error("there were " + stringutil::ToString(nwarn) + " warnings\n");
//    scan_end (); // no file to close!
  
    checkinit();

    return res;
}

void phDriver::error (const yy::location& l, const std::string& m) {
    fmt::print(std::cerr, "ParseHeader: {}: {}\n", l, m);
}
     
void phDriver::error (const std::string& m) {
    fmt::print(std::cerr, "ParseHeader: {}\n", m);
}

void phDriver::MESSAGE(std::string m, const yy::location& l) {
    if(Debug) {
        error(l, m);
    }
}

void phDriver::DEBUGOUT(std::string m, const yy::location& l) {
    if(Debug) {
        error(l, std::string("debug: ") + m);
    }
}

void phDriver::PHWARNING(std::string m, const yy::location& l) {
    if(Warnings) {
        nwarn++;
        error(l, std::string("Warning: ") + m);
    }
}

void phDriver::ERROR(std::string m, const yy::location& l) {
    (*yynerrs)++;
    error(l, std::string("ERROR ") + stringutil::ToString((*yynerrs))+ string(": ") + m);
    if(StopOnError) {
        error(string("there were errors and StopOnError was true; exiting..."));
        exit(1);
    }
}

void phDriver::scan_begin(const fs::path &fn) {
    extern int yy_flex_debug;
    yy_flex_debug = trace_scanning;
    if (fn == "-") 
        yyin = stdin;
    else if (!(yyin = fopen (fn.c_str (), "r")))
    {
        error (std::string ("scan_begin: cannot open file:\"") + filename + string("\"."));
        exit (1);
    }
}

void phDriver::scan_begin(const char *buffer, int len) {
    extern int yy_flex_debug;
    yy_flex_debug = trace_scanning;
    if(buffer == (char *) NULL) {
        error (std::string ("scan_begin: buffer is null."));
        exit (1);
    }
    if(len <= 0 ) {
        error (std::string ("scan_begin: buffer length <= 0."));
        exit (1);
    }
    yy_switch_to_buffer( yy_scan_buffer(const_cast<char *>(buffer), len) );
}

void phDriver::scan_end() {
    fclose (yyin);
}
   
// look up an id in the symbol table, adding one if necessary;
//   return a pointer to the entry if there
SYMENT *phDriver::lookup(const char *id, int check) {
    int i;
  
    for (i=0; i<NSYMS; i++) {
        if (symtab[i].name==0) goto newentry;
        if (strcmp(id,symtab[i].name)==0) return (&symtab[i]);
    }
  
    error(std::string("symbol table is full. Recompile with larger NSYMS."));
    exit(1);
  
newentry:
    if(check==0) {
        symtab[i].name = (char *) malloc(strlen(id)+1);
        symtab[i].init = (int *) malloc(sizeof(int));
        strcpy(symtab[i].name,id);
        symtab[i].pval.d = (double *) NULL;
        symtab[i].pvalorig.d = (double *) NULL;
        *symtab[i].init = 1;
        symtab[i].type = PHType::BAD;                   // don't care 
        return (&symtab[i]);
    }
    else
        return(NULL);
}

void phDriver::TypeValError(PHType type, ValType val, const yy::location& yylloc) {
    std::string typeStr;
    std::string valStr;

    switch (type) {
        case PHType::LONG: typeStr = "LONG"; break;
        case PHType::DOUBLE: typeStr = "DOUBLE"; break;
        case PHType::FLOAT: typeStr = "FLOAT"; break;
        case PHType::INTEGER: typeStr = "INTEGER"; break;
        case PHType::SIZE_T: typeStr = "SIZE_T"; break;
        case PHType::LOGICAL: typeStr = "LOGICAL"; break;
        case PHType::STRING: typeStr = "STRING"; break;
        case PHType::PATH: typeStr = "PATH"; break;
        case PHType::BAD: typeStr = "BAD"; break;
    }

    switch (val) {
        case ValType::INTEGER: valStr = "INTEGER"; break;
        case ValType::DOUBLE: valStr = "DOUBLE"; break;
        case ValType::LOGICAL: valStr = "LOGICAL"; break;
        case ValType::STRING: valStr = "STRING"; break;
    }

    ERROR("attempt to set variable of type " + typeStr + " to value of type " + valStr, yylloc);
}

#define STUFF(type, varsym, valsym, stride)                             \
    if(pf->mapped == (char *) NULL) {                                   \
        *pf->pval.varsym = (type) list[i].value.valsym;                 \
    }                                                                   \
    else {                                                              \
        *pf->pval.varsym += (type) list[i].value.valsym;                \
    }                                                                   \
    pf->pval.varsym += pf->stride;                                      \
    pf->nvals++;                                                        \
    if(pf->nvals>pf->dim) {                                             \
        ERROR(string("number of values (") +                            \
          stringutil::ToString(pf->nvals) + string(") exceeds dimension for ") +      \
          string(pf->name) + string("[") + stringutil::ToString(pf->dim) +            \
              string("]"), yylloc);                                     \
    }                                                                   \
    *pf->init = 1;                                                      \
    
// Put the values in list into the array pf
//   (a candidate for template metaprogramming if ever there was one...)
void phDriver::stuffit(SYMENT **spec, VAL_STACK *list, int len, int varzone,
                       const yy::location& yylloc) {
    int i;
    SYMENT *pf;

    for(i=0; i<len; i++) {
        if(varzone) 
            pf = spec[i];
        else
            pf = spec[0];


        if(pf->pval.d) {
            switch (list[i].type) {
        
            case ValType::INTEGER:
                switch (pf->type) {
                case PHType::LONG:
                    STUFF(LLint, l, l, stride);
                    break;
                case PHType::DOUBLE:
                    STUFF(double, d, l, stride);
                    break;
                case PHType::FLOAT:
                    STUFF(float, f, l, stride);
                    break;
                case PHType::INTEGER:
                    STUFF(int, i, l, stride);
                    if(llabs(list[i].value.l)>std::numeric_limits<int>::max())
                        ERROR(string("attempt to store too large a value: ") + 
                              stringutil::ToString(list[i].value.l) + string(" in an int variable: ") +
                              string(pf->name), yylloc);
                    break;
                case PHType::SIZE_T:
                    STUFF(size_t, z, l, stride);
                    if(list[i].value.l<0)
                        ERROR(string("attempt to store a negative value: ") + 
                              stringutil::ToString(list[i].value.l) + string(" in a size_t variable: ") +
                              string(pf->name), yylloc);
                    break;
                case PHType::LOGICAL:
                case PHType::STRING:
                case PHType::PATH:
                    TypeValError(pf->type,list[i].type, yylloc);
                    break;
                case PHType::BAD:
                    ERROR(string("stuffit: logic failure!"), yylloc);
                    exit(1);
                }
                break;
        
            case ValType::DOUBLE:
                switch (pf->type) {
                case PHType::LONG:
                    STUFF(LLint, l, d, stride);
                    PHWARNING(string("truncating a float: ") + stringutil::ToString(list[i].value.d) +
                            string(" to a long long int for \"") + string(pf->name) + 
                            string("\"."), yylloc);
                    break;
                case PHType::DOUBLE:
                    STUFF(double, d, d, stride);
                    break;
                case PHType::FLOAT:
                    STUFF(float, f, d, stride);
                    if(fabs(list[i].value.d)>std::numeric_limits<float>::max()) {
                        ERROR(string("attempt to store too large a value: ") +
                              stringutil::ToString(list[i].value.d) + string(" in an float variable: ") +
                              string(pf->name), yylloc);
                    }
                    break;
                case PHType::INTEGER:
                    STUFF(int, i, d, stride);
                    PHWARNING(string("truncating a float: ") + stringutil::ToString(list[i].value.d) + 
                            string(" to an int for \"") + string(pf->name) + string("\"."), yylloc);
                    break;
                case PHType::SIZE_T:
                    STUFF(size_t, z, d, stride);
                    if(list[i].value.d<0)
                        ERROR(string("attempt to store a negative value: ") + 
                              stringutil::ToString(list[i].value.d) + string(" in a size_t variable: ") +
                              string(pf->name), yylloc);
                    break;
                case PHType::LOGICAL:
                case PHType::STRING:
                case PHType::PATH:
                    TypeValError(pf->type,list[i].type, yylloc);
                    break;
                case PHType::BAD:
                    ERROR(string("stuffit: logic failure!"), yylloc);
                    exit(1);
                }
                break;
        
            case ValType::LOGICAL:
                switch (pf->type) {
                case PHType::LOGICAL:
                    STUFF(unsigned, u, u, stride);
                    break;
                case PHType::STRING:
                case PHType::PATH:
                case PHType::LONG:
                case PHType::DOUBLE:
                case PHType::FLOAT:
                case PHType::INTEGER:
                case PHType::SIZE_T:
                    TypeValError(pf->type,list[i].type, yylloc);
                    break;
                case PHType::BAD:
                    ERROR(string("stuffit: logic failure!"), yylloc);
                    exit(1);
                }
                break;
        
            case ValType::STRING:
                switch (pf->type) {
                case PHType::STRING:
                    // This is copying the char * into the std::string
                    *(pf->pval.s) = list[i].value.s;
                    pf->pval.s += pf->stride;
                    pf->nvals++;
                    if(pf->nvals>pf->dim) {
                        ERROR(string("number of values (") + 
                              stringutil::ToString(pf->nvals) + string(") exceeds dimension for ") +
                              string(pf->name) + string("[") + stringutil::ToString(pf->dim) +
                              string("]"), yylloc);
                    }
                    *pf->init = 1;
                    break;
                case PHType::PATH:
                    // This is copying the char * into the fs::path
                    *(pf->pval.p) = list[i].value.s;
                    pf->pval.s += pf->stride;
                    pf->nvals++;
                    if(pf->nvals>pf->dim) {
                        ERROR(string("number of values (") + 
                              stringutil::ToString(pf->nvals) + string(") exceeds dimension for ") +
                              string(pf->name) + string("[") + stringutil::ToString(pf->dim) +
                              string("]"), yylloc);
                    }
                    *pf->init = 1;
                    break;
                case PHType::LONG:
                case PHType::DOUBLE:
                case PHType::FLOAT:
                case PHType::INTEGER:
                case PHType::SIZE_T:
                case PHType::LOGICAL:
                    TypeValError(pf->type,list[i].type, yylloc);
                    break;
                case PHType::BAD:
                    ERROR(string("stuffit: logic failure!"), yylloc);
                    exit(1);
                }
                break;
            }
        }
    }
}
#undef STUFF

// generate a warning for variables which have not recieved a value and for
//   which a warning was desired
void phDriver::checkinit(void) {
    int i;
  
    for (i=0; i<NSYMS; i++) {
        if ((symtab[i].name!=0) && (*symtab[i].init!=1)) {
            fmt::print(stderr,"symtab[{:d}].name = {:s}, init={:d}\n",i,symtab[i].name, *symtab[i].init);
            if(Warnings)
                fmt::print(stderr, "symbol \"{}\" requires a value.\n", symtab[i].name);
        }
    }
}

void phDriver::InitSymtab(void) {
    int i;

    if(Debug) fmt::print(std::cerr, "initializing symtab...");
    init_symtab = 0;
    for(i=0; i<NSYMS; i++) {
        symtab[i].name = (char *) NULL;
        symtab[i].type = PHType::BAD;
        symtab[i].pval.d = (double *) NULL;
        symtab[i].pvalorig.d = (double *) NULL;
        symtab[i].stride = 0;
        symtab[i].init = (int *) NULL;
        symtab[i].dim = 0;
        symtab[i].nvals = 0;
        symtab[i].mapped = (char *) NULL;
    }
}


/**********************************************************************/
/*            External Interface                                      */
/**********************************************************************/

// reset the parser -- remove all entries from the symbol table
//   so that we can start anew...
void phDriver::ResetParser(void) {
    int i;
 
    init_symtab = 1;
    for (i=0; i<NSYMS; i++) {
        if(symtab[i].name != NULL){ free(symtab[i].name); symtab[i].name = NULL; }
        if(symtab[i].init != NULL){ free(symtab[i].init); symtab[i].init = NULL; }
        symtab[i].pval.d = (double *) NULL;
        symtab[i].pvalorig.d = (double *) NULL;
        symtab[i].type = PHType::BAD;                   // don't care
        symtab[i].dim = 0;
        symtab[i].nvals = 0;
    }

    check_flag=0;
    pp_len=0;
    init_symtab=1;

    nzones_prs = 0;
    nzprs_old = 0;
    nwarn = 0;
}

/**********************************************************************/
/* grammer actions                                                    */
/**********************************************************************/

int phDriver::InstallMapVar(char *val, const yy::location& yylloc) {
    int ret = 1;
    SYMENT *psym = lookup(val,1);
    if(psym != (SYMENT *) NULL) {
        DEBUGOUT(string("mapvar for \"") + string(val) + string("\" variables:"), yylloc);
        mapbase = psym;
        for(int i=0; i<pp_len; i++) {
            psym=lookup(id_stack[i],0);
            DEBUGOUT(string("      \"")+string(psym->name) + string("\"."), yylloc);
            if(psym->mapped != NULL) {
                ERROR(string("variable \"") + string(psym->name) + string("\" already mapped to \"") + 
                      string(psym->mapped) + string("\""), yylloc);
                ret = 0;
            }
            else {
                psym->type = mapbase->type;
                psym->pval.d = mapbase->pval.d;
                psym->stride = mapbase->stride;
                psym->mapped = mapbase->name;
                free(psym->init);
                psym->init = mapbase->init;
                psym->dim  = mapbase->dim;
            }
        }
    }
    else {
#if 0
        // this should not be an error...
        ERROR(string("mapvar: variable \"") + string(val) + string("\" not defined."),
              yylloc);
        ret = 0;
#endif
    }
    return ret;
}

int phDriver::InstallVCounter(char *val, const yy::location& yylloc) {
    SYMENT *psym = lookup(val,0);
    DEBUGOUT(string("installing vcounter: \"") + string(psym->name) + string("\"."), yylloc);
    if(psym->pval.i == (int *) NULL) {
        ERROR(string("unknown symbol for vcounter \"") + stringutil::ToString(psym->name) +
              string("\""), yylloc);
        return 0;
    }
    else {
        *psym->pval.i = 0;
        nzones_prs = psym->pval.i;
        nzones_init = psym->init;
        nzprs_old = *nzones_prs;
        zonecounter_psym_old = psym;
    }
    return 1;
}

int phDriver::ProcessEOS(const yy::location& yylloc) {
    if(check_flag == 1) {
        if(nzprs_old !=0) {
            if(nzprs_old != *nzones_prs) {
                ERROR(string("inconsistent number of elements in previous list:\n") + 
                      string("         element counter \"") + 
                      string(zonecounter_psym_old->name) + 
                      string("\", ") + stringutil::ToString(nzprs_old) + string(" elements previously, ") + 
                      stringutil::ToString(*nzones_prs) + string(" found here."), yylloc);
                return 0;
            }
        }
        check_flag = 0;
    }
    return 1;
}

int phDriver::ProcessVec_Start(const yy::location& yylloc) {

    DEBUGOUT(string("vector declaration variables: "), yylloc);
    if(nzones_prs == (int *) NULL) {
        ERROR(string("must specify a vcounter variable before vector statement."), yylloc);
    }
    for(int i=0; i<pp_len; i++) {
        DEBUGOUT(string("       \"") + string(id_stack[i]) + string("\""), yylloc);
        zonespec[i]=lookup(id_stack[i],1);
        if(no_undefined && zonespec[i] == (SYMENT *) NULL) {
            ERROR(string("Unknown variable in assignment: \"") + string(id_stack[i]) + 
                  string("\"."), yylloc);
            return 0;
        }
        zonespec[i]=lookup(id_stack[i],0);
    }
    nzprs_old = *nzones_prs;
    *nzones_prs = 0;
    check_flag = 1;
    return 1;
}

void phDriver::ProcessValue_List(const yy::location& yylloc) {
    (*nzones_prs)++;
    stuffit(zonespec,val_stack,pp_len,1, yylloc);
    pp_len = 0;
    *nzones_init = 1;
}

int phDriver::ProcessEquals(char *val, const yy::location& yylloc) {
    SYMENT *psym = lookup(val,1);
    if(no_undefined && psym == (SYMENT *) NULL) {
        ERROR(string("Unknown variable in assignment: \"") + 
              string("\"."), yylloc);
        return 0;
    }
    else {
        psym = lookup(val,0);
        stuffit(&psym,val_stack,pp_len,0, yylloc);
    }
    return 1;
}

void phDriver::SetID_List(char *val) {
    id_stack[pp_len] = strdup(val);
    pp_len++;
}

void phDriver::ReSetID_List(char *val) {
    pp_len = 0;
    id_stack[pp_len] = strdup(val);
    pp_len++;
}

void phDriver::AddToValStack(LLint c) {
    val_stack[pp_len].value.l = c;
    val_stack[pp_len].type = ValType::INTEGER;
    pp_len++;
}

void phDriver::AddToValStack(double c) {
    val_stack[pp_len].value.d = c;
    val_stack[pp_len].type = ValType::DOUBLE;
    pp_len++;
}

void phDriver::AddToValStack(bool c) {
    if(c)
        val_stack[pp_len].value.u = 1;
    else
        val_stack[pp_len].value.u = 0;
    val_stack[pp_len].type = ValType::LOGICAL;
    pp_len++;
}

void phDriver::AddToValStack(char *c) {
    val_stack[pp_len].value.s = strdup(c);
    val_stack[pp_len].type = ValType::STRING;
    pp_len++;
}


void phDriver::ReSetValStack(LLint c) {
    pp_len = 0;
    val_stack[pp_len].value.l = c;
    val_stack[pp_len].type = ValType::INTEGER;
    pp_len++;
}

void phDriver::ReSetValStack(double c) {
    pp_len = 0;
    val_stack[pp_len].value.d = c;
    val_stack[pp_len].type = ValType::DOUBLE;
    pp_len++;
}

void phDriver::ReSetValStack(bool c) {
    pp_len = 0;
    if(c)
        val_stack[pp_len].value.u = 1;
    else
        val_stack[pp_len].value.u = 0;
    val_stack[pp_len].type = ValType::LOGICAL;
    pp_len++;
}

void phDriver::ReSetValStack(char *c) {
    pp_len = 0;
    val_stack[pp_len].value.s = strdup(c);
    val_stack[pp_len].type = ValType::STRING;
    pp_len++;
}

void phDriver::DumpAllValues(void) {
    for (int i=0; i<NSYMS; i++) {
        if (symtab[i].name==0) break;
        fmt::print(stderr,"symbol {:s} ", symtab[i].name);
        if(symtab[i].dim==1) 
            fmt::print(stderr,"(scalar)\n");
        else
            fmt::print(stderr,"(vector)\n");
        switch(symtab[i].type) {
        case PHType::LOGICAL:
            for(int j=0; j<symtab[i].nvals; j++)
                fmt::print(stderr,"       {:d}\n", symtab[i].pvalorig.u[j]);
            break;
        case PHType::INTEGER:
            for(int j=0; j<symtab[i].nvals; j++)
                fmt::print(stderr,"       {:d}\n", symtab[i].pvalorig.i[j]);
            break;
        case PHType::SIZE_T:
            for(int j=0; j<symtab[i].nvals; j++)
                fmt::print(stderr,"       {:d}\n", symtab[i].pvalorig.z[j]);
            break;
        case PHType::FLOAT:
            for(int j=0; j<symtab[i].nvals; j++)
                fmt::print(stderr,"       {:e}\n", symtab[i].pvalorig.f[j]);
            break;
        case PHType::DOUBLE:
            for(int j=0; j<symtab[i].nvals; j++)
                fmt::print(stderr,"       {:e}\n", symtab[i].pvalorig.d[j]);
            break;
        case PHType::LONG:
            for(int j=0; j<symtab[i].nvals; j++)
                fmt::print(stderr,"       {:d}\n", symtab[i].pvalorig.l[j]);
            break;
        case PHType::STRING:
            for(int j=0; j<symtab[i].nvals; j++)
                fmt::print(stderr,"       {:s}\n", *(symtab[i].pvalorig.s + j*symtab[i].stride));
            break;
        case PHType::PATH:
            for(int j=0; j<symtab[i].nvals; j++)
                fmt::print(stderr,"       {:s}\n", (symtab[i].pvalorig.p + j*symtab[i].stride)->string());
            break;
        case PHType::BAD:
            fmt::print(stderr,"       BAD\n");
            break;
        }
    }
}

