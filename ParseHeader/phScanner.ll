/* -*- C++ -*- */
%option outfile="phScanner.cc"
%option header-file="phScanner.hh"
%{
# include <cstdlib>
# include <cerrno>
# include <climits>
# include <string>
# include <assert.h>
# include "stringutil.hh"
# include "phDriver.hh"
# include "phParser.tab.hh"
     
    /* Work around an incompatibility in flex (at least versions
       2.5.31 through 2.5.33): it generates code that does
       not conform to C89.  See Debian bug 333231
       <http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=333231>.  */
# undef yywrap
# define yywrap() 1
    
    /* By default yylex returns int, we use token_type.
       Unfortunately yyterminate by default returns 0, which is
       not of token_type.  */
#define yyterminate() return token::END

/*
 The rules are simple, just note the use of the driver to report errors.
 It is convenient to use a typedef to shorten yy::phParser::token::identifier
 into token::identifier for instance. 
*/

    typedef yy::phParser::token token;
    typedef yy::phParser::token_type token_type;


//PAP
  // Stacks for included files and filenames 
#define MAX_INCLUDE_DEPTH 100
  YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
  int include_stack_ptr = 0;
  char *current_input_name;
  char *file_name_stack[MAX_INCLUDE_DEPTH];
  int file_line_num[MAX_INCLUDE_DEPTH];

  double myatod(char *a);
    char *remove_quotes(char *text, int len);

  // reserved word scanner used to avoid searching for all keywords via lex
    token_type screen(char *, int);

  // reserved keyword table
  static struct key_table {
      char const * key_name;                       // keyword name
      token_type key_yylex;                         // value returned by yylex
  }  key_table [] = {                      // NB: in sorted order
      {"eos",          token::EOS},
      {"false",        token::FALSEval},
      {"float",        token::FLOATval},
      {"id",           token::IDval},
      {"int",          token::INTval},
      {"mapvar",       token::MAPVARkey},
      {"true",         token::TRUEval},
      {"vector",       token::VECTOR},
      {"vcounter",     token::VCOUNTER}
  };

    // stack for bison locations
#define MAX_INCLUDE_STACK 1024
    yy::location location_stack[MAX_INCLUDE_STACK];
    std::string file_name_string;
    int include_level = 0;

%}

/*
 Because there is no #include-like feature we don't need yywrap.
 We don't need unput either, and we parse an actual file.
 This is not an interactive session with the user.
 Finally we enable the scanner tracing features.
*/
%option noyywrap nounput batch debug

 /*
  The following paragraph suffices to track locations accurately.
  Each time yylex is invoked, the begin position is moved onto the end position.
  Then when a pattern is matched, the end position is advanced of its width.
  In case it matched ends of lines, the end cursor is adjusted, 
  and each time blanks are matched, the begin cursor is moved onto the end cursor
  to effectively ignore the blanks preceding tokens. 
  Comments would be treated equally.
 */
%{
# define YY_USER_ACTION  yylloc->columns(yyleng);

%}

/* PAP */
/* the "incl" state is used for picking up the name of an include file */
%x incl
 /* the "IN_COMMENT" is used for multi-line comments */
%s IN_COMMENT

/* comments follow semicolons and extend to the end of the line */
comment        #[^\n]*

/* multi-line comments are delimited by this token on a line by itself */
mlcomment      "##"

/* blanks and tabs are "whitespace" which separate tokens */
white_space    [ \t]+

/* a backslash, followed by any amount of whitespace, at the end of a line
   signifies continuation onto the next line */
continuation   \\[ \t]*\n

/* the end of a statement if not in a continuation */
eos            [\n]

/* a variable name or unquoted-string value; 
     anything beginning with a letter or one of (_.$#&)
     followed by any of the same and numbers. */
id             [a-zA-Z_.$][a-zA-Z_.$0-9]*

/* this complex mess matches all forms of floating point numbers known to C and Fortran:
     (+ or - or nothing)
   followed by one of:
     1.234
     .1234
     12.34
   followed by one of:
     (D or E) 1234
     (+ or -) 1234  ( To match some Fortran output which omits the D or E )
     (D or E) (+ or -) 1234
     nothing
   OR (for the lazy; i.e. "1e21"
     (+ or - or nothing)
   followed by
     1213
   followed by 
     (D or E) 1234
     (D or E) (+ or -) 1234
   Clearly, 1+24 should be illegal!
*/

sign           ("+"|"-")
mant1          (([0-9]+"."[0-9]*)|([0-9]*"."[0-9]+))
exp1           ((((D|d|E|e)?("+"|"-"))|((D|e|E|e)("+"|"-")?))[0-9]+)
mant2          ([0-9]+)
exp2           (D|d|E|e)("+"|"-")?([0-9]+)
float          ({sign}?{mant1}{exp1}?)|({sign}?{mant2}{exp2})

/* matches integers; an optional sign followed by any
     string of only numbers. No embedded whitespace is allowed. */
int            {sign}?[0-9]+

/* quoted strings; may be enclosed in either matching single-
     or double quotes. */
quoted_string  (\"[^\n\"]*\"|\'[^\n\']*\')

/* matches anything else */
other          [a-zA-Z_+=0-9.,><?/\|{}!@#$%^&*()]*


/* PAP */



%%

%{
    yylloc->step ();
%}

<INITIAL>{
^{mlcomment} BEGIN(IN_COMMENT);
}
<IN_COMMENT>{
^{mlcomment} BEGIN(INITIAL);
[^#\n]+                   { yylloc->step (); }
"#"
\n                        { yylloc->lines (yyleng); yylloc->step (); }
}

{comment}               {  /* do nothing for comments; this must come after ml comment */
                           yylloc->step ();
                        }

include                 BEGIN(incl);

<incl>{quoted_string}   { /* got the include file name */
                            char *tmp = remove_quotes(yytext, yyleng);
                            yyin = fopen( tmp, "r" );
                            if ( yyin == (FILE*) NULL ) {
                                driver.MESSAGE(std::string("failed to open include file \"") + std::string(tmp) + 
                                               std::string("\". exiting..."), *yylloc);
                                exit(1);
                            }
                            else {
                                driver.DEBUGOUT(std::string("including file: \"") + std::string(tmp) +
                                                            std::string("\"."), *yylloc);
                            }

                            yypush_buffer_state(yy_create_buffer( yyin, YY_BUF_SIZE ));

                            location_stack[include_level++] = *yylloc;
                            if(include_level>=MAX_INCLUDE_STACK) {
                                driver.MESSAGE(std::string("ERROR: exceeded maximum include stack depth: ") +
                                               ToString(MAX_INCLUDE_STACK) + std::string("."), *yylloc);
                                exit(1);
                            }
                            file_name_string.assign(tmp);
                            yylloc->initialize(&file_name_string);
                            BEGIN(INITIAL);
                        }

<<EOF>>                 { /* At end of file, return to previous file if any */
                            yypop_buffer_state();
                            if ( !YY_CURRENT_BUFFER ) yyterminate();
                            *yylloc = location_stack[--include_level];
                        }

{id}                    { /* look up an id in the keyword list. */
                            yylval->y_str = strdup(yytext); return screen(yytext, yyleng);
                        }

{quoted_string}         { /* get length of string and add trailing blank */
                            yylval->y_str = remove_quotes(yytext, yyleng);
                            return token::IDval;
                        }

{float}                 { /* return a float as a value; use our own atod function
                             which groks the D's fortran uses for double precision
                             values. */
                          yylval->y_float = myatod(yytext); return token::FLOATval;
                        }

{int}                   { /* the system atol is ok for integers */
                          yylval->y_lli = atoll(yytext); return token::INTval;
                        }

{continuation}          { /* keep on going for continuations */ 
                          yylloc->lines(1); yylloc->step();
                        }

{eos}                   { 
                          yylloc->lines(yyleng); yylloc->step();
                          return token::EOS; /* return an end of statement token */ 
                        }

{white_space}           { /* the meaning of whitespace is only to separate tokens */ 
                          yylloc->step();
                        }

{other}                 { return token_type (yytext[0]); /* return the value */ }


%%

#include <cmath>

// Match the token to the keyword; all keywords are case-insensitive.
// non-keywords are returned as IDval's
token_type screen(char* intext, int len) {
    int len_key_table, i;
    char *text;
    text = strdup(intext);
    // text = strndup(intext, len);
    strnlwr(text, len);
    len_key_table = (int) (sizeof(key_table)/sizeof(key_table[0]));
    
    for (i = 0; i < len_key_table; i++) {
        if ((strcmp(key_table[i].key_name, text)) == 0) {
            free(text);
            return(key_table[i].key_yylex);
        }
    }
    free(text);
    return token::IDval;
}

// string to double conversion (better than the system's,
//   as it handles some weird Fortran constructions -- see lex file for details)
double myatod(char *s) {
    double val, power, eval;
    int i, sign, esign;
  
    for(i=0;isspace(s[i]); i++) ;               // skip over white space 
  
    sign = (s[i] == '-') ? -1 : 1;
    if( (s[i]=='-') | (s[i]=='+') ) i++;        // skip over sign 
    for(val = 0.0; isdigit(s[i]); i++)
        val = 10.0 * val + (s[i]-'0');
    if(s[i] == '.') i++;                        // skip over decimal point 
    for(power=1.0; isdigit(s[i]); i++) {
        val = 10.0 * val + (s[i]-'0');
        power *= 10.0;
    }
    if( (s[i] == 'e') | (s[i] == 'E') | (s[i] == 'd') | (s[i] == 'D') ) {
        // there is an exponent 
        i++;
        esign = (s[i] == '-') ? -1 : 1;
        // skip over sign 
        if((s[i]=='-') | (s[i]=='+') ) i++;    
        for(eval = 0.0; isdigit(s[i]); i++)
            eval = 10.0 * eval + (s[i]-'0');
    }
    else
    { esign = 1; eval = 0.0;}
  
    return (sign*val/power*pow(10.0,esign*eval));
}


char *remove_quotes(char *text, int len) {
    char *ctext = new char[len+1];
    char c = text[0];
    assert(c=='"' || c=='\'');
    int i;
    for(i=1; i<len; i++) {
        if(text[i]==c) break;
        ctext[i-1] = text[i];
    }
    ctext[i-1] = 0x0;
    return ctext;
}
