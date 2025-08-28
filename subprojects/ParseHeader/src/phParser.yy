/* -*- C++ -*- */
%skeleton "lalr1.cc"
%require "2.4"
%defines
%define parser_class_name { phParser }
%expect 3

// declarations/inclusions needed to define the %union. 
// forward declaration of driver
%code requires {
# include <string>
# include "detail/stringutil.hh"
class phDriver;
}

// driver is passed by reference to the parser and to the scanner
%parse-param { phDriver& driver }
%lex-param   { phDriver& driver }

// request the location tracking feature, and initialize the first location's file name.
// Afterward new locations are computed relatively to the previous locations: the file 
//    name will be automatically propagated.
%locations
%initial-action
{
    // Initialize the initial location.
    @$.begin.filename = @$.end.filename = &driver.filename;
    driver.yynerrs = &(yynerrs_);
};

// enable parser tracing and verbose error messages
%debug
%error-verbose

// code between ‘%code {’ and ‘}’ is output in the *.cc file; 
//   it needs detailed knowledge about the driver
%code {
# include "detail/phDriver.hh"
using std::string;
}

// possible types for semantic values
%union	{
    struct symtab *y_sim;	// indentifier
    char *y_str;		// string
    long long int y_lli;	// int; all integers are long long ints until stuffed
    double y_float;	        // float; all floats are doubles until stuffed
}

// The token numbered as 0 corresponds to end of file; 
// the following line allows for nicer error messages referring to “end of file” instead of “$end”.
// Similarly user friendly named are provided for each symbol.
// Note that the tokens names are prefixed by TOKEN_ to avoid name clashes. 

%token        END      0 "end of file"
%token EOS
%token FALSEval
%token <y_float> FLOATval "floating-point value"
%token <y_str> IDval      "string"
%token <y_lli> INTval     "integer value"
%token MAPVARkey
%token TRUEval
%token VECTOR
%token VCOUNTER

 // To enable memory deallocation during error recovery, use %destructor. 
 //%printer    { debug_stream () << *$$; } IDval
 //%destructor { delete $$; } IDval

%%

 // the grammer describing a general input file format,

gen_stmt_list
		: gen_stmt_list gen_stmt
		| gen_stmt
                | error EOS          // allows parsing of rest of file to find more errors...
		;

gen_stmt
                : id_stmt

                | vec_start vec_list

		| MAPVARkey IDval id_list EOS
                  { driver.InstallMapVar($2, @2); }

                | VCOUNTER IDval EOS
                  { driver.InstallVCounter($2, @2); }

                | IDval value_list EOS
                  {
                    driver.ERROR(string("syntax error, unexpected string \"") + string($1) + 
                            string("\", expecting '='"), @1);
                  }

                | FLOATval value_list EOS
                  {
                    driver.ERROR(string("syntax error, unexpected value \"") + 
                                 stringutil::ToString($1) + string("\", expecting 'identifier ='"), @1);
                  }

                | INTval value_list EOS
                  {
                    driver.ERROR(string("syntax error, unexpected value \"") + 
                                 stringutil::ToString($1) + string("\", expecting 'identifier ='"), @1);
                  }

		| EOS
                  { driver.ProcessEOS(@1); }
		;

vec_start      : VECTOR id_list EOS
                  { driver.ProcessVec_Start(@1); }
                ;

vec_list       : vec_list vec_stmt
                | vec_stmt
                ;

vec_stmt       : value_list EOS
                  { driver.ProcessValue_List(yyla.location); }
                ;

id_stmt
		: IDval '=' value_list EOS
                  { driver.ProcessEquals($1,@1); }
                ;

id_list
		: id_list IDval
		  { driver.SetID_List($2); }
		| IDval
		  { driver.ReSetID_List($1); }
                | id_list INTval
                  { driver.ERROR(string("numeric value in id list"), @2); }
                | id_list FLOATval
                  { driver.ERROR(string("numeric value in id list"), @2); }
                | id_list TRUEval
                  { driver.ERROR(string("logical value in id list"), @2); }
                | id_list FALSEval
                  { driver.ERROR(string("logical value in id list"), @2); }
		;

value_list
                : value_list INTval
                  { driver.AddToValStack($2); }
		| value_list FLOATval
                  { driver.AddToValStack($2); }
		| value_list TRUEval
                  { driver.AddToValStack(true); }
		| value_list FALSEval
                  { driver.AddToValStack(false); }
		| value_list IDval
                  { driver.AddToValStack($2); }
		| INTval 
                  { driver.ReSetValStack($1); }
		| FLOATval
                  { driver.ReSetValStack($1); }
		| TRUEval
                  { driver.ReSetValStack(true); }
		| FALSEval
                  { driver.ReSetValStack(false); }
		| IDval
                  { driver.ReSetValStack($1); }
		;

%%

// the error member function registers the errors to the driver. 
void yy::phParser::error (const yy::phParser::location_type& l, const std::string& m) {
    (*driver.yynerrs)--; // ERROR increments yynerrs so don't do twice
    driver.ERROR(m,l);
}

