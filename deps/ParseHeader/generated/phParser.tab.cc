// A Bison parser, made by GNU Bison 3.0.4.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.


// First part of user declarations.

#line 37 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:404

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

#include "phParser.tab.hh"

// User implementation prologue.

#line 51 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:412
// Unqualified %code blocks.
#line 37 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:413

# include "detail/phDriver.hh"
using std::string;

#line 58 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:413


#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (/*CONSTCOND*/ false)
# endif


// Suppress unused-variable warnings by "using" E.
#define YYUSE(E) ((void) (E))

// Enable debugging if requested.
#if YYDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << std::endl;                  \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yystack_print_ ();                \
  } while (false)

#else // !YYDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YYUSE(Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void>(0)
# define YY_STACK_PRINT()                static_cast<void>(0)

#endif // !YYDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)


namespace yy {
#line 144 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:479

  /* Return YYSTR after stripping away unnecessary quotes and
     backslashes, so that it's suitable for yyerror.  The heuristic is
     that double-quoting is unnecessary unless the string contains an
     apostrophe, a comma, or backslash (other than backslash-backslash).
     YYSTR is taken from yytname.  */
  std::string
   phParser ::yytnamerr_ (const char *yystr)
  {
    if (*yystr == '"')
      {
        std::string yyr = "";
        char const *yyp = yystr;

        for (;;)
          switch (*++yyp)
            {
            case '\'':
            case ',':
              goto do_not_strip_quotes;

            case '\\':
              if (*++yyp != '\\')
                goto do_not_strip_quotes;
              // Fall through.
            default:
              yyr += *yyp;
              break;

            case '"':
              return yyr;
            }
      do_not_strip_quotes: ;
      }

    return yystr;
  }


  /// Build a parser object.
   phParser :: phParser  (phDriver& driver_yyarg)
    :
#if YYDEBUG
      yydebug_ (false),
      yycdebug_ (&std::cerr),
#endif
      driver (driver_yyarg)
  {}

   phParser ::~ phParser  ()
  {}


  /*---------------.
  | Symbol types.  |
  `---------------*/

  inline
   phParser ::syntax_error::syntax_error (const location_type& l, const std::string& m)
    : std::runtime_error (m)
    , location (l)
  {}

  // basic_symbol.
  template <typename Base>
  inline
   phParser ::basic_symbol<Base>::basic_symbol ()
    : value ()
  {}

  template <typename Base>
  inline
   phParser ::basic_symbol<Base>::basic_symbol (const basic_symbol& other)
    : Base (other)
    , value ()
    , location (other.location)
  {
    value = other.value;
  }


  template <typename Base>
  inline
   phParser ::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, const semantic_type& v, const location_type& l)
    : Base (t)
    , value (v)
    , location (l)
  {}


  /// Constructor for valueless symbols.
  template <typename Base>
  inline
   phParser ::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, const location_type& l)
    : Base (t)
    , value ()
    , location (l)
  {}

  template <typename Base>
  inline
   phParser ::basic_symbol<Base>::~basic_symbol ()
  {
    clear ();
  }

  template <typename Base>
  inline
  void
   phParser ::basic_symbol<Base>::clear ()
  {
    Base::clear ();
  }

  template <typename Base>
  inline
  bool
   phParser ::basic_symbol<Base>::empty () const
  {
    return Base::type_get () == empty_symbol;
  }

  template <typename Base>
  inline
  void
   phParser ::basic_symbol<Base>::move (basic_symbol& s)
  {
    super_type::move(s);
    value = s.value;
    location = s.location;
  }

  // by_type.
  inline
   phParser ::by_type::by_type ()
    : type (empty_symbol)
  {}

  inline
   phParser ::by_type::by_type (const by_type& other)
    : type (other.type)
  {}

  inline
   phParser ::by_type::by_type (token_type t)
    : type (yytranslate_ (t))
  {}

  inline
  void
   phParser ::by_type::clear ()
  {
    type = empty_symbol;
  }

  inline
  void
   phParser ::by_type::move (by_type& that)
  {
    type = that.type;
    that.clear ();
  }

  inline
  int
   phParser ::by_type::type_get () const
  {
    return type;
  }


  // by_state.
  inline
   phParser ::by_state::by_state ()
    : state (empty_state)
  {}

  inline
   phParser ::by_state::by_state (const by_state& other)
    : state (other.state)
  {}

  inline
  void
   phParser ::by_state::clear ()
  {
    state = empty_state;
  }

  inline
  void
   phParser ::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  inline
   phParser ::by_state::by_state (state_type s)
    : state (s)
  {}

  inline
   phParser ::symbol_number_type
   phParser ::by_state::type_get () const
  {
    if (state == empty_state)
      return empty_symbol;
    else
      return yystos_[state];
  }

  inline
   phParser ::stack_symbol_type::stack_symbol_type ()
  {}


  inline
   phParser ::stack_symbol_type::stack_symbol_type (state_type s, symbol_type& that)
    : super_type (s, that.location)
  {
    value = that.value;
    // that is emptied.
    that.type = empty_symbol;
  }

  inline
   phParser ::stack_symbol_type&
   phParser ::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
    value = that.value;
    location = that.location;
    return *this;
  }


  template <typename Base>
  inline
  void
   phParser ::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);

    // User destructor.
    YYUSE (yysym.type_get ());
  }

#if YYDEBUG
  template <typename Base>
  void
   phParser ::yy_print_ (std::ostream& yyo,
                                     const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YYUSE (yyoutput);
    symbol_number_type yytype = yysym.type_get ();
    // Avoid a (spurious) G++ 4.8 warning about "array subscript is
    // below array bounds".
    if (yysym.empty ())
      std::abort ();
    yyo << (yytype < yyntokens_ ? "token" : "nterm")
        << ' ' << yytname_[yytype] << " ("
        << yysym.location << ": ";
    YYUSE (yytype);
    yyo << ')';
  }
#endif

  inline
  void
   phParser ::yypush_ (const char* m, state_type s, symbol_type& sym)
  {
    stack_symbol_type t (s, sym);
    yypush_ (m, t);
  }

  inline
  void
   phParser ::yypush_ (const char* m, stack_symbol_type& s)
  {
    if (m)
      YY_SYMBOL_PRINT (m, s);
    yystack_.push (s);
  }

  inline
  void
   phParser ::yypop_ (unsigned int n)
  {
    yystack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
   phParser ::debug_stream () const
  {
    return *yycdebug_;
  }

  void
   phParser ::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


   phParser ::debug_level_type
   phParser ::debug_level () const
  {
    return yydebug_;
  }

  void
   phParser ::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // YYDEBUG

  inline  phParser ::state_type
   phParser ::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - yyntokens_] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - yyntokens_];
  }

  inline bool
   phParser ::yy_pact_value_is_default_ (int yyvalue)
  {
    return yyvalue == yypact_ninf_;
  }

  inline bool
   phParser ::yy_table_value_is_error_ (int yyvalue)
  {
    return yyvalue == yytable_ninf_;
  }

  int
   phParser ::parse ()
  {
    // State.
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

    // FIXME: This shoud be completely indented.  It is not yet to
    // avoid gratuitous conflicts when merging into the master branch.
    try
      {
    YYCDEBUG << "Starting parse" << std::endl;


    // User initialization code.
    #line 25 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:741
{
    // Initialize the initial location.
    yyla.location.begin.filename = yyla.location.end.filename = &driver.filename;
    driver.yynerrs = &(yynerrs_);
}

#line 524 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:741

    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, yyla);

    // A new symbol was pushed on the stack.
  yynewstate:
    YYCDEBUG << "Entering state " << yystack_[0].state << std::endl;

    // Accept?
    if (yystack_[0].state == yyfinal_)
      goto yyacceptlab;

    goto yybackup;

    // Backup.
  yybackup:

    // Try to take a decision without lookahead.
    yyn = yypact_[yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token: ";
        try
          {
            yyla.type = yytranslate_ (yylex (&yyla.value, &yyla.location, driver));
          }
        catch (const syntax_error& yyexc)
          {
            error (yyexc);
            goto yyerrlab1;
          }
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.type_get ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.type_get ())
      goto yydefault;

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", yyn, yyla);
    goto yynewstate;

  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;

  /*-----------------------------.
  | yyreduce -- Do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_(yystack_[yylen].state, yyr1_[yyn]);
      /* If YYLEN is nonzero, implement the default value of the
         action: '$$ = $1'.  Otherwise, use the top of the stack.

         Otherwise, the following line sets YYLHS.VALUE to garbage.
         This behavior is undocumented and Bison users should not rely
         upon it.  */
      if (yylen)
        yylhs.value = yystack_[yylen - 1].value;
      else
        yylhs.value = yystack_[0].value;

      // Compute the default @$.
      {
        slice<stack_symbol_type, stack_type> slice (yystack_, yylen);
        YYLLOC_DEFAULT (yylhs.location, slice, yylen);
      }

      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
      try
        {
          switch (yyn)
            {
  case 7:
#line 86 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.InstallMapVar((yystack_[2].value.y_str), yystack_[2].location); }
#line 634 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 8:
#line 89 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.InstallVCounter((yystack_[1].value.y_str), yystack_[1].location); }
#line 640 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 9:
#line 92 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    {
                    driver.ERROR(string("syntax error, unexpected string \"") + string((yystack_[2].value.y_str)) + 
                            string("\", expecting '='"), yystack_[2].location);
                  }
#line 649 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 10:
#line 98 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    {
                    driver.ERROR(string("syntax error, unexpected value \"") + 
                                 stringutil::ToString((yystack_[2].value.y_float)) + string("\", expecting 'identifier ='"), yystack_[2].location);
                  }
#line 658 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 11:
#line 104 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    {
                    driver.ERROR(string("syntax error, unexpected value \"") + 
                                 stringutil::ToString((yystack_[2].value.y_lli)) + string("\", expecting 'identifier ='"), yystack_[2].location);
                  }
#line 667 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 12:
#line 110 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ProcessEOS(yystack_[0].location); }
#line 673 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 13:
#line 114 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ProcessVec_Start(yystack_[2].location); }
#line 679 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 16:
#line 122 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ProcessValue_List(yyla.location); }
#line 685 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 17:
#line 127 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ProcessEquals((yystack_[3].value.y_str),yystack_[3].location); }
#line 691 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 18:
#line 132 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.SetID_List((yystack_[0].value.y_str)); }
#line 697 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 19:
#line 134 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ReSetID_List((yystack_[0].value.y_str)); }
#line 703 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 20:
#line 136 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ERROR(string("numeric value in id list"), yystack_[0].location); }
#line 709 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 21:
#line 138 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ERROR(string("numeric value in id list"), yystack_[0].location); }
#line 715 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 22:
#line 140 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ERROR(string("logical value in id list"), yystack_[0].location); }
#line 721 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 23:
#line 142 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ERROR(string("logical value in id list"), yystack_[0].location); }
#line 727 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 24:
#line 147 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.AddToValStack((yystack_[0].value.y_lli)); }
#line 733 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 25:
#line 149 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.AddToValStack((yystack_[0].value.y_float)); }
#line 739 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 26:
#line 151 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.AddToValStack(true); }
#line 745 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 27:
#line 153 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.AddToValStack(false); }
#line 751 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 28:
#line 155 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.AddToValStack((yystack_[0].value.y_str)); }
#line 757 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 29:
#line 157 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ReSetValStack((yystack_[0].value.y_lli)); }
#line 763 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 30:
#line 159 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ReSetValStack((yystack_[0].value.y_float)); }
#line 769 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 31:
#line 161 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ReSetValStack(true); }
#line 775 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 32:
#line 163 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ReSetValStack(false); }
#line 781 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;

  case 33:
#line 165 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:859
    { driver.ReSetValStack((yystack_[0].value.y_str)); }
#line 787 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
    break;


#line 791 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:859
            default:
              break;
            }
        }
      catch (const syntax_error& yyexc)
        {
          error (yyexc);
          YYERROR;
        }
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;
      YY_STACK_PRINT ();

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, yylhs);
    }
    goto yynewstate;

  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        error (yyla.location, yysyntax_error_ (yystack_[0].state, yyla));
      }


    yyerror_range[1].location = yyla.location;
    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.type_get () == yyeof_)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:

    /* Pacify compilers like GCC when the user code never invokes
       YYERROR and the label yyerrorlab therefore never appears in user
       code.  */
    if (false)
      goto yyerrorlab;
    yyerror_range[1].location = yystack_[yylen - 1].location;
    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    goto yyerrlab1;

  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    {
      stack_symbol_type error_token;
      for (;;)
        {
          yyn = yypact_[yystack_[0].state];
          if (!yy_pact_value_is_default_ (yyn))
            {
              yyn += yyterror_;
              if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == yyterror_)
                {
                  yyn = yytable_[yyn];
                  if (0 < yyn)
                    break;
                }
            }

          // Pop the current state because it cannot handle the error token.
          if (yystack_.size () == 1)
            YYABORT;

          yyerror_range[1].location = yystack_[0].location;
          yy_destroy_ ("Error: popping", yystack_[0]);
          yypop_ ();
          YY_STACK_PRINT ();
        }

      yyerror_range[2].location = yyla.location;
      YYLLOC_DEFAULT (error_token.location, yyerror_range, 2);

      // Shift the error token.
      error_token.state = yyn;
      yypush_ ("Shifting", error_token);
    }
    goto yynewstate;

    // Accept.
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;

    // Abort.
  yyabortlab:
    yyresult = 1;
    goto yyreturn;

  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack"
                 << std::endl;
        // Do not try to display the values of the reclaimed symbols,
        // as their printer might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
  }

  void
   phParser ::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what());
  }

  // Generate an error message.
  std::string
   phParser ::yysyntax_error_ (state_type yystate, const symbol_type& yyla) const
  {
    // Number of reported tokens (one for the "unexpected", one per
    // "expected").
    size_t yycount = 0;
    // Its maximum.
    enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
    // Arguments of yyformat.
    char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];

    /* There are many possibilities here to consider:
       - If this state is a consistent state with a default action, then
         the only way this function was invoked is if the default action
         is an error action.  In that case, don't check for expected
         tokens because there are none.
       - The only way there can be no lookahead present (in yyla) is
         if this state is a consistent state with a default action.
         Thus, detecting the absence of a lookahead is sufficient to
         determine that there is no unexpected or expected token to
         report.  In that case, just report a simple "syntax error".
       - Don't assume there isn't a lookahead just because this state is
         a consistent state with a default action.  There might have
         been a previous inconsistent state, consistent state with a
         non-default action, or user semantic action that manipulated
         yyla.  (However, yyla is currently not documented for users.)
       - Of course, the expected token list depends on states to have
         correct lookahead information, and it depends on the parser not
         to perform extra reductions after fetching a lookahead from the
         scanner and before detecting a syntax error.  Thus, state
         merging (from LALR or IELR) and default reductions corrupt the
         expected token list.  However, the list is correct for
         canonical LR with one exception: it will still contain any
         token that will not be accepted due to an error action in a
         later state.
    */
    if (!yyla.empty ())
      {
        int yytoken = yyla.type_get ();
        yyarg[yycount++] = yytname_[yytoken];
        int yyn = yypact_[yystate];
        if (!yy_pact_value_is_default_ (yyn))
          {
            /* Start YYX at -YYN if negative to avoid negative indexes in
               YYCHECK.  In other words, skip the first -YYN actions for
               this state because they are default actions.  */
            int yyxbegin = yyn < 0 ? -yyn : 0;
            // Stay within bounds of both yycheck and yytname.
            int yychecklim = yylast_ - yyn + 1;
            int yyxend = yychecklim < yyntokens_ ? yychecklim : yyntokens_;
            for (int yyx = yyxbegin; yyx < yyxend; ++yyx)
              if (yycheck_[yyx + yyn] == yyx && yyx != yyterror_
                  && !yy_table_value_is_error_ (yytable_[yyx + yyn]))
                {
                  if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                    {
                      yycount = 1;
                      break;
                    }
                  else
                    yyarg[yycount++] = yytname_[yyx];
                }
          }
      }

    char const* yyformat = YY_NULLPTR;
    switch (yycount)
      {
#define YYCASE_(N, S)                         \
        case N:                               \
          yyformat = S;                       \
        break
        YYCASE_(0, YY_("syntax error"));
        YYCASE_(1, YY_("syntax error, unexpected %s"));
        YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
        YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
        YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
        YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
#undef YYCASE_
      }

    std::string yyres;
    // Argument number.
    size_t yyi = 0;
    for (char const* yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount)
        {
          yyres += yytnamerr_ (yyarg[yyi++]);
          ++yyp;
        }
      else
        yyres += *yyp;
    return yyres;
  }


  const signed char  phParser ::yypact_ninf_ = -6;

  const signed char  phParser ::yytable_ninf_ = -1;

  const signed char
   phParser ::yypact_[] =
  {
      15,     1,    -6,    81,    23,    81,    -1,     6,     9,     3,
      -6,    81,    -6,    -6,    -6,    -6,    -6,    -6,    -6,    33,
      81,    40,    47,     6,    -6,    54,     4,    -6,    -6,    81,
      -6,    61,    -6,    -6,    -6,    -6,    -6,    -6,    68,    -6,
      -6,    75,    -6,    -6,    -6,    -6,    -6,    -6,    -6,    -6,
      -6,    -6,    -6
  };

  const unsigned char
   phParser ::yydefact_[] =
  {
       0,     0,    12,     0,     0,     0,     0,     0,     0,     0,
       3,     0,     5,     4,    32,    30,    33,    29,    31,     0,
       0,     0,     0,     0,    19,     0,     0,     1,     2,     6,
      15,     0,    10,    27,    25,    28,    24,    26,     0,     9,
      11,     0,    13,    23,    21,    18,    20,    22,     8,    14,
      16,    17,     7
  };

  const signed char
   phParser ::yypgoto_[] =
  {
      -6,    -6,    10,    -6,    -6,    -5,    -6,     8,    -3
  };

  const signed char
   phParser ::yydefgoto_[] =
  {
      -1,     9,    10,    11,    29,    30,    12,    25,    31
  };

  const unsigned char
   phParser ::yytable_[] =
  {
      19,    21,    22,    27,    13,    23,     2,    48,     3,     4,
       5,     6,    24,     7,     8,    26,     1,    38,     2,    28,
       3,     4,     5,     6,    49,     7,     8,    14,    15,    16,
      17,    41,    18,     0,     0,    20,    32,    33,    34,    35,
      36,     0,    37,    39,    33,    34,    35,    36,     0,    37,
      40,    33,    34,    35,    36,     0,    37,    42,    43,    44,
      45,    46,     0,    47,    50,    33,    34,    35,    36,     0,
      37,    51,    33,    34,    35,    36,     0,    37,    52,    43,
      44,    45,    46,     0,    47,    14,    15,    16,    17,     0,
      18
  };

  const signed char
   phParser ::yycheck_[] =
  {
       3,     4,     5,     0,     3,     6,     3,     3,     5,     6,
       7,     8,     6,    10,    11,     6,     1,    20,     3,     9,
       5,     6,     7,     8,    29,    10,    11,     4,     5,     6,
       7,    23,     9,    -1,    -1,    12,     3,     4,     5,     6,
       7,    -1,     9,     3,     4,     5,     6,     7,    -1,     9,
       3,     4,     5,     6,     7,    -1,     9,     3,     4,     5,
       6,     7,    -1,     9,     3,     4,     5,     6,     7,    -1,
       9,     3,     4,     5,     6,     7,    -1,     9,     3,     4,
       5,     6,     7,    -1,     9,     4,     5,     6,     7,    -1,
       9
  };

  const unsigned char
   phParser ::yystos_[] =
  {
       0,     1,     3,     5,     6,     7,     8,    10,    11,    14,
      15,    16,    19,     3,     4,     5,     6,     7,     9,    21,
      12,    21,    21,     6,     6,    20,     6,     0,    15,    17,
      18,    21,     3,     4,     5,     6,     7,     9,    21,     3,
       3,    20,     3,     4,     5,     6,     7,     9,     3,    18,
       3,     3,     3
  };

  const unsigned char
   phParser ::yyr1_[] =
  {
       0,    13,    14,    14,    14,    15,    15,    15,    15,    15,
      15,    15,    15,    16,    17,    17,    18,    19,    20,    20,
      20,    20,    20,    20,    21,    21,    21,    21,    21,    21,
      21,    21,    21,    21
  };

  const unsigned char
   phParser ::yyr2_[] =
  {
       0,     2,     2,     1,     2,     1,     2,     4,     3,     3,
       3,     3,     1,     3,     2,     1,     2,     4,     2,     1,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     1,
       1,     1,     1,     1
  };



  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a yyntokens_, nonterminals.
  const char*
  const  phParser ::yytname_[] =
  {
  "\"end of file\"", "error", "$undefined", "EOS", "FALSEval",
  "\"floating-point value\"", "\"string\"", "\"integer value\"",
  "MAPVARkey", "TRUEval", "VECTOR", "VCOUNTER", "'='", "$accept",
  "gen_stmt_list", "gen_stmt", "vec_start", "vec_list", "vec_stmt",
  "id_stmt", "id_list", "value_list", YY_NULLPTR
  };

#if YYDEBUG
  const unsigned char
   phParser ::yyrline_[] =
  {
       0,    75,    75,    76,    77,    81,    83,    85,    88,    91,
      97,   103,   109,   113,   117,   118,   121,   126,   131,   133,
     135,   137,   139,   141,   146,   148,   150,   152,   154,   156,
     158,   160,   162,   164
  };

  // Print the state stack on the debug stream.
  void
   phParser ::yystack_print_ ()
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator
           i = yystack_.begin (),
           i_end = yystack_.end ();
         i != i_end; ++i)
      *yycdebug_ << ' ' << i->state;
    *yycdebug_ << std::endl;
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void
   phParser ::yy_reduce_print_ (int yyrule)
  {
    unsigned int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):" << std::endl;
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // YYDEBUG

  // Symbol number corresponding to token number t.
  inline
   phParser ::token_number_type
   phParser ::yytranslate_ (int t)
  {
    static
    const token_number_type
    translate_table[] =
    {
     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    12,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11
    };
    const unsigned int user_token_number_max_ = 266;
    const token_number_type undef_token_ = 2;

    if (static_cast<int>(t) <= yyeof_)
      return yyeof_;
    else if (static_cast<unsigned int> (t) <= user_token_number_max_)
      return translate_table[t];
    else
      return undef_token_;
  }


} // yy
#line 1246 "subprojects/ParseHeader/phParser.tab.cc" // lalr1.cc:1167
#line 168 "../subprojects/ParseHeader/src/phParser.yy" // lalr1.cc:1168


// the error member function registers the errors to the driver. 
void yy::phParser::error (const yy::phParser::location_type& l, const std::string& m) {
    (*driver.yynerrs)--; // ERROR increments yynerrs so don't do twice
    driver.ERROR(m,l);
}

