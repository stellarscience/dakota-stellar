/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with nidr or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum nidrtokentype {
     KEY_WORD = 258,
     IDENTIFIER = 259,
     REAL = 260,
     STRING = 261,
     QUOTED_STRING = 262,
     SEPARATOR = 263,
     EQUALS = 264,
     KEYWORDERROR = 265,
     END = 266,
     EXIT = 267
   };
#endif
/* Tokens.  */
#define KEY_WORD 258
#define IDENTIFIER 259
#define REAL 260
#define STRING 261
#define QUOTED_STRING 262
#define SEPARATOR 263
#define EQUALS 264
#define KEYWORDERROR 265
#define END 266
#define EXIT 267




/* Copy the first part of user declarations.  */
/* #line 34 "nidrparse.y" */


#define const

#ifndef NIDR_H
#include "nidr.h"
#endif
#include <string.h>

#include <stdlib.h>
#include <string.h>
extern int nidrlex(void);
#include <math.h>
#include <stdio.h>

extern int		nidrLineNumber;
extern void		(*nidr_bufr)(Real, int);
extern void		(*nidr_bufs)(const char*);
extern int		nidr_cleanup(void);
extern const char*	nidr_keyword_name(void);
extern void		nidr_reinit(void);
extern void		nidr_reset(void);
extern void		nidr_setup(const char*, FILE*);
extern void		nidr_signal_parse_error(void);
extern void		reset_nidrlex_state(void);
extern char		nidr_please_refer[];



#if  defined __STDC__  ||  defined __cplusplus
void	nidrerror( char* s );
#else
void	nidrerror();
#endif

#ifdef IDR_DEBUG
#define	PRINT_DEBUG_NL	putc('\n', stderr);
#else
#define	PRINT_DEBUG_NL
#endif

 static void
please_refer(void)
{
	fprintf(stderr, "%s\n", nidr_please_refer);
	PRINT_DEBUG_NL;
	nidr_signal_parse_error();
	reset_nidrlex_state();
	}

 static int str_lineno;
 static char garbled_fmt[] =
	"Unrecognized data.\n\tIdentifier value garbled for identifier '%s'.\n";


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
/* #line 91 "nidrparse.y" */
{
	KeyWord *keyword;
	KeyWord	*identifier;
	Real	real;
	int	integer;
	char*	string;
	char*	qstring;
}
/* Line 193 of yacc.c.  */
/* #line 182 "nidr-parser.c" */
	YYSTYPE;
# define nidrstype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
/* #line 195 "nidr-parser.c" */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 nidrtype_uint8;
#else
typedef unsigned char nidrtype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 nidrtype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char nidrtype_int8;
#else
typedef short int nidrtype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 nidrtype_uint16;
#else
typedef unsigned short int nidrtype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 nidrtype_int16;
#else
typedef short int nidrtype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined nidroverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined nidroverflow || YYERROR_VERBOSE */


#if (! defined nidroverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union nidralloc
{
  nidrtype_int16 nidrss;
  YYSTYPE nidrvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union nidralloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (nidrtype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T nidri;				\
	  for (nidri = 0; nidri < (Count); nidri++)	\
	    (To)[nidri] = (From)[nidri];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T nidrnewbytes;						\
	YYCOPY (&nidrptr->Stack, Stack, nidrsize);				\
	Stack = &nidrptr->Stack;						\
	nidrnewbytes = nidrstacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	nidrptr += nidrnewbytes / sizeof (*nidrptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   33

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  15
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  8
/* YYNRULES -- Number of rules.  */
#define YYNRULES  23
/* YYNRULES -- Number of states.  */
#define YYNSTATES  33

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   267

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? nidrtranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const nidrtype_uint8 nidrtranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,    13,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    14,     2,
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
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const nidrtype_uint8 nidrprhs[] =
{
       0,     0,     3,     4,     7,     9,    12,    13,    19,    20,
      25,    28,    30,    32,    35,    39,    41,    44,    46,    48,
      51,    54,    56,    58
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const nidrtype_int8 nidrrhs[] =
{
      16,     0,    -1,    -1,    16,    17,    -1,    11,    -1,    20,
      11,    -1,    -1,    20,     6,    18,     1,    11,    -1,    -1,
       6,    19,     1,    11,    -1,     1,    11,    -1,    10,    -1,
       3,    -1,     3,    21,    -1,     3,     9,    21,    -1,    22,
      -1,    21,    22,    -1,     8,    -1,     5,    -1,    13,     5,
      -1,    14,     5,    -1,     7,    -1,     4,    -1,     4,     9,
      -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const nidrtype_uint8 nidrrline[] =
{
       0,   123,   123,   124,   129,   130,   131,   131,   138,   138,
     144,   153,   161,   162,   163,   167,   168,   172,   173,   174,
     175,   176,   177,   178
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const nidrtname[] =
{
  "$end", "error", "$undefined", "KEY_WORD", "IDENTIFIER", "REAL",
  "STRING", "QUOTED_STRING", "SEPARATOR", "EQUALS", "KEYWORDERROR", "END",
  "EXIT", "'*'", "':'", "$accept", "statements", "statement", "@1", "@2",
  "command", "data", "datum", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const nidrtype_uint16 nidrtoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,    42,    58
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const nidrtype_uint8 nidrr1[] =
{
       0,    15,    16,    16,    17,    17,    18,    17,    19,    17,
      17,    17,    20,    20,    20,    21,    21,    22,    22,    22,
      22,    22,    22,    22
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const nidrtype_uint8 nidrr2[] =
{
       0,     2,     0,     2,     1,     2,     0,     5,     0,     4,
       2,     1,     1,     2,     3,     1,     2,     1,     1,     2,
       2,     1,     1,     2
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const nidrtype_uint8 nidrdefact[] =
{
       2,     0,     1,     0,    12,     8,    11,     4,     3,     0,
      10,    22,    18,    21,    17,     0,     0,     0,    13,    15,
       0,     6,     5,    23,    14,    19,    20,    16,     0,     0,
       9,     0,     7
};

/* YYDEFGOTO[NTERM-NUM].  */
static const nidrtype_int8 nidrdefgoto[] =
{
      -1,     1,     8,    29,    20,     9,    18,    19
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -17
static const nidrtype_int8 nidrpact[] =
{
     -17,     0,   -17,    -6,     8,   -17,   -17,   -17,   -17,    -2,
     -17,     5,   -17,   -17,   -17,    19,     2,    13,    19,   -17,
      18,   -17,   -17,   -17,    19,   -17,   -17,   -17,     9,    24,
     -17,    17,   -17
};

/* YYPGOTO[NTERM-NUM].  */
static const nidrtype_int8 nidrpgoto[] =
{
     -17,   -17,   -17,   -17,   -17,   -17,    14,   -16
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const nidrtype_uint8 nidrtable[] =
{
       2,     3,    27,     4,    21,    10,     5,    25,    27,    22,
       6,     7,    11,    12,    23,    13,    14,    15,    26,    28,
      30,    16,    17,    11,    12,    31,    13,    14,    32,    24,
       0,     0,    16,    17
};

static const nidrtype_int8 nidrcheck[] =
{
       0,     1,    18,     3,     6,    11,     6,     5,    24,    11,
      10,    11,     4,     5,     9,     7,     8,     9,     5,     1,
      11,    13,    14,     4,     5,     1,     7,     8,    11,    15,
      -1,    -1,    13,    14
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const nidrtype_uint8 nidrstos[] =
{
       0,    16,     0,     1,     3,     6,    10,    11,    17,    20,
      11,     4,     5,     7,     8,     9,    13,    14,    21,    22,
      19,     6,    11,     9,    21,     5,     5,    22,     1,    18,
      11,     1,    11
};

#define nidrerrok		(nidrerrstatus = 0)
#define nidrclearin	(nidrchar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto nidracceptlab
#define YYABORT		goto nidrabortlab
#define YYERROR		goto nidrerrorlab


/* Like YYERROR except do call nidrerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto nidrerrlab

#define YYRECOVERING()  (!!nidrerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (nidrchar == YYEMPTY && nidrlen == 1)				\
    {								\
      nidrchar = (Token);						\
      nidrlval = (Value);						\
      nidrtoken = YYTRANSLATE (nidrchar);				\
      YYPOPSTACK (1);						\
      goto nidrbackup;						\
    }								\
  else								\
    {								\
      nidrerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `nidrlex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX nidrlex (YYLEX_PARAM)
#else
# define YYLEX nidrlex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (nidrdebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (nidrdebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      nidr_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
nidr_symbol_value_print (FILE *nidroutput, int nidrtype, YYSTYPE const * const nidrvaluep)
#else
static void
nidr_symbol_value_print (nidroutput, nidrtype, nidrvaluep)
    FILE *nidroutput;
    int nidrtype;
    YYSTYPE const * const nidrvaluep;
#endif
{
  if (!nidrvaluep)
    return;
# ifdef YYPRINT
  if (nidrtype < YYNTOKENS)
    YYPRINT (nidroutput, nidrtoknum[nidrtype], *nidrvaluep);
# else
  YYUSE (nidroutput);
# endif
  switch (nidrtype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
nidr_symbol_print (FILE *nidroutput, int nidrtype, YYSTYPE const * const nidrvaluep)
#else
static void
nidr_symbol_print (nidroutput, nidrtype, nidrvaluep)
    FILE *nidroutput;
    int nidrtype;
    YYSTYPE const * const nidrvaluep;
#endif
{
  if (nidrtype < YYNTOKENS)
    YYFPRINTF (nidroutput, "token %s (", nidrtname[nidrtype]);
  else
    YYFPRINTF (nidroutput, "nterm %s (", nidrtname[nidrtype]);

  nidr_symbol_value_print (nidroutput, nidrtype, nidrvaluep);
  YYFPRINTF (nidroutput, ")");
}

/*------------------------------------------------------------------.
| nidr_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
nidr_stack_print (nidrtype_int16 *bottom, nidrtype_int16 *top)
#else
static void
nidr_stack_print (bottom, top)
    nidrtype_int16 *bottom;
    nidrtype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (nidrdebug)							\
    nidr_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
nidr_reduce_print (YYSTYPE *nidrvsp, int nidrrule)
#else
static void
nidr_reduce_print (nidrvsp, nidrrule)
    YYSTYPE *nidrvsp;
    int nidrrule;
#endif
{
  int nidrnrhs = nidrr2[nidrrule];
  int nidri;
  unsigned long int nidrlno = nidrrline[nidrrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     nidrrule - 1, nidrlno);
  /* The symbols being reduced.  */
  for (nidri = 0; nidri < nidrnrhs; nidri++)
    {
      fprintf (stderr, "   $%d = ", nidri + 1);
      nidr_symbol_print (stderr, nidrrhs[nidrprhs[nidrrule] + nidri],
		       &(nidrvsp[(nidri + 1) - (nidrnrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (nidrdebug)				\
    nidr_reduce_print (nidrvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int nidrdebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef nidrstrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define nidrstrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
nidrstrlen (const char *nidrstr)
#else
static YYSIZE_T
nidrstrlen (nidrstr)
    const char *nidrstr;
#endif
{
  YYSIZE_T nidrlen;
  for (nidrlen = 0; nidrstr[nidrlen]; nidrlen++)
    continue;
  return nidrlen;
}
#  endif
# endif

# ifndef nidrstpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define nidrstpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
nidrstpcpy (char *nidrdest, const char *nidrsrc)
#else
static char *
nidrstpcpy (nidrdest, nidrsrc)
    char *nidrdest;
    const char *nidrsrc;
#endif
{
  char *nidrd = nidrdest;
  const char *nidrs = nidrsrc;

  while ((*nidrd++ = *nidrs++) != '\0')
    continue;

  return nidrd - 1;
}
#  endif
# endif

# ifndef nidrtnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for nidrerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from nidrtname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
nidrtnamerr (char *nidrres, const char *nidrstr)
{
  if (*nidrstr == '"')
    {
      YYSIZE_T nidrn = 0;
      char const *nidrp = nidrstr;

      for (;;)
	switch (*++nidrp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++nidrp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (nidrres)
	      nidrres[nidrn] = *nidrp;
	    nidrn++;
	    break;

	  case '"':
	    if (nidrres)
	      nidrres[nidrn] = '\0';
	    return nidrn;
	  }
    do_not_strip_quotes: ;
    }

  if (! nidrres)
    return nidrstrlen (nidrstr);

  return nidrstpcpy (nidrres, nidrstr) - nidrres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
nidrsyntax_error (char *nidrresult, int nidrstate, int nidrchar)
{
  int nidrn = nidrpact[nidrstate];

  if (! (YYPACT_NINF < nidrn && nidrn <= YYLAST))
    return 0;
  else
    {
      int nidrtype = YYTRANSLATE (nidrchar);
      YYSIZE_T nidrsize0 = nidrtnamerr (0, nidrtname[nidrtype]);
      YYSIZE_T nidrsize = nidrsize0;
      YYSIZE_T nidrsize1;
      int nidrsize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *nidrarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int nidrx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *nidrfmt;
      char const *nidrf;
      static char const nidrunexpected[] = "syntax error, unexpected %s";
      static char const nidrexpecting[] = ", expecting %s";
      static char const nidror[] = " or %s";
      char nidrformat[sizeof nidrunexpected
		    + sizeof nidrexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof nidror - 1))];
      char const *nidrprefix = nidrexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int nidrxbegin = nidrn < 0 ? -nidrn : 0;

      /* Stay within bounds of both nidrcheck and nidrtname.  */
      int nidrchecklim = YYLAST - nidrn + 1;
      int nidrxend = nidrchecklim < YYNTOKENS ? nidrchecklim : YYNTOKENS;
      int nidrcount = 1;

      nidrarg[0] = nidrtname[nidrtype];
      nidrfmt = nidrstpcpy (nidrformat, nidrunexpected);

      for (nidrx = nidrxbegin; nidrx < nidrxend; ++nidrx)
	if (nidrcheck[nidrx + nidrn] == nidrx && nidrx != YYTERROR)
	  {
	    if (nidrcount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		nidrcount = 1;
		nidrsize = nidrsize0;
		nidrformat[sizeof nidrunexpected - 1] = '\0';
		break;
	      }
	    nidrarg[nidrcount++] = nidrtname[nidrx];
	    nidrsize1 = nidrsize + nidrtnamerr (0, nidrtname[nidrx]);
	    nidrsize_overflow |= (nidrsize1 < nidrsize);
	    nidrsize = nidrsize1;
	    nidrfmt = nidrstpcpy (nidrfmt, nidrprefix);
	    nidrprefix = nidror;
	  }

      nidrf = YY_(nidrformat);
      nidrsize1 = nidrsize + nidrstrlen (nidrf);
      nidrsize_overflow |= (nidrsize1 < nidrsize);
      nidrsize = nidrsize1;

      if (nidrsize_overflow)
	return YYSIZE_MAXIMUM;

      if (nidrresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *nidrp = nidrresult;
	  int nidri = 0;
	  while ((*nidrp = *nidrf) != '\0')
	    {
	      if (*nidrp == '%' && nidrf[1] == 's' && nidri < nidrcount)
		{
		  nidrp += nidrtnamerr (nidrp, nidrarg[nidri++]);
		  nidrf += 2;
		}
	      else
		{
		  nidrp++;
		  nidrf++;
		}
	    }
	}
      return nidrsize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
nidrdestruct (const char *nidrmsg, int nidrtype, YYSTYPE *nidrvaluep)
#else
static void
nidrdestruct (nidrmsg, nidrtype, nidrvaluep)
    const char *nidrmsg;
    int nidrtype;
    YYSTYPE *nidrvaluep;
#endif
{
  YYUSE (nidrvaluep);

  if (!nidrmsg)
    nidrmsg = "Deleting";
  YY_SYMBOL_PRINT (nidrmsg, nidrtype, nidrvaluep, nidrlocationp);

  switch (nidrtype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int nidrparse (void *YYPARSE_PARAM);
#else
int nidrparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int nidrparse (void);
#else
int nidrparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int nidrchar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE nidrlval;

/* Number of syntax errors so far.  */
int nidrnerrs;



/*----------.
| nidrparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
nidrparse (void *YYPARSE_PARAM)
#else
int
nidrparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
nidrparse (void)
#else
int
nidrparse ()

#endif
#endif
{
  
  int nidrstate;
  int nidrn;
  int nidrresult;
  /* Number of tokens to shift before error messages enabled.  */
  int nidrerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int nidrtoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char nidrmsgbuf[128];
  char *nidrmsg = nidrmsgbuf;
  YYSIZE_T nidrmsg_alloc = sizeof nidrmsgbuf;
#endif

  /* Three stacks and their tools:
     `nidrss': related to states,
     `nidrvs': related to semantic values,
     `nidrls': related to locations.

     Refer to the stacks thru separate pointers, to allow nidroverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  nidrtype_int16 nidrssa[YYINITDEPTH];
  nidrtype_int16 *nidrss = nidrssa;
  nidrtype_int16 *nidrssp;

  /* The semantic value stack.  */
  YYSTYPE nidrvsa[YYINITDEPTH];
  YYSTYPE *nidrvs = nidrvsa;
  YYSTYPE *nidrvsp;



#define YYPOPSTACK(N)   (nidrvsp -= (N), nidrssp -= (N))

  YYSIZE_T nidrstacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE nidrval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int nidrlen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  nidrstate = 0;
  nidrerrstatus = 0;
  nidrnerrs = 0;
  nidrchar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  nidrssp = nidrss;
  nidrvsp = nidrvs;

  goto nidrsetstate;

/*------------------------------------------------------------.
| nidrnewstate -- Push a new state, which is found in nidrstate.  |
`------------------------------------------------------------*/
 nidrnewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  nidrssp++;

 nidrsetstate:
  *nidrssp = nidrstate;

  if (nidrss + nidrstacksize - 1 <= nidrssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T nidrsize = nidrssp - nidrss + 1;

#ifdef nidroverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *nidrvs1 = nidrvs;
	nidrtype_int16 *nidrss1 = nidrss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if nidroverflow is a macro.  */
	nidroverflow (YY_("memory exhausted"),
		    &nidrss1, nidrsize * sizeof (*nidrssp),
		    &nidrvs1, nidrsize * sizeof (*nidrvsp),

		    &nidrstacksize);

	nidrss = nidrss1;
	nidrvs = nidrvs1;
      }
#else /* no nidroverflow */
# ifndef YYSTACK_RELOCATE
      goto nidrexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= nidrstacksize)
	goto nidrexhaustedlab;
      nidrstacksize *= 2;
      if (YYMAXDEPTH < nidrstacksize)
	nidrstacksize = YYMAXDEPTH;

      {
	nidrtype_int16 *nidrss1 = nidrss;
	union nidralloc *nidrptr =
	  (union nidralloc *) YYSTACK_ALLOC (YYSTACK_BYTES (nidrstacksize));
	if (! nidrptr)
	  goto nidrexhaustedlab;
	YYSTACK_RELOCATE (nidrss);
	YYSTACK_RELOCATE (nidrvs);

#  undef YYSTACK_RELOCATE
	if (nidrss1 != nidrssa)
	  YYSTACK_FREE (nidrss1);
      }
# endif
#endif /* no nidroverflow */

      nidrssp = nidrss + nidrsize - 1;
      nidrvsp = nidrvs + nidrsize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) nidrstacksize));

      if (nidrss + nidrstacksize - 1 <= nidrssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", nidrstate));

  goto nidrbackup;

/*-----------.
| nidrbackup.  |
`-----------*/
nidrbackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  nidrn = nidrpact[nidrstate];
  if (nidrn == YYPACT_NINF)
    goto nidrdefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (nidrchar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      nidrchar = YYLEX;
    }

  if (nidrchar <= YYEOF)
    {
      nidrchar = nidrtoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      nidrtoken = YYTRANSLATE (nidrchar);
      YY_SYMBOL_PRINT ("Next token is", nidrtoken, &nidrlval, &nidrlloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  nidrn += nidrtoken;
  if (nidrn < 0 || YYLAST < nidrn || nidrcheck[nidrn] != nidrtoken)
    goto nidrdefault;
  nidrn = nidrtable[nidrn];
  if (nidrn <= 0)
    {
      if (nidrn == 0 || nidrn == YYTABLE_NINF)
	goto nidrerrlab;
      nidrn = -nidrn;
      goto nidrreduce;
    }

  if (nidrn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (nidrerrstatus)
    nidrerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", nidrtoken, &nidrlval, &nidrlloc);

  /* Discard the shifted token unless it is eof.  */
  if (nidrchar != YYEOF)
    nidrchar = YYEMPTY;

  nidrstate = nidrn;
  *++nidrvsp = nidrlval;

  goto nidrnewstate;


/*-----------------------------------------------------------.
| nidrdefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
nidrdefault:
  nidrn = nidrdefact[nidrstate];
  if (nidrn == 0)
    goto nidrerrlab;
  goto nidrreduce;


/*-----------------------------.
| nidrreduce -- Do a reduction.  |
`-----------------------------*/
nidrreduce:
  /* nidrn is the number of a rule to reduce with.  */
  nidrlen = nidrr2[nidrn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  nidrval = nidrvsp[1-nidrlen];


  YY_REDUCE_PRINT (nidrn);
  switch (nidrn)
    {
        case 3:
/* #line 124 "nidrparse.y" */
    { nidr_reset(); }
    break;

  case 4:
/* #line 129 "nidrparse.y" */
    {}
    break;

  case 5:
/* #line 130 "nidrparse.y" */
    {}
    break;

  case 6:
/* #line 131 "nidrparse.y" */
    {str_lineno = nidrLineNumber;}
    break;

  case 7:
/* #line 132 "nidrparse.y" */
    {
					  fprintf(stderr,"\n\tunrecognized identifier '%s'\n", (nidrvsp[(2) - (5)].string));
					  fprintf(stderr,"\tinput line %d, within %s keyword.\n",
						str_lineno, nidr_keyword_name());
					  please_refer();
					}
    break;

  case 8:
/* #line 138 "nidrparse.y" */
    {str_lineno = nidrLineNumber;}
    break;

  case 9:
/* #line 139 "nidrparse.y" */
    {
					  fprintf(stderr,"\nMisplaced '%s' on input line %d.\n",
						(nidrvsp[(1) - (4)].string), str_lineno);
					  please_refer();
					}
    break;

  case 10:
/* #line 145 "nidrparse.y" */
    {
					  fprintf(stderr,
						"\nearly keyword termination, input line %d.\n",
						nidrLineNumber );
					  PRINT_DEBUG_NL;
					  nidr_signal_parse_error();
					  reset_nidrlex_state();
					}
    break;

  case 11:
/* #line 153 "nidrparse.y" */
    {
					fprintf(stderr,
						"\nUnrecognized keyword \"%s\", input line %d.\n",
						(nidrvsp[(1) - (1)].string), nidrLineNumber);
					}
    break;

  case 17:
/* #line 172 "nidrparse.y" */
    {}
    break;

  case 18:
/* #line 173 "nidrparse.y" */
    { nidr_bufr((nidrvsp[(1) - (1)].real),0); }
    break;

  case 19:
/* #line 174 "nidrparse.y" */
    { nidr_bufr((nidrvsp[(2) - (2)].real),1); }
    break;

  case 20:
/* #line 175 "nidrparse.y" */
    { nidr_bufr((nidrvsp[(2) - (2)].real),2); }
    break;

  case 21:
/* #line 176 "nidrparse.y" */
    { nidr_bufs((nidrvsp[(1) - (1)].qstring)); }
    break;


/* Line 1267 of yacc.c.  */
/* #line 1497 "nidr-parser.c" */
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", nidrr1[nidrn], &nidrval, &nidrloc);

  YYPOPSTACK (nidrlen);
  nidrlen = 0;
  YY_STACK_PRINT (nidrss, nidrssp);

  *++nidrvsp = nidrval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  nidrn = nidrr1[nidrn];

  nidrstate = nidrpgoto[nidrn - YYNTOKENS] + *nidrssp;
  if (0 <= nidrstate && nidrstate <= YYLAST && nidrcheck[nidrstate] == *nidrssp)
    nidrstate = nidrtable[nidrstate];
  else
    nidrstate = nidrdefgoto[nidrn - YYNTOKENS];

  goto nidrnewstate;


/*------------------------------------.
| nidrerrlab -- here on detecting error |
`------------------------------------*/
nidrerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!nidrerrstatus)
    {
      ++nidrnerrs;
#if ! YYERROR_VERBOSE
      nidrerror (YY_("syntax error"));
#else
      {
	YYSIZE_T nidrsize = nidrsyntax_error (0, nidrstate, nidrchar);
	if (nidrmsg_alloc < nidrsize && nidrmsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T nidralloc = 2 * nidrsize;
	    if (! (nidrsize <= nidralloc && nidralloc <= YYSTACK_ALLOC_MAXIMUM))
	      nidralloc = YYSTACK_ALLOC_MAXIMUM;
	    if (nidrmsg != nidrmsgbuf)
	      YYSTACK_FREE (nidrmsg);
	    nidrmsg = (char *) YYSTACK_ALLOC (nidralloc);
	    if (nidrmsg)
	      nidrmsg_alloc = nidralloc;
	    else
	      {
		nidrmsg = nidrmsgbuf;
		nidrmsg_alloc = sizeof nidrmsgbuf;
	      }
	  }

	if (0 < nidrsize && nidrsize <= nidrmsg_alloc)
	  {
	    (void) nidrsyntax_error (nidrmsg, nidrstate, nidrchar);
	    nidrerror (nidrmsg);
	  }
	else
	  {
	    nidrerror (YY_("syntax error"));
	    if (nidrsize != 0)
	      goto nidrexhaustedlab;
	  }
      }
#endif
    }



  if (nidrerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (nidrchar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (nidrchar == YYEOF)
	    YYABORT;
	}
      else
	{
	  nidrdestruct ("Error: discarding",
		      nidrtoken, &nidrlval);
	  nidrchar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto nidrerrlab1;


/*---------------------------------------------------.
| nidrerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
nidrerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label nidrerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto nidrerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (nidrlen);
  nidrlen = 0;
  YY_STACK_PRINT (nidrss, nidrssp);
  nidrstate = *nidrssp;
  goto nidrerrlab1;


/*-------------------------------------------------------------.
| nidrerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
nidrerrlab1:
  nidrerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      nidrn = nidrpact[nidrstate];
      if (nidrn != YYPACT_NINF)
	{
	  nidrn += YYTERROR;
	  if (0 <= nidrn && nidrn <= YYLAST && nidrcheck[nidrn] == YYTERROR)
	    {
	      nidrn = nidrtable[nidrn];
	      if (0 < nidrn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (nidrssp == nidrss)
	YYABORT;


      nidrdestruct ("Error: popping",
		  nidrstos[nidrstate], nidrvsp);
      YYPOPSTACK (1);
      nidrstate = *nidrssp;
      YY_STACK_PRINT (nidrss, nidrssp);
    }

  if (nidrn == YYFINAL)
    YYACCEPT;

  *++nidrvsp = nidrlval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", nidrstos[nidrn], nidrvsp, nidrlsp);

  nidrstate = nidrn;
  goto nidrnewstate;


/*-------------------------------------.
| nidracceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
nidracceptlab:
  nidrresult = 0;
  goto nidrreturn;

/*-----------------------------------.
| nidrabortlab -- YYABORT comes here.  |
`-----------------------------------*/
nidrabortlab:
  nidrresult = 1;
  goto nidrreturn;

#ifndef nidroverflow
/*-------------------------------------------------.
| nidrexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
nidrexhaustedlab:
  nidrerror (YY_("memory exhausted"));
  nidrresult = 2;
  /* Fall through.  */
#endif

nidrreturn:
  if (nidrchar != YYEOF && nidrchar != YYEMPTY)
     nidrdestruct ("Cleanup: discarding lookahead",
		 nidrtoken, &nidrlval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (nidrlen);
  YY_STACK_PRINT (nidrss, nidrssp);
  while (nidrssp != nidrss)
    {
      nidrdestruct ("Cleanup: popping",
		  nidrstos[*nidrssp], nidrvsp);
      YYPOPSTACK (1);
    }
#ifndef nidroverflow
  if (nidrss != nidrssa)
    YYSTACK_FREE (nidrss);
#endif
#if YYERROR_VERBOSE
  if (nidrmsg != nidrmsgbuf)
    YYSTACK_FREE (nidrmsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (nidrresult);
}


/* #line 181 "nidrparse.y" */


 void
nidrerror( char* s )
{
  PRINT_DEBUG_NL;
}

#include <setjmp.h>

static jmp_buf *nidr_jb;

 int
nidr_parse(const char *parser, FILE *df) {
	int rv;
	jmp_buf jb;

	nidr_reinit();
	nidr_jb = &jb;
	if (setjmp(jb))
		return 1;
	nidr_setup(parser, df);
	rv = nidrparse();
	nidr_jb = 0;
	rv += nidr_cleanup();
	nidr_reinit();
	return rv;
	}

 void
nidr_abort(void)
{
	if (nidr_jb)
		longjmp(*nidr_jb, 1);
	}

