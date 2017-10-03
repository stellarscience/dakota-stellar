/*********************************************************************
Copyright 2008, 2010 Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
retains certain rights in this software.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of Sandia Corporation nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
***********************************************************************/

%{

#define const

#include "nidr.h"
#include <string.h>

#include <stdlib.h>
#include <string.h>
extern int yylex(void);
#include <math.h>
#include <stdio.h>

extern int		yyLineNumber;
extern void		(*nidr_bufr)(Real, int);
extern void		(*nidr_bufs)(const char*);
extern int		nidr_cleanup(void);
extern const char*	nidr_keyword_name(void);
extern void		nidr_reinit(void);
extern void		nidr_reset(void);
extern void		nidr_setup(const char*, FILE*);
extern void		nidr_signal_parse_error(void);
extern void		reset_lex_state(void);
extern char		nidr_please_refer[];



#if  defined __STDC__  ||  defined __cplusplus
void	yyerror( char* s );
#else
void	yyerror();
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
	reset_lex_state();
	}

 static int str_lineno;
 static char garbled_fmt[] =
	"Unrecognized data.\n\tIdentifier value garbled for identifier '%s'.\n";
%}




%union
{
	KeyWord *keyword;
	KeyWord	*identifier;
	Real	real;
	int	integer;
	char*	string;
	char*	qstring;
}


%token	<keyword>	KEY_WORD
%token	<identifier>	IDENTIFIER
%token	<real>		REAL
%token	<string>	STRING
%token	<qstring>	QUOTED_STRING
%token			SEPARATOR
%token			EQUALS
%token	<string>	KEYWORDERROR
%token			END
%token			EXIT
%type	<keyword>	command






%start statements


%%

statements	:
				| statements statement { nidr_reset(); }
				;



statement	: END {}
				| command END {}
				| command STRING {str_lineno = nidrLineNumber;} error END
					{
					  fprintf(stderr,"\n\tunrecognized identifier '%s'\n", $2);
					  fprintf(stderr,"\tinput line %d, within %s keyword.\n",
						str_lineno, nidr_keyword_name());
					  please_refer();
					}
				| STRING {str_lineno = nidrLineNumber;}  error END
					{
					  fprintf(stderr,"\nMisplaced '%s' on input line %d.\n",
						$1, str_lineno);
					  please_refer();
					}
				| error END
					{
					  fprintf(stderr,
						"\nearly keyword termination, input line %d.\n",
						yyLineNumber );
					  PRINT_DEBUG_NL;
					  nidr_signal_parse_error();
					  reset_lex_state();
					}
				| KEYWORDERROR {
					fprintf(stderr,
						"\nUnrecognized keyword \"%s\", input line %d.\n",
						$1, yyLineNumber);
					}
				;


command		: KEY_WORD
				| KEY_WORD data
				| KEY_WORD EQUALS data
				;


data		: datum
			| data datum
			;


datum		: SEPARATOR {}
			| REAL			{ nidr_bufr($1,0); }
			| '*' REAL		{ nidr_bufr($2,1); }
			| ':' REAL		{ nidr_bufr($2,2); }
			| QUOTED_STRING		{ nidr_bufs($1); }
			| IDENTIFIER
			| IDENTIFIER EQUALS
			;

%%

 void
yyerror( char* s )
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
