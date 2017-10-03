%{
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

#include "nidr.h"
#include "nidr-parser.h"
#include <stdlib.h>	/* possibly redundant, but necessary for strtod on some systems */

static int my_yy_input(char*, int);
extern char *nidr_KWscopy(const char*);

#ifdef IDR_DEBUG
static void print_token( char*, char*, char* );
#else
#define print_token(a,b,c) /*nothing*/
#endif


#define YY_INPUT(buf,result,maxsize) result = my_yy_input(buf,maxsize)

int	yyLineNumber = 1;

extern void	nidr_abort(void);
extern void	(*nidr_comment)(const char*);
extern KeyWord *(*nidr_identifier)(const char*);
extern KeyWord *nidr_keyword(const char*);
extern const char *nidr_keyword_name(void);
/* Hardwire the now rarely used help text.  It truly is something the
   application should provide, not NIDR, but had trouble with the
   declaration/definition w.r.t. MSVS when building DAKOTA source.
   Ultimately would want to pass it into nidr_parse(...) 
extern char	nidr_please_refer[];
*/
char nidr_please_refer[] = 
  "Please refer to the input specification reference for this application.";
extern void	nidr_tolower(char*);
static void	bad_number(const char*);
static void	unput_string(char *, int);
%}

s	([+-]?)
d	([0-9])
e	([eE]{s}{d}+)
integer	{s}{d}+([a-zA-Z_+\-.]+)?
NaN	[nN][aA][nN]
Inf	[iI][nN][fF]([iI][nN][iI][tT][yY])?
real	{s}(((({d}+\.{d}*{e}?)|({d}*\.{d}+{e}?)|({d}+{e}))([a-zA-Z_+\-.]+)?)|{NaN}|{Inf})
word	[a-zA-Z][0-9a-zA-Z_@~.-]*
cmnt	[ \t]*#[^\\\n]*
cmnt_e	{cmnt}\\[ \t]*
nl	[ \t]*\\[ \t]*\n
sep	[ \t]*\,
eq	[ \t]*=[ \t]*
ws	[ \t]+
exit	[eE][xX][iI][tT]
end	[;]

%s DATA
%s DQ_STRING
%s Q_STRING

%%

%{
	KeyWord	*identifier, *keyword;
	char *se;
%}

<INITIAL>{exit}	{
			print_token("exit", NULL, "<INITIAL>" );
			return 0;
			}

<INITIAL,DATA>\n {
			if (nidr_comment)
				nidr_comment(yytext);
			yyLineNumber++;
			}

<INITIAL>{word}	{
			nidr_tolower(yytext);
			keyword = nidr_keyword(yytext);

			if  (keyword) {
				print_token( "keyword", (char*) yytext, "<INITIAL>" );
				if (keyword->kind & KWKind_Dynmult)
					return END;
				yylval.keyword = keyword;
				BEGIN DATA;
				return KEY_WORD;
				}
			print_token( "error", (char*) yytext, "<INITIAL>" );
			yylval.string = yytext;
			return KEYWORDERROR;
			}

<INITIAL>#*{end} {
			if (*yytext == '#') {
				if (nidr_comment)
					nidr_comment(yytext);
				}
			else {
				print_token( "end", NULL, "<INITIAL>" );
				return END;
				}
			}

<INITIAL>{ws}		{ print_token( "whitespace", (char*) yytext, "<INITIAL>" ); }

<INITIAL>[^;#] {
			print_token( "word", (char*) yytext, "<INITIAL>" );
			yylval.string = nidr_KWscopy(yytext);
			return STRING;
			}
<DATA>{eq} {
			print_token( "equals", NULL, "<DATA>" );
			return EQUALS;
			}
<DATA>{sep} {
			print_token( "separator", NULL, "<DATA>" );
			return SEPARATOR;
			}
<DATA>{integer}	{
			print_token( "integer", (char*) yytext, "<DATA>" );
			/* The following, old way has trouble with */
			/* big floating-point numbers written without */
			/* a decimal point */
			/* yylval.integer = (int) atoi( (char*) yytext ); */
			/* return INTEGER; */
			yylval.real = strtod(yytext,&se);
			if (*se)
				bad_number(yytext);
			return REAL;
			}
<DATA>{real} {
			print_token( "real", (char*) yytext, "<DATA>" );
			yylval.real = strtod(yytext,&se);
			if (*se)
				bad_number(yytext);
			return REAL;
			}
<DATA>{word} {
			nidr_tolower(yytext);
			if (nidr_keyword(yytext)) {
				unput_string(yytext, yyleng);
				BEGIN INITIAL;
				return END;
				}
			identifier = nidr_identifier(yytext);

			if (identifier) {
				print_token( "identifier", (char*) yytext, "<DATA>");

				yylval.identifier = identifier;
				return IDENTIFIER;
				}
			else {
				print_token( "word", (char*) yytext, "<DATA>" );
				fprintf(stderr,"\n\tunrecognized identifier '%s'\n", yytext);
				fprintf(stderr,"\tinput line %d, within %s keyword.\n",
						yyLineNumber, nidr_keyword_name());
				fprintf(stderr, "%s\n", nidr_please_refer);
				nidr_abort();
				yylval.string = nidr_KWscopy(yytext);
				return STRING;
				}
			}

<DATA>{ws}		{ print_token( "whitespace", (char*) yytext, "<DATA>" ); }

<DATA>'	{
			print_token( "opening quote", NULL, "<DATA>" );
			BEGIN Q_STRING;
			}
<DATA>{end} {
			print_token( "end", NULL, "<DATA>" );
			BEGIN INITIAL;
			return END;
			}
<DATA>\" {
			print_token( "double quote", NULL, "<DATA>" );
			BEGIN DQ_STRING;
			}
<DATA>[^#;*:]	{
			print_token( "word", (char*) yytext, "<DATA>" );
			fprintf(stderr, "Input line %d:  '%s' unexpected.\n",
				yyLineNumber, yytext);
			nidr_abort();
			yylval.string = nidr_KWscopy(yytext);
			return STRING;
			}
<Q_STRING>'' {
			print_token( "q:repeated quote", NULL, "<Q_STRING>" );
			yytext[--yyleng] = 0;
			yymore();
			}
<Q_STRING>[^'\n]* {
			yymore();
			print_token( "q:text", (char*) yytext, "<Q_STRING>" );
			}
<Q_STRING>('|\n) {
			print_token( "closing quote", NULL, "<Q_STRING>" );
			if (yytext[--yyleng] == '\n')
				unput('\n');
			yytext[yyleng] = 0;
			print_token( "quoted string", (char*) yytext, "<Q_STRING>" );
			yylval.qstring = nidr_KWscopy(yytext);
			BEGIN DATA;
			return QUOTED_STRING;
			}
<DQ_STRING>\"\" {
			print_token( "q:repeated double quote", NULL, "<DQ_STRING>" );
			yytext[--yyleng] = 0;
			yymore();
			}
<DQ_STRING>[^\"\n]* {
			yymore();
			print_token( "q:text", (char*) yytext, "<DQ_STRING>" );
			}
<DQ_STRING>(\"|\n) {
			print_token( "closing double quote", NULL, "<DQ_STRING>" );
			if (yytext[--yyleng] == '\n')
				unput('\n');
			yytext[yyleng] = 0;
			print_token( "quoted string", (char*) yytext, "<DQ_STRING>" );
			yylval.qstring = nidr_KWscopy(yytext);
			BEGIN DATA;
			return QUOTED_STRING;
			}
{nl} {
			print_token( "newline", NULL, NULL );
			if (nidr_comment)
				nidr_comment(yytext);
			yyLineNumber++;
			}
{cmnt_e} {
			if (nidr_comment)
				nidr_comment(yytext);
			print_token( "comment - escaped newline", (char*) yytext, "<INITIAL>" );
			unput( '\\' );
			}
{cmnt}			{
			if (nidr_comment)
				nidr_comment(yytext);
			print_token( "comment", (char*) yytext, "<INITIAL>" );
			}
<DATA>[:*]		{ return yytext[0]; }

%%

 static int nidr_ateof = 0;
 static const char *nidr_instring = 0;

 void
bad_number(const char *s)
{
	fprintf(stderr, "\n*** Invalid number \"%s\" on input line %d within "
			"%s keyword.\n", s, yyLineNumber, nidr_keyword_name());
	nidr_abort();
	}

 void
nidr_set_input_string(const char *s)
{
	nidr_instring = s;
	}

 static int
my_yy_input(char *buf, int bmax)
{
	int c, n;

	if (nidr_ateof)
		return 0;
	n = 0;
	if (nidr_instring) {
		while(n < bmax && (c = *nidr_instring)) {
			++nidr_instring;
			if (c != '\r')
				buf[n++] = c;
			if (c == '\n')
				break;
			}
		}
	else if (nidr_current_buffer->nidr_is_interactive)
		while(n < bmax) {
		if ((c = getc(yyin)) == EOF)
			break;
		if (c != '\r')
			buf[n++] = c;
		if (c == '\n')
			break;
		}
	else
		while(n < bmax) {
		if ((c = getc(yyin)) == EOF)
			break;
		if (c != '\r')
			buf[n++] = c;
		}
	if (n <= 0) {
		nidr_ateof = n = 1;
		if (c == EOF) {
			fclose(yyin);
			yyin = NULL;
			}
		buf[0] = ';';
		}
	return n;
	}

 void
reset_lex_state(void)
{
	BEGIN INITIAL;
	nidr_abort();
	}

 void
unput_string(char *s, int len)
{
	while(--len >= 0)
		unput(s[len]);
	}

 void
nidr_reinit(void)
{
	nidr_init = 1;
	nidr_ateof = 0;
	if (nidr_current_buffer) {
		free(nidr_current_buffer->nidr_ch_buf);
		free(nidr_current_buffer);
		nidr_current_buffer = 0;
		}
	}

 int
yywrap(void)
{ return 1; }
