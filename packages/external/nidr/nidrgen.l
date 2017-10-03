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

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nidr.h"	/* for KeyWordKind */

#ifndef NO_NIDR_DYNLIB
typedef KeyWord *(*KW_ADD)(void);
#ifdef _WIN32
#include <windows.h>
#define dlopen(x,y) LoadLibrary(x)
#define find_dlsym(a,b,c) (a = (KW_ADD)GetProcAddress((HINSTANCE)(b),c))
#define dlclose(x) FreeLibrary((HMODULE)x)
#define NO_DLERROR
#else
#include <dlfcn.h>
#define find_dlsym(a,b,c) (a = (KW_ADD)dlsym(b,c))
#undef NO_DLERROR
#endif
#endif

#ifndef DMBLK_GULP
#define DMBLK_GULP 1024
#endif
#ifndef DMLBK_SGULP
#define DMBLK_SGULP 16000
#endif
#ifndef MBLK_GULP
#define MBLK_GULP 8191
#endif
#ifndef KWSTACKLEN
#define KWSTACKLEN 100	/*should be overkill*/
#endif
#ifndef DEFTAGLEN
#define DEFTAGLEN 2048	/* should be overkill*/
#endif

 typedef union
YYSTYPE { char *s; } YYSTYPE;

 enum LastSeen {
	Saw_kwname	= 1,
	Saw_type	= 2,
	Saw_func	= 3,
	Saw_dflt	= 4,
	Saw_dfltname	= 5,
	Saw_lt		= 6,
	Saw_le		= 7,
	Saw_ge		= 8,
	Saw_gt		= 9,
	Saw_desc	= 10,
	Saw_group	= 11,
	Saw_tag		= 12,
	Saw_len		= 13,
	Saw_dy		= 14
	};

 static const char *tokname[] = {
	"<bug>",
	"kwname",
	"type",
	"{funcspec}",
	":=",
	"#define name for :=",
	"<",
	"<=",
	">=",
	">",
	"DESC",
	"GROUP",
	"TAG",
	"LEN"
	};

 enum ParState {
	/* parstate values... */
	PPar_emtpy	= 0,
	PPar_starting	= 1,
	PPar_started	= 2,
	PPar_leftpar	= 3,
	Dupchecking	= 4,
	Dupchecked	= 5,
	Printing	= 6,
	KwStacked	= 7,
	Sorted		= 8,
	Leftput		= 9,
	Parchecked	= 10
	};

 /* values for lastseen also include single characters ( [ { } ] ) ; */

 /* enums to make constants visible to debuggers */
 enum {
	DMblk_gulp	= DMBLK_GULP,
	Dmblk_sgulp	= DMBLK_SGULP,
	Mblk_gulp	= MBLK_GULP,
	KW_stacklen	= KWSTACKLEN,
	Br_stacklen	= 2*KW_stacklen,
	DHTlen		= 997,
	DHTfactor	= 43,
	HTlen		= 127,
	KWPlen		= 32,
	KWOUTlen0	= 128 };

 typedef struct Defhash Defhash;
 typedef struct DefhashHead DefhashHead;
 typedef struct DefhashTab DefhashTab;
 typedef struct HashHead HashHead;
 typedef struct KWctx KWctx;
 typedef struct Mblk Mblk;
 typedef struct Taghash Taghash;
 typedef struct StrList StrList;

 struct
DefhashHead {
	DefhashHead *hnext, *thread;
	char *name;
	};

 struct
Defhash {
	DefhashHead hh;
	char *value;
	int qkeep;
	};

 struct
Taghash {
	DefhashHead hh;
	KWctx *kw;
	int uses;
	};

 struct
DefhashTab {
	DefhashHead *Tab[DHTlen];
	DefhashHead *lastdh;
	size_t nalloc;	/* size of new entry allocations in units of sizeof(void*) */
	int nent;
	};

 struct
Mblk {
	Mblk *next;
	void *stuff[Mblk_gulp];
	};

 struct
HashHead { HashHead *hnext, *hprev; };

 struct
StrList {
	HashHead h;
	StrList *next;
	char *val;
	int Lineno;
	int primary;
	int seqno;
	int needrev;
	};

 struct
KWctx {
	KWctx *next;
	StrList *names;	/* list of name and aliases */
	char *name1;	/* first name (for tagout) */
	KWctx *kw;	/* list of contained keywords */
	KWctx *master;	/* non-null ==> this is an alias; use master's kwkno */
	KWctx *pmaster;	/* used in computing paoff */
	KWctx *pparent;	/* paren parent, for ((...)|(...) ...) */
	KWctx *rparent;	/* "real" parent, for stuff inside ((...)) */
	StrList *funcs;	/* pre- and post order functions and (void*) args */
	char *Br_top;	/* Br_top when this keyword began */
	char *lb;	/* lower bound (if kind & KWKind_Lb is nonzero) */
	char *ub;	/* upper bound (if kind & KWKind_Ub is nonzero) */
	char *init;	/* initial value (if kind & KWKind_init is nonzero) and */
			/* kind & KWKind_Mask is 1 or 2 (KWKind_Int or KWKind_Real) */
	char *cinit;	/* initial value for STRING-valued keyword,	*/
			/* ==> (kind & (KWKind_Mask | KWKind_init))	*/
			/*	== (KWKind_Str | KWKind_init)		*/
	char *defname;	/* name of preprocessor token to be #defined */
	char *desc;
	char *dylib;
	char *group;
	char *tag;
	char *len;	/* keyword giving array length */
	size_t dtnext;	/* for adjusting deftag when this keyword goes out of scope */
	size_t name1len;
	int alt;	/* alternative group number of this keyword */
	int kind;	/* kind of values for this keyword */
	int kno;	/* number of this keyword, for generating kw_n names */
	int kwkno;	/* number of this keyword's kw array (if any) */
	int kwknop;	/* if this is a pparent and has a pparent, kwkno for this and parents */
	int nalt;	/* number of alternative groups in contained keywords */
	int nfcn;	/* number of function and arg entries in funcs */
	int nkw;	/* number of contained keywords */
	int nreq;	/* number of required elements in contained keywords */
	int req;	/* required group for this keyword (0 if optional) */
	enum ParState parstate;	/* state in handling ((...)) */
	int altoffset;	/* for expand */
	int reqoffset;	/* for expand */
	int level;
	int objno;	/* for kwcomp2 */
	int agroup;	/* for GuiKeyWord */
	int paoff;	/* offset of primary keyword */
	};

 static int lastagroup, lastobjno, nalias, seqno;

 Defhash *lastdef;
 DefhashTab DHTab, TagTab;
 FILE *mtagout, *tagin, *tagout;
 HashHead HTab[HTlen];
 KWctx KWbase, *lastkw, *lastkw1;
 KWctx *KWStack[KW_stacklen], **KWtop = KWStack;
 Mblk FirstMblk, *CurMblk = &FirstMblk;
 YYSTYPE yyval;
 char Br_stack[Br_stacklen], *Br_top = Br_stack;
 char deftag[DEFTAGLEN], *dtend = deftag + DEFTAGLEN - 2, *dtnext = deftag;
 char *infname, *lasttag, *progname, *specfile, *tagfname;
 int Lineno = 1, dfltgroup = 1, startwarn = 1, wantfuncs = 1;
 int brace, btype, dupnames, dynlibseen, dynline, expand, guikeywds, ignoretags;
 int kwmult, lastkno, lastseen, liberror, libexpand, missing, needcomma, nsquawk;
 int or_mode, saweof;
 size_t mbavail = Mblk_gulp;
 void **mbnext = FirstMblk.stuff;
 static void Eat_comment(void);
 static void Keywd(const char *);
 static void Saw_num(const char*);
 static void Saw_quote(void);
 static void Saw_tok(int, int);

 /* Possible values for *Br_top:
 //	 0  ==> empty (for Br_stack[0] only)
 //	'(' ==> open paren seen, awaiting ')'
 //	'[' ==> open square bracket, awaiting ']'
 */

 static void
botch(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	fprintf(stderr, "\n%s:\n\tbotch on line %d", progname, Lineno);
	if (infname)
		fprintf(stderr, " of \"%s\"", infname);
	fprintf(stderr, ":\n\t");
	vfprintf(stderr, fmt, ap);
	fputs(".\n", stderr);
	va_end(ap);
	exit(1);
	}

 static void
squawk(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	fprintf(stderr, "\n%s: line %d", progname, Lineno);
	if (infname)
		fprintf(stderr, " of \"%s\"", infname);
	fprintf(stderr, ":\n\t");
	vfprintf(stderr, fmt, ap);
	fputs(".\n", stderr);
	va_end(ap);
	++nsquawk;
	}

 static void*
Malloc(size_t L, const char *where)
{
	void *rv = malloc(L);
	if (!rv)
		botch("malloc(%lu) failure in %s()", (unsigned long)L);
	return rv;
	}

 static void*
Realloc(void *rv, size_t L, const char *where)
{
	rv = realloc(rv, L);
	if (!rv)
		botch("ralloc(%lu) failure in %s()", (unsigned long)L);
	return rv;
	}

 void
unexpected(const char *what)
{ botch("unexpected \"%s\"", what); }

 static void *
Alloc(size_t len)
{
	Mblk *mb;
	size_t L;
	void *rv;

	L = (len + sizeof(void*) - 1)/sizeof(void*);
	if (L >= Mblk_gulp)
		botch("Alloc(%lu) failure", (unsigned long)len);
	if (L > mbavail) {
		if (!(mb = CurMblk->next)) {
			mb = (Mblk*)Malloc(sizeof(Mblk), "Alloc");
			CurMblk->next = mb;
			}
		CurMblk = mb;
		mbavail = Mblk_gulp;
		mbnext = mb->stuff;
		mb->next = 0;
		memset(mbnext, 0, sizeof(mb->stuff));
		}
	rv = mbnext;
	mbnext += L;
	mbavail -= L;
	return rv;
	}

 StrList *
new_StrList(const char *s, StrList *nxt)
{
	size_t L = strlen(s) + 1;
	StrList *rv = (StrList*)Alloc(sizeof(StrList) + L);
	strcpy(rv->val = (char*)(rv + 1), s);
	if ((rv->next = nxt))
		rv->needrev = 1;
	rv->Lineno = Lineno;
	return rv;
	}

 static void
Alias(const char *s)
{
	KWctx *kw;
	StrList *sl;

	if (lastseen != Saw_kwname)
		botch("ALIAS out of place");
	kw = lastkw;
	for(s += 6; *s <= ' '; s++);
	sl = kw->names = new_StrList(s, kw->names);
	++nalias;
	sl->seqno = ++seqno;
	lastseen = Saw_kwname;
	}

 static void
reqchk(const char *a, const char *s)
{
	if (a)
		while(*a)
			if (*a++ == lastseen)
				return;
	botch("\"%s\" out of place", s);
	}

 static void
reqchk1(const char *a, const char *s)
{
	int c;

	if (*Br_top != *s) {
		if ((c = *Br_top)) {
			switch(c) {
			 case '[': c = ']'; break;
			 case '(': c = ')'; break;
			 default:  c = '?';
			 }
			botch("Missing '%c'", c);
			}
		else
			a = 0;
		}
	reqchk(a, s+1);
	}

 static void
Bar(void)
{
	static char allowed[] = { /*(*/ ')', Saw_kwname, Saw_type, 0 };
	reqchk(KWtop - KWStack > 1 || KWtop[0]->kw ? allowed : "", "|");
	lastseen = '|';
	lastkw = lastkw1 = 0;
	or_mode = 1;
	}

 static void
Lpar(void)
{
	KWctx *kw;
	int olastseen = lastseen;
	static char allowed[] = { '(', '[', '|', ']', ')', Saw_kwname, Saw_type, 0 };

	reqchk(allowed, "("/*)*/);
	lastseen = '('; /*)*/
	if (KWtop > KWStack) {
		kw = KWtop[0];
		if (or_mode && kw->parstate == PPar_leftpar && !kw->kw) {
			dtnext = deftag + kw->dtnext;
			--KWtop;
			Keywd(0);
			kw->pparent = KWtop[0];
			}
		else if (olastseen == '(' || olastseen == '[' /*])*/)
			Keywd(0);
		lastseen = '(' /*)*/;
		}
	if (++Br_top - Br_stack >= Br_stacklen)
		botch("Br_stack overflow");
	*Br_top = '('; /*)*/
	lastkw = lastkw1 = 0;
	}

 static void
Lsqb(void)
{
	static char allowed[] = { /*([*/ ']', ')', Saw_kwname, Saw_type, 0 };
	reqchk(allowed, "["/*]*/);
	if (++Br_top - Br_stack >= Br_stacklen)
		botch("Br_stack overflow");
	lastseen = *Br_top = '['; /*]*/
	lastkw = lastkw1 = 0;
	}

 static void
Rpar(void)
{
	static char allowed[] =  {/*[(*/ ')', ']', Saw_kwname, Saw_type, 0 };
	reqchk1(allowed, "()");
	if (KWtop <= KWStack)
		botch("bug: KW stack error");
	--Br_top;
	if (!KWtop[0]->kw && KWtop - KWStack > 2
	 && KWtop[-1]->parstate == PPar_starting)
		KWtop[-1]->Br_top = Br_top;
	--KWtop;
	lastseen = /*(*/ ')';
	lastkw = lastkw1 = 0;
	}

 static void
Rsqb(void)
{
	static char allowed[] = {/*(*/')', '[', ']', Saw_kwname, Saw_type, 0 };
	reqchk1(allowed, "[]");
	if (KWtop <= KWStack)
		botch("bug: KW stack error");
	--Br_top;
	--KWtop;
	lastseen = /*[*/ ']';
	lastkw = lastkw1 = 0;
	}

 static void
Lbrace(void)
{
	if (brace++ || (lastseen != Saw_kwname && lastseen != Saw_type))
		unexpected("{");
	lastseen = '{'; /*}}*/
	}

 static void
Rbrace(void)
{
	if (--brace || (lastseen != Saw_func && lastseen != ';'))
		/*{{*/unexpected("}");
	lastseen = Saw_type;
	}

 static void
Type(int n, const char *s)
{
	if (lastseen != Saw_kwname || !lastkw)
		botch("type %s out of place", s);
	lastkw->kind |= n;
	lastkw = 0;
	lastseen = Saw_type;
	}

 static char *
def_name(const char *s0)
{
	char *s;
	const char *s1;
	size_t L, L1;
	static char *slast, *snext;

	s1 = s0;
	while(*s1++);
	L = s1 - s0;
	if (snext + L > slast) {
		snext = (char*)Malloc(L1 = Dmblk_sgulp + L, "def_name");
		slast = snext + L1;
		}
	strcpy(s = snext, s0);
	snext += L;
	return s;
	}

static void **DHHlast, **DHHnext;	/* Could be private to DHHlookup, but here */
					/* to make them visible to debuggers. */
 static DefhashHead *
DHHlookup(const char *s0, DefhashTab *Htab, int add)
{
	DefhashHead *h, **hp;
	const char *s;
	size_t L;
	unsigned int c, x;

	for(x = 0, s = s0; (c = *(unsigned char*)s); ++s)
		x += DHTfactor*x + c;
	x %= DHTlen;
	for(hp = &Htab->Tab[x]; (h = *hp); hp = &h->hnext) {
		if (!strcmp(s0, h->name))
			return h;
		}
	if (!add)
		return h;
	Htab->nent++;
	if (DHHnext + Htab->nalloc >= DHHlast) {
		DHHnext = (void**)Malloc(L = DMblk_gulp*sizeof(void*), "DHHlookup");
		DHHlast = DHHnext + DMblk_gulp;
		memset(DHHnext, 0, L);
		}
	*hp = h = (DefhashHead*)DHHnext;
	DHHnext += Htab->nalloc;
	h->name = def_name(s0);
	h->thread = Htab->lastdh;
	return Htab->lastdh = h;
	}

 static void
bothplaces(KWctx *kw, const char *name)
{ squawk("For %s, %s appears both in specfile and tagfile", kw->names->val, name); }

 static void
tagimport(KWctx *kw, KWctx *tw)
{
	typedef struct KWimport KWimport;
	struct KWimport { const char *name; size_t off; };
	static KWimport Imp[] = {
		{ "lower bound",	offsetof(KWctx,lb) },
		{ "upper bound",	offsetof(KWctx,ub) },
		{ "default (:= value)",	offsetof(KWctx,init) },
		{ "default (:= value)",	offsetof(KWctx,cinit) },
		{ "defname",		offsetof(KWctx,defname) },
		{ "DESC",		offsetof(KWctx,desc) },
		{ "GROUP",		offsetof(KWctx,group) },
		{ "LEN",		offsetof(KWctx,len) },
		{0}};
	KWimport *ki;
#define field(k,n) *(char**)((char*)k+n)
	for(ki = Imp; ki->name; ++ki) {
		if (field(tw,ki->off)) {
			if (field(kw,ki->off))
				bothplaces(kw, ki->name);
			else
				field(kw,ki->off) = field(tw,ki->off);
			}
		}
#undef field
	}

#ifndef NO_NIDR_DYNLIB

static const char *typename[8] = {
	"no value", "INTEGER", "REAL", "STRING",
	"(Bug! no value LIST)", "INTEGERLIST", "REALLIST", "STRINGLIST" };

 static int
kwseqnocomp(const void *a, const void *b)
{
	const KeyWordx *ka, *kb;

	ka = *(KeyWordx**)a;
	kb = *(KeyWordx**)b;
	return (int)ka->seqno - (int)kb->seqno;
	}

 static void
Lsqbpar(void)
{
	Lsqb();
	Lpar();
	}

 static void
Rparsqb(void)
{
	Rpar();
	Rsqb();
	}

 static void
Nobrak(void) {}

 typedef struct
KW_Cstak { KeyWordx **Cstak; Uint n, nmax; } KW_Cstak;

 static void kw_finish(int);
 static void kw_insert(KeyWord*, KW_Cstak*, int);

 static KeyWord*
kw_stakadj(KeyWord *kw, KW_Cstak *Stk, int needbar, int n0, KeyWord *okw3)
{
	KeyWordx *kw1, *kw2, **x;
	Uint i, j, k, nst, nz;

	kw1 = (KeyWordx*)kw;
	nz = 0;
	while(!kw1[++nz].kw.name);
	if (nz == 1 && *kw1[nz].kw.name)
		okw3 = 0;
	nst = Stk->n;
	x = Stk->Cstak;
	kw2 = (KeyWordx*)kw1->kw.kw;
	for(i = 0; i < nst; ++i) {
		if (x[i]->kw.kw == &kw2->kw) {
			j = i;
			++kw1;
			/* No need to test kw1->kw.name in the following loop */
			/* as the == test will fail if kw1->kw.name is nonzero. */
			while(++j < nst && x[j]->kw.kw == kw1->kw.kw)
				++kw1;
			k = nz - (j-i);
			goto have_jk;
			}
		}
	k = nz;
	j = n0;
 have_jk:
	while(nst > j) {
		kw_insert(&x[--nst]->kw, Stk, 0);
		Rpar();
		}
	if (needbar)
		Bar();
	if (k) {
		if (nst + k > Stk->nmax)
			x = Stk->Cstak = (KeyWordx**)Realloc(x,
				(Stk->nmax = 2)*sizeof(KeyWordx*), "kw_stakadj");
		do {
			if (kw1->kw.kw != okw3)
				Lpar();
			x[j++] = kw1++;
			}
			while(--k);
		nst = j;
		}
	Stk->n = nst;
	if (!*kw1->kw.name)
		return 0;
	return &kw1->kw;
	}

 static int
same_funcs(KeyWordx *kw1, KeyWordx *kw2)
{
	KeyWordx *kw3, *kw4;

	if (!kw1->kw.nkw || !kw2->kw.nkw)
		return 0;
	kw3 = (KeyWordx*)kw1->kw.kw;
	kw4 = (KeyWordx*)kw2->kw.kw;
	if (kw3 == kw4)
		return 1;
	if (kw4->kw.name)
		return 0;
	kw4 = (KeyWordx*)kw4->kw.kw;
	return kw3 == kw4;
	}

 static char *
StrCpy(const char *s)
{
	char *t;
	size_t L;

	if (!s)
		return 0;
	L = strlen(s);
	t = (char*)Alloc(L+1);
	strcpy(t,s);
	return t;
	}

 static int Saw_name(const char*);
 static void Saw_stok(int, size_t);

 static void
kw_insert(KeyWord *kw, KW_Cstak *Stk, int top)
{
	enum { NKW0 = 32 };
	KeyWord *kw3, *kw30, *okw3;
	KeyWordx *kw0, *kw1, *kw2, *kw4, **x, *x0[NKW0];
	Uint i, j, k, knd, nkw, nst0, nst1;
	char buf[32];
	int alt2, incc, needbar, needbar1, samef0, samef1;
	typedef void (*Brak)(void);
	Brak Lbrak, Rbrak;

	kw1 = (KeyWordx*)kw->kw;
	while(!kw1->kw.name)
		++kw1;
	nkw = kw->nkw;
	x = x0;
	if (nkw > NKW0)
		x = (KeyWordx**)Malloc(nkw * sizeof(KeyWordx*), "kw_insert");
	for(i = 0; i < nkw; ++i)
		x[i] = kw1++;
	if (nkw > 1)
		qsort(x, nkw, sizeof(KeyWordx*), kwseqnocomp);
	kw0 = 0;
	kw30 = 0;
	nst0 = Stk->n;
	for(i = j = samef0 = 0; i < nkw; i = j) {
		kw1 = x[i];
		kw2 = 0;
		while(++j < nkw) {
			kw2 = x[j];
			if (!kw2->kw.paoff)
				break;
			kw2 = 0;
			}
		Lbrak = Rbrak = 0;
		okw3 = kw30;
		if (top) {
			kw_finish(kw1->kw.kind & KWKind_12);
			kw3 = kw1->kw.kw;
			Rbrak = Nobrak;
			goto nobrak;
			}
		needbar = needbar1 = kw0 && kw0->kw.alt == kw1->kw.alt;
		nst1 = Stk->n;
		alt2 = incc = 0;
		kw4 = 0;
		if (kw2) {
			alt2 = kw1->kw.alt == kw2->kw.alt;
			kw4 = (KeyWordx*)kw2->kw.kw;
			}
		if ((kw3 = kw30 = kw1->kw.kw)) {
			if (!kw3->name) {
				kw3 = kw_stakadj(kw3, Stk, needbar, nst0, okw3);
				needbar1 = 0;
				}
			else
				while(Stk->n > nst0) {
					if (Stk->Cstak[Stk->n-1]->kw.kw == kw3) {
						if (kw4 && kw4->kw.kw == kw3)
							incc = 1; /* in contained context */
						break;
						}
					kw_insert(&Stk->Cstak[--Stk->n]->kw, Stk, 0);
					Rpar();
					}
			}
		else while(Stk->n > nst0) {
			kw_insert(&Stk->Cstak[--Stk->n]->kw, Stk, 0);
			Rpar();
			}
		samef0 = needbar && same_funcs(kw1,kw0);
		samef1 = alt2    && same_funcs(kw1,kw2);
		if (kw3) {
			if (samef1) {
				if (!needbar)
					Lbrak = kw1->kw.req ? Lpar : Lsqbpar;
				else if (!samef0 && kw3->name)
					Lbrak = Lpar;
				kw3 = 0;
				goto barchk;
				}
			else if (needbar) {
				if (!incc)
					Rbrak = alt2 || kw1->kw.req ? Rpar : Rparsqb;
				if (!samef0
				 && (Stk->n <= nst1
				    || !kw0 || kw0->kw.kw != Stk->Cstak[Stk->n-1]->kw.kw))
					Lbrak = Lpar;
				}
			else if (kw1->kw.req) {
				Lbrak = Lpar;
				Rbrak = Rpar;
				}
			else {
				Lbrak = Lsqb;
				Rbrak = Rsqb;
				if (alt2) {
					Lbrak = Lsqbpar;
					Rbrak = Rpar;
					}
				}
			if (kw2) {
				if (kw4 && !kw4->kw.name && kw4->kw.kw == kw3)
					goto zap_kw3;
				if (Stk->n > 0 && Stk->Cstak[Stk->n-1]->kw.kw == kw3)
					--Stk->n;
				}
			else if (Stk->n > 0 && Stk->Cstak[Stk->n-1]->kw.kw == kw3) {
 zap_kw3:
				Rbrak = 0;
				kw3 = 0;
				}
			}
		else if (!kw1->kw.req) {
			if (!needbar)
				Lbrak = Lsqb;
			if (!alt2)
				Rbrak = Rsqb;
			}
 barchk:
		if (needbar1)
			Bar();
		if (Lbrak)
			Lbrak();
 nobrak:
		Keywd(kw1->kw.name);
		while(++i < j) {
			kw2 = x[i];
			Alias(kw2->kw.name - 6);	/* Alias will increment by 6 */
			}
		knd = kw1->kw.kind;
		if ((k = knd & (KWKind_Mask|KWKind_List)))
			Type(k, typename[k]);
		if ((k = knd & KWKind_Lb)) {
			Saw_tok(knd & KWKind_strictLb ? Saw_gt : Saw_ge, KWKind_Lb);
			btype = k;
			snprintf(buf, sizeof(buf), "%g", kw1->kw.Lb);
			Saw_num(buf);
			}
		if ((k = knd & KWKind_Ub)) {
			Saw_tok(knd & KWKind_strictUb ? Saw_lt : Saw_le, KWKind_Ub);
			btype = k;
			snprintf(buf, sizeof(buf), "%g", kw1->kw.Ub);
			Saw_num(buf);
			}
		if (kw1->init || kw1->cinit) {
			Saw_tok(Saw_dflt, KWKind_init);
			lastdef = 0;
			if (kw1->defname)
				Saw_name(kw1->defname);
			if (kw1->init)
				Saw_num(kw1->init);
			else
				lastkw1->cinit = StrCpy(kw1->cinit);
			}
		if (kw1->funcs)
			lastkw1->funcs = new_StrList(kw1->funcs, 0);
		else if (kw1->kw.f.vf && !kw1->kw.f.final) {
			lastkw1->dylib = StrCpy((const char*)kw1->kw.f.vf);
			dynlibseen = 1;
			}
		else
			lastkw1->funcs = new_StrList("0", 0);
		lastkw1->desc  = StrCpy(kw1->desc);
		lastkw1->group = StrCpy(kw1->group);
		if (kw1->alen) {
			Saw_stok(Saw_len, offsetof(KWctx,len));
			lastkw1->kind |= knd & KWKind_Len1OK;
			Saw_name(kw1->alen);
			}
		if (Rbrak) {
			if (kw3)
				kw_insert((KeyWord*)kw1, Stk, 0);
			Rbrak();
			}
		kw0 = kw1;
		if (top) {
			while(Stk->n > 0) {
				kw_insert(&Stk->Cstak[--Stk->n]->kw, Stk, 0);
				Rpar();
				}
			if (!kw2)
				kw_finish(-1);
			}
		}
	while(Stk->n > nst0) {
		kw_insert(&Stk->Cstak[--Stk->n]->kw, Stk, 0);
		Rpar();
		}
	if (x != x0)
		free(x);
	}

 static void
kw_insert1(KeyWord *kw, int top)
{
	KW_Cstak Stk;

	Stk.n = 0;
	Stk.Cstak = (KeyWordx**)Malloc((Stk.nmax = 64)*sizeof(KeyWordx*), "kw_insert1");
	if (top) {
		KWtop = KWStack;
		KWbase.kw = 0;
		}
	kw_insert(kw, &Stk, top);
	while(Stk.n > 0) {
		kw_insert(&Stk.Cstak[--Stk.n]->kw, &Stk, 0);
		Rpar();
		}
	free(Stk.Cstak);
	}

 static int
libread1(const char *kname, const char *lname, Uint n, int allow_top)
{
	KW_ADD kwa;
	KeyWord *kw, kw0;
	Uint i;
	const char *s;
	void *h, (*Botch)(const char*, ...);

	h = dlopen(lname, RTLD_NOW);
	if (!h) {
		Botch = botch;
		if (kname && libexpand)
			Botch = squawk;
#ifndef _WIN32
		if ((s = dlerror()))
			if (kname)
				Botch("Cannot open library \"%s\" for %s:\n\t%s",
					lname, kname, s);
			else
				Botch("Cannot open library \"%s\":\n\t%s",
					lname, s);
		else
#endif
			if (kname)
				Botch("Cannot open library \"%s\" for %s",
					lname, kname);
			else
				Botch("Cannot open library \"%s\"", lname);
		return ++liberror;
		}
	if (!find_dlsym(kwa, h, "keyword_add"))
		botch("dlsym(\"keyword_add\") failed for %s", lname);
	kw = (*kwa)();
	if (!(s = kw->name)) {
		s = "<NULL>";
		goto namebotch;
		}
	if (allow_top) {
		if (strcmp(s, "KeywordTop")) {
			if (allow_top == 2) {
				memset(&kw0, 0, sizeof(KeyWord));
				kw0.kw = kw;
				kw = &kw0;
				kw0.nkw = 1;
				}
			else
				allow_top = 0;
			}
		else
			kname = 0;
		}
	else if (kname && strcmp(s, kname)) {
 namebotch:
		botch("Library %s: expected top keyword to be %s but got %s",
			lname, kname, s);
		}
	if (kname) {
		i = kw->kind & (KWKind_Mask|KWKind_List);
		if (i != n)
			botch("Library %s: expected %s to be %s, but got %s",
				lname, kname, typename[n], typename[i]);
		}
	if (kw->kw)
		kw_insert1(kw, allow_top);
	dlclose(h);
	return 0;
	}

 static int
lib_read(KWctx *kc, int allow_top)
{
	const char *kname, *lname;

	kname = kc->name1;
	lname = kc->dylib;
	kc->dylib = 0;
	return libread1(kname, lname, kc->kind & (KWKind_Mask|KWKind_List), allow_top);
	}
#else

 static int
lib_read(KWctx *kc, int allow_top)
{
	botch("dlsym(\"keyword_add\") failed for %s", kc->dylib);
	return 1;
	}
#endif

 static void
Keywd(const char *s)
{
	KWctx *kw, *kw0, *kw1, *kw2, *kwm, **kwp;
	StrList *nam;
	Taghash *th;
	char *s1;
	int alt, lsave, req;
	size_t L, L1;
	static char allowed[] = { '(',')','[',']','|', Saw_kwname, Saw_type, 0 };

	if (tagin)
		botch("Inappropriate text \"%s\" in tagfile \"%s\"", s, tagfname);
	kw1 = *KWtop;
#ifndef NO_NIDR_DYNLIB
	if (kw1->dylib && libexpand && !lib_read(kw1,0))
		kw1 = *KWtop;
#endif
	kw0 = kw1;
	kwm = 0;
	if (!kw1->names && s) {
		switch(kw1->parstate) {
		  case PPar_emtpy:
			do {
				kw1->parstate = PPar_starting;
				kw1 = kw1->pparent;
				}
				while(kw1 && kw1->parstate == PPar_emtpy);
			kw1 = kw0->rparent;
			break;
		  case PPar_starting:
			if (or_mode)
				kw1 = kw0->rparent;
			else
				kw1->parstate = PPar_started;
			break;
		  }
		}
	if (or_mode && !kw1->kw) {
		kw2 = KWtop[-1];
		if (!kw2->names && kw2->parstate == PPar_starting)
			kw2 = kw2->rparent;
		if (s)
			kw1->master = kwm = kw1;
		kw1 = kw2;
		}
	L = alt = req = 0;
	if (s) {
		kwp = KWStack + kw1->level + 1;
		if (kwp > KWtop)
			kwp = KWtop;
		while(kwp > KWStack) {
			kw2 = *kwp--;
			if (kw2->name1) {
				dtnext = deftag + kw2->dtnext;
				L = kw2->name1len;
				strcpy(dtnext-L, kw2->name1);
				break;
				}
			}
		alt = ++kw1->nalt;
		if (*Br_top != '['/*]*/ || Br_top == kw1->Br_top) {
			req = ++kw1->nreq;
			if (*Br_top == '(' /*)*/) /* allow [(...)|...] */
				for(s1 = Br_top; --s1 > kw1->Br_top; ) {
					if (*s1 != '(' /*)*/) {
						if (*s1 == '[' /*]*/) {
							--kw1->nreq;
							req = 0;
							}
						break;
						}
					}
			}
		if (or_mode) {
			alt = --kw1->nalt;
			if (req)
				req = --kw1->nreq;
			}
		if (lastseen)
			reqchk(allowed, s);
		or_mode = 0;
		}
	lastkw = lastkw1 = kw = (KWctx*) Alloc(sizeof(KWctx));
	if (!kw0->names)
		kw->pparent = kw0;
	if (s) {
		kw->names = nam = new_StrList(s,0);
		kw->name1 = nam->val;
		kw->name1len = L = strlen(s);
		for(kwp = KWtop; kwp > KWStack; --kwp) {
			kw2 = *kwp;
			if (kw2->names) {
				if (kwp < KWtop) {
					L1 = kw2->dtnext + L + 1;
					while(++kwp <= KWtop) {
						kw2 = *kwp;
						if (!kw2->name1) {
							kw2->name1 = kw->name1;
							kw2->name1len = L;
							kw2->dtnext = L1;
							}
						}
					}
				break;
				}
			}
		if (kw1->dylib) {
			lsave = Lineno;
			Lineno = dynline; /* Line where DYNLIB appeared */
			squawk("Suppressing DYNLIB for %s to parse contained keywords",
				kw1->names ? kw1->names->val : "???? (bug!)");
			Lineno = lsave;
			kw1->dylib = 0;
			}
		nam->primary = KWKind_primary;
		nam->seqno = ++seqno;
		kw->next = kw1->kw;
		kw1->kw = kw;
		if (lastseen == '(' /*)*/)
			kw->parstate = PPar_leftpar;
		kw->agroup = lastagroup++;
		}
	else {
		kw->objno = ++lastobjno;
		if (!kw1->kw && (kw->rparent = kw1->rparent))
			kw1 = kw->rparent;
		else
			kw->rparent = kw1;
		alt = kw1->alt;
		req = kw1->req;
		}
	kw->alt = alt;
	kw->req = req;
	kw->Br_top = Br_top;
	kw->master = kwm;
	kw->level = KWtop - KWStack;
	if (lastseen == '(' || lastseen == '[' /*])*/ || KWtop == KWStack) {
		if (++KWtop - KWStack >= KW_stacklen)
			botch("KW stack overflow");
		*KWtop = kw;
		}
	lasttag = 0;
	if (s) {
		if (dtnext > deftag)
			*dtnext++ = '/';
		if (dtnext + L >= dtend)
			botch("deftag overflow");
		memcpy(dtnext, s, L);
		*(dtnext += L) = 0;
		if (tagout)
			fprintf(tagout, "TAG \"%s\"\n", deftag);
		if (tagfname) {
			if ((th = (Taghash*)DHHlookup(deftag, &TagTab, ignoretags))) {
				++th->uses;
				lasttag = th->hh.name;
				if (!ignoretags) {
					lastkw1->tag = lasttag;
					tagimport(lastkw1, th->kw);
					}
				}
			else if (mtagout)
				fprintf(mtagout, "TAG \"%s\" #missing\n", deftag);
			}
		}
	kw->dtnext = dtnext - deftag;
	if (lastseen == 0)
		kw->kind = kwmult;
	lastseen = Saw_kwname;
	}

 static void
kw_finish(int nextkind)
{
	KWctx *kw;

	if (*Br_top)
		botch("unmatched '%c'", *Br_top);
	if (KWtop - KWStack > 1)
		botch("bug: KW stack botch in kw_finish");
	if (KWtop > KWStack) {
		kw = *KWtop;
		if (kw->dylib && libexpand)
			lib_read(kw, nextkind == -1);
		}
	KWtop = KWStack;
	dtnext = deftag;
	lastseen = 0;
	kwmult = nextkind < 0 ? 0 : nextkind;
	}

 static void
funcbotch(KWctx *kw)
{
	botch("Keyword %s:\n\t\t\"start\" and \"final\" functions cannot be used with DYNLIB"
		".\n\t\tThey must be supplied by the dynamically loaded library",
		kw->names->val);
	}

 static void
Fnc(char *s)
{
	KWctx *kw;

	if (!(kw = lastkw1))
		botch("misplaced {...}");
	if (++kw->nfcn > 4)
		botch("too many functions in {...}");
	else if (kw->dylib)
		funcbotch(kw);
	else if (kw->kind & KWKind_Libname)
		botch("Keyword %s:\n\t\t\"start\" and \"final\" functions "
			"cannot be used with LIBNAME.",
			kw->names->val);
	kw->funcs = wantfuncs ? new_StrList(s, kw->funcs) : 0;
	lastseen = Saw_func;
	}

 static const char*
tname_adj(const char *tname, int what, int kbits)
{
	Uint k = lastkw1->kind;

	switch(what) {
	 case Saw_lt:
		if (k & KWKind_caneqUb)
			tname = "<=";
		break;
	 case Saw_le:
		if (k & KWKind_strictUb)
			tname = "<";
		break;
	 case Saw_ge:
		if (k & KWKind_strictLb)
			tname = ">";
		break;
	 case Saw_gt:
		if (k & KWKind_caneqLb)
			tname = ">=";
	 }
	return tname;
	}

 static void
Saw_tok(int what, int kbits)
{
	const char *tname;

	if ((lastseen == Saw_kwname || lastseen == Saw_type) && lastkw1) {
		if (lastkw1->kind & kbits) {
			tname = tokname[what];
			if (kbits & (KWKind_Lb | KWKind_Ub))
				tname = tname_adj(tname, what, kbits);
			botch("%s already seen for keyword %s",
				tname, lastkw1->names->val);
			}
		lastseen = what;
		}
	else
		unexpected(tokname[what]);
	}

 static void
Saw_stok(int what, size_t woff)
{
	StrList *sl;

	if ((lastseen == Saw_kwname || lastseen == Saw_type) && lastkw1) {
		if (*(void**)((char*)lastkw1 + woff)) {
			if ((sl = lastkw1->names))
				botch("%s already seen for keyword %s",
					tokname[what], sl->val);
			else
				botch("%s already seen", tokname[what]);
			}
		lastseen = what;
		}
	else
		unexpected(tokname[what]);
	}

 static KWctx *
lenfind(const char *s)
{
	KWctx *kw, *kw0, **pkw;
	StrList *sl;

	kw0 = lastkw1;
	for(pkw = KWtop; pkw > KWStack; --pkw) {
		kw = *pkw;
		for(kw = kw->kw; kw; kw = kw->next) {
			for(sl = kw->names; sl; sl = sl->next)
				if (!strcmp(sl->val, s))
					return kw;
			}
		}
	return 0;
	}

 static void
alreadyseen(const char *what)
{
	squawk("%s already given for %s", what, lastkw1->names->val);
	}

 static int
Saw_name(const char *s)
{
	KWctx *kw;
	StrList *sl;
	Taghash *th;
	const char *s1;
	static const char notfound[] =
		"keyword \"%s\" not found in this or an enclosing context";

	if (lastkw1) {
		switch(lastseen) {
		  case Saw_dflt:
			lastseen = Saw_dfltname;
			lastdef = (Defhash*)DHHlookup(s, &DHTab, 1);
			if (!lastkw1->defname)
				lastkw1->defname = def_name(s);
			else if (lastdef->value
				 && ((s1 = lastkw1->init) || (s1 = lastkw1->cinit))
				 && strcmp(s1, lastdef->value))
					squawk("Inconsistent definitions for %s: %s and %s\n",
						lastdef->hh.name, lastdef->value, s1);
			return 0;
		  case Saw_len:
			if (lastkw1->len) {
				alreadyseen("LEN");
				goto ret1;
				}
			kw = lenfind(s);
			if (!kw)
				squawk(notfound,s);
			else {
				if ((kw->kind & KWKind_Mask) != KWKind_Int)
					squawk("%s is not of type INTEGER", kw->names->val);
				for(sl = kw->names; sl->next; sl = sl->next);
				lastkw1->len = sl->val;
				}
			goto ret1;
		  case Saw_group:
			if (lastkw1->group)
				alreadyseen("GROUP");
			else
				strcpy(lastkw1->group = (char*)Alloc(strlen(s)+1), s);
			goto ret1;
		  case Saw_tag:
			if (tagin)
				goto new_tag;
			if (lasttag) {
				alreadyseen("TAG");
				goto ret1;
				}
			if (!tagfname && !ignoretags) {
				squawk("TAG \"%s\" not found: no tag file given", s);
				goto ret1;
				}
			th = (Taghash*)DHHlookup(s, &TagTab, ignoretags);
			if (!th)
				squawk("TAG \"%s\" not found in tag file \"%s\"", s, tagfname);
			else {
				++th->uses;
				lasttag = th->hh.name;
				if (!ignoretags) {
					lastkw1->tag = lasttag;
					tagimport(lastkw1, th->kw);
					}
				}
 ret1:
			lastseen = Saw_type;
			return 1;
		  }
		}
	else if (tagin) {
 new_tag:
		th = (Taghash*)DHHlookup(s, &TagTab, 1);
		if (th->kw)
			squawk("Tag \"%s\" already seen in tagfile %s\n", s, tagfname);
		else
			th->kw = lastkw1 = (KWctx*) Alloc(sizeof(KWctx));
		goto ret1;
		}
	botch("unexpected name \"%s\"", s);
	return 0; /* not reached */
	}

 static void
boundalready(const char *what)
{
	squawk("%s bound already given for %s", what, lastkw1->names->val);
	}

 static void
set_lastdef(const char *s, int qkeep)
{
	if (!lastdef->value) {
		lastdef->value = def_name(s);
		lastdef->qkeep = qkeep;
		}
	else if (strcmp(s, lastdef->value))
		squawk("Two := values for %s: %s and %s", lastdef->hh.name,
			lastdef->value, s);
	}

 static void
Saw_num(const char *s)
{
	Defhash *ld;
	double x, xa;
	char *se, **sp;
	int k;
	static double b[3];

	static char allowed[] = { Saw_dflt, Saw_dfltname, Saw_lt, Saw_le, Saw_ge, Saw_gt, 0 };

	reqchk(allowed, "number");
	x = strtod(s,&se);
	if (*se)
		botch("invalid number: \"%s\"", s);
	ld = 0;
	sp = 0;
	switch(lastseen) {
		case Saw_dflt:
		case Saw_dfltname:
			sp = &lastkw1->init;
			lastkw1->kind |= KWKind_init;
			ld = lastdef;
			if ((lastkw1->kind & KWKind_Mask) == KWKind_Int)
				x = x >= 0. ? floor(x) : ceil(x);
			b[2] = x;
			if (lastkw1->lb && x < b[0])
				squawk("default value inconsistent with lower bound");
			if (lastkw1->ub && x > b[1])
				squawk("default value inconsistent with upper bound");
			break;
		case Saw_gt:
		case Saw_ge:
			if (lastkw1->lb) {
				boundalready("lower");
				goto ret;
				}
			if ((lastkw1->kind & KWKind_Mask) == KWKind_Int) {
				xa = ceil(x);
				if (xa > x) {
					if (btype == KWKind_strictLb)
						btype = KWKind_caneqLb;
					x = xa;
					}
				}
			sp = &lastkw1->lb;
			lastkw1->kind |= btype;
			b[0] = x;
			if (lastkw1->ub && x >= b[1])
				squawk("inconsistent bounds");
			if (lastkw1->init && x > b[2])
				squawk("lower bound inconsistent with default value");
			break;
		case Saw_le:
		case Saw_lt:
			if (lastkw1->ub) {
				boundalready("upper");
				goto ret;
				}
			if ((lastkw1->kind & KWKind_Mask) == KWKind_Int) {
				xa = floor(x);
				if (xa < x) {
					if (btype == KWKind_strictUb)
						btype = KWKind_caneqUb;
					x = xa;
					}
				}
			sp = &lastkw1->ub;
			lastkw1->kind |= btype;
			b[1] = x;
			if (lastkw1->lb && x <= b[0])
				squawk("inconsistent bounds");
			if (lastkw1->init && x < b[2])
				squawk("upper bound inconsistent with default value");
			break;
		default:
			unexpected(s);
		}
	k = lastkw1->kind & KWKind_Mask;
	if ((k == KWKind_Void && !tagin) || k == KWKind_Str)
		squawk("cannot accept a %s value", tokname[lastseen]);
	else {
		strcpy(*sp = (char*)Alloc(strlen(s)+1), s);
		if (ld)
			set_lastdef(s, 0);
		}
 ret:
	lastseen = Saw_type;
	}

%}

ALias	"ALIAS"[ \t]+[a-z][a-z_0-9]*
D	[0-9]
Kword	[a-z][a-z_0-9]*
Func	[a-zA-Z_0-9()][,a-zA-Z_0-9<>():.+*/\-\[\] \t]*
Num	[\-+]?({D}+("."{D}*)?|"."{D}+)([eE]([\-+]?){D}+)?
PPname	[a-zA-Z][a-zA-Z0-9_]*
ws	[ \t]+

%s FNC VAL

%%

\n		{ ++Lineno; }
{ws}		{}
<INITIAL>{ALias}	{ Alias(yytext); }
<INITIAL>"INTEGER"	{ Type(KWKind_Int, yytext); }
<INITIAL>"INTEGERLIST"	{ Type(KWKind_Int | KWKind_List, yytext); }
<INITIAL>"KEYWORD"	{ kw_finish(0); }
<INITIAL>"KEYWORD1"	{ kw_finish(KWKind_1); }
<INITIAL>"KEYWORD01"	{ kw_finish(KWKind_01); }
<INITIAL>"KEYWORD12"	{ kw_finish(KWKind_12); }
<INITIAL>"LIBNAME"	{ Type(KWKind_Str|KWKind_Libname, yytext); }
<INITIAL>"REAL"		{ Type(KWKind_Real, yytext); }
<INITIAL>"REALLIST"	{ Type(KWKind_Real | KWKind_List, yytext); }
<INITIAL>"STRING"	{ Type(KWKind_Str, yytext); }
<INITIAL>"STRINGLIST"	{ Type(KWKind_Str | KWKind_List, yytext); }
<INITIAL>{Kword}	{ Keywd(yytext); }
<INITIAL>";"		{ kw_finish(-1); }
<INITIAL>"DESC"		{ Saw_stok(Saw_desc, offsetof(KWctx,desc));	BEGIN VAL; }
<INITIAL>"TAG"		{ Saw_stok(Saw_tag, offsetof(KWctx,tag));	BEGIN VAL; }
<INITIAL>"GROUP"	{ Saw_stok(Saw_group, offsetof(KWctx,group));	BEGIN VAL; }
<INITIAL>"DYNLIB"	{ Saw_tok(Saw_dy, offsetof(KWctx,dylib));
			  dynlibseen = 1; dynline = Lineno;		BEGIN VAL; }
<INITIAL>"LEN1"		{ Saw_stok(Saw_len, offsetof(KWctx,len));
			  lastkw1->kind |= KWKind_Len1OK;		BEGIN VAL; }
<INITIAL>"LEN"		{ Saw_stok(Saw_len, offsetof(KWctx,len));	BEGIN VAL; }
<INITIAL>":="		{ Saw_tok(Saw_dflt, KWKind_init); lastdef = 0;	BEGIN VAL; }
<INITIAL>">="		{ Saw_tok(Saw_ge, KWKind_Lb); btype = KWKind_caneqLb;	BEGIN VAL; }
<INITIAL>">"		{ Saw_tok(Saw_gt, KWKind_Lb); btype = KWKind_strictLb;	BEGIN VAL; }
<INITIAL>"<="		{ Saw_tok(Saw_le, KWKind_Ub); btype = KWKind_caneqUb;	BEGIN VAL; }
<INITIAL>"<"		{ Saw_tok(Saw_lt, KWKind_Ub); btype = KWKind_strictUb;	BEGIN VAL; }
<VAL>{PPname}		{ if (Saw_name(yytext))	BEGIN INITIAL; }
<VAL>{Num}		{ Saw_num(yytext);	BEGIN INITIAL; }
<VAL>\"			{ Saw_quote();		BEGIN INITIAL; }

"{"		{ Lbrace(); BEGIN FNC; }
"}"		{ Rbrace(); BEGIN INITIAL; }
"|"		{ Bar(); }
"("		{ Lpar(); }
")"		{ Rpar(); }
"["		{ Lsqb(); }
"]"		{ Rsqb(); }
<FNC>{Func}	{ Fnc(yytext); }
<FNC>";"	{}
"#".*\n		{++Lineno;}
"/*"		{ Eat_comment(); }
.	{ botch("\"%s\" unexpected", yytext); }

%%

 static void
cantopen(char *what)
{ fprintf(stderr, "%s: cannot open \"%s\"\n", progname, what); }

 int
yywrap(void)
{
	if (tagin) {
		if (tagin != stdin)
			fclose(tagin);
		tagin = 0;
		Lineno = 1;
		saweof = 1;
		if (*specfile == '-' && !specfile[1]) {
			yyin = stdin;
			infname = 0;
			return 0;
			}
		if ((yyin = fopen(infname = specfile, "r")))
			return 0;
		cantopen(specfile);
		exit(1);
		}
	return 1;
	}

 static void
Eat_comment(void)
{
	int c;

	for(;;) {
		c = input();
 c_check:
		if (c == EOF)
			botch("EOF within /* comment */");
		if (c == '*') {
			c = input();
			if (c == '/')
				return;
			goto c_check;
			}
		if (c == '\n')
			++Lineno;
		}
	}

 static void
Saw_quote(void)
{
	KWctx *kw, *kw1, *kw2;
	char buf[1024];
	char *b, *b0, *be, *lname;
	int c, c0;
	size_t L, L0, L1;
	unsigned char *us;

	static char escname[32];
	static char hexdig[16] = "0123456789abcdef";
	static unsigned char esetup[] = "\aa\bb\ff\nn\rr\tt\vv";

	if (!escname['\a'])
		for(us = esetup; *us; us += 2)
			escname[*us] = us[1];	/* assuming ASCII */
	b0 = b = buf;
	L = sizeof(buf);
	be = b + L - 3;	/* allow room for /xhh */
	c0 = 0;
	saweof = 0;
	for(;;) {
		c = input();
		if (c == EOF)
			botch("EOF within quoted string");
		if (b >= be) {
			L1 = b - b0;
			L0 = L;
			L <<= 1;
			be = (char*)Malloc(L, "Saw_quote");
			memcpy(be, b0, L1);
			if (b0 != buf)
				free(b0);
			b0 = be;
			b = b0 + L1;
			be = b0 + L - 3;
			}
		if (c < ' ') {
			*b++ = '\\';
			if (escname[c]) {
				if (c == '\n')
					++Lineno;
				c = escname[c];
				}
			else {
				*b++ = 'x';
				*b++ = '0' + (c >> 4);
				c = hexdig[c & 0xf];
				}
			}
		*b++ = c;
		if (c0 == '\\')
			c0 = 0;
		else {
			c0 = c;
			if (c == '"') {
				while((c = input()) <= ' ') {
					if (c == EOF || !c) {
						c = ';';
						break;
						}
					if (c == '\n')
						++Lineno;
					}
				if (c == '"' && !saweof) {
					--b;
					continue;
					}
				unput(c);
				break;
				}
			}
		}
	b[-1] = 0;	/* replacing " */
	if (lastseen == Saw_tag) {
		if (saweof) /* kludge around poor flex design */ {
			tagin = stdin;
			Saw_name(b0);
			tagin = 0;
			}
		else
			Saw_name(b0);
		}
	else {
		L = b - b0;
		be = (char*)Alloc(L);
		memcpy(be, b0, L);
		switch(lastseen) {
		  case Saw_dflt:
		  case Saw_dfltname:
			if ((lastkw1->kind & KWKind_Mask) != KWKind_Str)
				botch("%s cannot have a string-valued := value",
					lastkw1->names->val);
			lastkw1->cinit = be;
			if (lastdef)
				set_lastdef(b0, 1);
			break;
		  case Saw_desc:
			lastkw1->desc = be;
			break;
		  case Saw_dy:
			lastkw1->dylib = be;
			if (lastkw1->nfcn > 0)
				funcbotch(lastkw1);
			if (libexpand) {
				kw = *KWtop;
				lastseen = Saw_kwname;
				if (kw == lastkw1)
					kw = *--KWtop;
				else if ((kw2 = kw->rparent))
					kw = kw2;
				kw1 = lastkw1;
				if ((kw2 = kw->kw) == kw1) {
					if ((kw2 = kw->kw = kw2->next)
					 && kw2->alt == kw1->alt)
						Bar();
					else  {
						--kw->alt;
						if (kw1->req)
							--kw->nreq;
						}
					}
				Lpar();
				Keywd(kw1->name1);
				lname = kw1->dylib;
				if (lib_read(kw1,0))
					lastkw1->dylib = lname;
				Rpar();
				if (kw1->parstate == PPar_leftpar)
					++KWtop;
				}
			break;
		  case Saw_group:
			lastkw1->group = be;
			break;
		  default:
			unexpected("string");
		  }
		}
	if (b0 != buf)
		free(b0);
	lastseen = Saw_type;
	}

 static int
usage(int rc)
{
	fprintf(rc ? stderr : stdout,
	"Usage: %s [options] specfile [tagfile [defs_out [keywds_out]]]\n\
	where the last two are output files, \".\" means \"omit\", and \"-\"\n\
	means write to stdout.  Missing trailing arguments are not read\n\
	or written, except that if all outputs are missing, \". -\" is assumed.\n\
	Options may modify this behavior.\noptions:\n\
	-D		{ just write definitions to stdout }\n\
	-d		{ suppress test for ambiguity from duplicate keywords }\n\
	-e		{ expand ((...)) for debugging }\n\
	-f		{ omit functions, i.e., {...} with -p; without -p,\n\
			  supply NULL function pointers }\n\
	-G		{ for keywords with a TAG but not GROUP, do not assume\n\
			  that implictly GROUP == TAG }\n\
	-g		{ write GuiKeyWord rather than KeyWord header;\n\
			  with -p, include GUI details: TAG, DESC, GROUP, LEN }\n\
	-h headername	{ Specify header file(s), e.g., \"-h nidr.h\";\n\
			  several -h options can be given;\n\
			  default = -h- ==> no headers }\n\
	-j specsum	{ write file specsum for use by \"dakreord specsum ...\"\n\
			  and \"dakreorder -j specsum ...\" }\n\
	-k keyname	{ library mode with input containing several keywords\n\
			  to be contained in keyword keyname; implies -l and\n\
			  has no other effect if -D, -g, -j, -p, or -s is given.\n\
			  If -l is given without -k and several keywords appear\n\
			  in the input, assume \"-k KeywordTop\". }\n\
	-L libname	{ obtain keyword details from library libname (compiled\n\
			  from source produced by an earlier \"nidrgen -l ...\"\n\
			  invocation); non-option arguments in this case are\n\
			  [defs_out [keywds_out]] }\n\
	-l		{ library mode: with -D, -g, -j, -p, or -s, attempt to\n\
			  load libraries specified with \"DYNLIB libname\" and\n\
			  incorporate keyword specifications provided by the\n\
			  loaded libraries; otherwise generate source for\n\
			  compilation into a shared library for loading when a\n\
			  keyword marked with \"DYNLIB libname\" (in a separate\n\
			  run of nidrgen) is seen.  The library provides\n\
			  contained keywords and any needed final routine for\n\
			  the marked keyword. }\n\
	-m mtagfile	{ just write to mtagfile tags not found in tagfile }\n\
	-n namespname	{ assume namespace namespname;\n\
			  default = -n Dakota\n\
			  -n- ==> no namespace }\n\
	-p[w]		{ pretty-print rather than produce KeyWord header;\n\
			  add w spaces (default 2) for each nested context. }\n\
	-s[w]		{ like -p[w], but sort the keywords }\n\
	-T outtagfile	{ write dummy tagfile to outtagfile (and quit) }\n\
	-t		{ ignore TAG fields if tagfile is not given }\n\
	-w		{ cancel warn mode: do not report (on stderr)\n\
			  the number of implicitly named start routines;\n\
			  still flag their names  with /*!!*/ }\n",
		progname);
	return rc;
	}

 static StrList *
reverse(StrList **x)
{
	StrList *s, *t, *u;

	if (!(s = *x) || !s->next)
		return s;
	for(t = 0; s; s = u) {
		u = s->next;
		s->next = t;
		t = s;
		}
	return *x = t;
	}

 static int
kwcomp(const void *a, const void *b)
{
	return strcmp((*(const KWctx**)a)->names->val, (*(const KWctx**)b)->names->val);
	}

 static KWctx *
kwsort(KWctx *kw0)
{
	/* Sort keywords after expanding kw list with aliases; compute nkw, */
	/* including aliases. */

	KWctx *kw, *kw1, *kw2, *kwm, **x, **x0, **xe, *x00[512];
	StrList *dname, *minname, *nam, **pmin, **pnam;
	char *descsave, *groupsave;
	int agadj, nk;
	size_t L, nkw, nkwa, onkwa;

	reverse(&kw0->funcs);
	nkw = nkwa = 0;
	for(kw = kw0->kw; kw; kw = kw->next) {
		++nkw;
		pmin = &kw->names;
		minname = *pmin;
		pnam = &minname->next;
		onkwa = nkwa;
		while((nam = *pnam)) {
			dname = nam;
			++nkwa;
			if (strcmp(nam->val, minname->val) < 0) {
				pmin = pnam;
				minname = nam;
				}
			pnam = &nam->next;
			}
		agadj = 2*kw->agroup;
		kw->agroup = agadj + 1;
		if (nkwa > onkwa) {
			*pmin = minname->next;
			minname->next = kw->names;
			nam = kw->names = minname;
			nk = kw->kno;	/* ALias() made room for nk use below */
			descsave = kw->desc;
			groupsave = kw->group;
			if (dname != minname) {
				kw->desc = kw->group = 0;
				kw->agroup = agadj;
				}
			kwm = kw;
			kw = kw2 = 0;
			if (nam->primary)
				kw2 = kwm;
			while((nam = nam->next)) {
				kw1 = (KWctx*)Alloc(sizeof(KWctx));
				if (!kw)
					kw = kw1;	/* so kw->next is right */
				kw1->next = kwm->next;
				kwm->next = kw1;
				kw1->names = nam;
				kw1->master = kwm;
				kw1->kind = kwm->kind;
				kw1->kno = ++nk;
				kw1->agroup = agadj;
				if (nam == dname) {
					kw1->desc = descsave;
					kw1->group = groupsave;
					kw1->agroup = agadj + 1;
					}
				if (nam->primary)
					kw2 = kw1;
				}
			if (kw2 != kwm && !guikeywds) {
				kw2->defname = kwm->defname;
				kwm->defname = 0;
				kw2->len = kwm->defname;
				kwm->defname = 0;
				kw2->len = kwm->len;
				kwm->len = 0;
				kw2->init = kwm->init;
				kwm->init = 0;
				kw2->cinit = kwm->cinit;
				kwm->cinit = 0;
				kw2->dylib = kwm->dylib;
				/* no need to change kwm->dylib */
				}
			for(kw1 = kwm;; kw1 = kw1->next) {
				if (kw1 != kw2)
					kw1->pmaster = kw2;
				if (kw1 == kw)
					break;
				}
			}
		}
	kw0->nkw = (int)(nkw += nkwa);
	if (nkw <= 1)
		return kw0->kw;
	x0 = x00;
	if (nkw > sizeof(x00)/sizeof(x00[0]))
		x0 = (KWctx**)Malloc(L = nkw*sizeof(KWctx*), "kwsort");
	for(x = x0, kw = kw0->kw; kw; kw = kw->next)
		*x++ = kw;
	qsort(x0, nkw, sizeof(KWctx*), kwcomp);
	for(xe = x - 1, x = x0; x < xe; ++x) {
		if (!strcmp(x[0]->names->val, x[1]->names->val)) {
			fprintf(stderr, "%s: duplicate appearance of \"%s\"\n",
				progname, x[0]->names->val);
			++dupnames;
			}
		x[0]->next = x[1];
		}
	x[0]->next = 0;
	for(x = x0, ++xe; x < xe; ++x)
		(*x)->paoff = x - x0;
	for(x = x0; x < xe; ++x)
		if ((kw1 = (kw = *x)->pmaster))
			kw->paoff = kw1->paoff - kw->paoff;
	for(x = x0; x < xe; ++x)
		if (!(kw = *x)->pmaster)
			kw->paoff = 0;
	kw0->kw = kw = x0[0];
	if (x0 != x00)
		free(x0);
	return kw;
	}

 static int
alphanum(char *s)
{
	for(; *s; ++s)
		if (!isalnum(*s))
			return 0;
	return 1;
	}

 enum CommaState {
	CS_Comma = 0,
	CS_Clbra = 1,
	CS_rbraC = 2,
	CS_lrbrC = 3,
	CS_rbr2C = 4
	};

 static const char
	*Commasx[5] = { ",", ",{", "},", ",{},",/*{{*/ "}}," },
	*Commas0[5] = { ",", ",", ",", ",", "," };
 static const char **Comma;

 typedef struct Zbuf Zbuf;
 struct
Zbuf {
	char buf[80];	/* overkill */
	char *b;
	int nnz;
	enum CommaState cs;
	};

 static void
iput(int i, Zbuf *z)
{
	char *s;
	const char *t;

	t = Comma[z->cs];
	if (i == 0) {
		s = z->b;
		for(; (*s = *t); ++s, ++t);
		*s++ = '0';
		z->b = s;
		}
	else {
		if (z->b > z->buf) {
			*z->b = 0;
			printf(z->b = z->buf);
			}
		printf("%s%d", t, i);
		++z->nnz;
		}
	z->cs = CS_Comma;
	}

 static void
nput(char *s, Zbuf *z)
{
	char *fmt, *s1;
	const char *t;

	t = Comma[z->cs];
	if (!s || (strtod(s,&s1) == 0. && !*s1)) {
		s = z->b;
		for(; (*s = *t); ++s, ++t);
		*s++ = '0';
		*s++ = '.';
		z->b = s;
		}
	else {
		fmt = "%s%s.";
		for(s1 = s; *s1; ++s1)
			switch(*s1) {
				case '.':
				case 'e':
				case 'E':
					fmt = "%s%s";
					goto have_fmt;
				}
 have_fmt:
		if (z->b > z->buf) {
			*z->b = 0;
			printf(z->b = z->buf);
			}
		printf(fmt, t, s);
		++z->nnz;
		}
	z->cs = CS_Comma;
	}

 static void
sput(char *s, Zbuf *z)
{
	const char *t;

	t = Comma[z->cs];
	if (!s) {
		s = z->b;
		for(; (*s = *t); ++s, ++t);
		*s++ = '0';
		z->b = s;
		}
	else {
		if (z->b > z->buf) {
			*z->b = 0;
			printf(z->b = z->buf);
			}
		printf("%s\"%s\"", t, s);
		++z->nnz;
		}
	z->cs = CS_Comma;
	}

 static void
zbclear(Zbuf *z)
{
	if (z->b > z->buf) {
		*z->b = 0;
		printf("%s", z->b = z->buf);
		}
	z->cs = CS_Comma;
	}

 typedef struct
FCat { char *b, *b0, *be; } FCat;

 static char *
func_cat(FCat *fc, StrList *f)
{
	char *b, *b1, *be;
	int k, n;
	size_t L, L1, L2, L3;

	if (!f)
		return 0;
	if (!f->next)
		return f->val;
	b = fc->b;
	be = fc->be;
	for(n = k = 0; f; f = f->next, ++n) {
		L = L1 = strlen(f->val);
		if (n) {
			++L1;
			if (n & 1) {
				if (alphanum(f->val)) {
					k = 2;
					L1 += 7;
					}
				else {
					k = 3;
					L1 += 9;
					}
				}
			else
				k = 1;
			}
		if (L1 > be - b) {
			L2 = b - fc->b;
			L3 = 2*(be - fc->b + L1);
			b1 = (char*)Malloc(L3, "func_cat");
			if (L2)
				memcpy(b1, fc->b, L2);
			if (fc->b != fc->b0)
				free(fc->b);
			fc->b = b1;
			be = fc->be = b1 + L3;
			b = b1 + L2;
			}
		switch(k) {
			case 0: break;
			case 1: *b++ = ',';
				break;
			case 2: memcpy(b, ",(void*)", 8);
				b += 8;
				break;
			case 3: memcpy(b, ",(void*)(", 9);
				b += 9;
			}
		strcpy(b, f->val);
		b += L;
		if (k == 3)
			*b++ = ')';
		*b = 0;
		}
	return fc->b;
	}

 static void
kwout(KWctx *kw, KWctx **kwtop, Uint xflag)
{
	FCat fc;
	KWctx *kw0, *kw1, *kw2, *kwm, *kwp, **kws, **x, *x0[KWOUTlen0];
	StrList *S, *nam;
	Zbuf z;
	char fbuf[256], *funcs, *s;
	const char *fmt, *ffmt;
	int aoff, k, i, m, m1, n, roff, sentinnel, skipfuncs;

	if (expand && !kw->names)
		return;
	if (needcomma) {
		printf(",\n");
		needcomma = 0;
		}
	kw1 = kw->kw;
	n = 0;
	fc.b = fc.b0 = fbuf;
	fc.be = fbuf + sizeof(fbuf);
	if (kw->names) {
		kw0 = kw1;
		if (kw->parstate != KwStacked)
			for(kw2 = kw->pparent; kw2; kw2 = kw2->pparent) {
				kw2->next = kw1;
				kw1 = kw2;
				++n;
				}
		if (!kw0) {
			if (!--n) {
				kw->nkw = kw1->nkw;
				if (kw1->kwkno) {
					kw->kwkno = kw1->kwkno;
					return;
					}
				kw->kwkno = lastkno + 1;
				kw = kw1;
				kw1 = kw->kw;
				}
			else {
				kw2 = kw->pparent;
				if (kw2->kwknop) {
					kw->kwkno = kw2->kwknop;
					return;
					}
				kw2->kwknop = lastkno + 1;
				kw1 = kw2;
				while((kw2 = kw2->pparent)) {
					kw2->next = kw1;
					kw1 = kw2;
					}
				}
			}
		}
	for(m = 0, kw0 = kw1; kw0; kw0 = kw0->next)
		++m;
	if (expand)
		for(kw0 = kw1; kw0 && !kw0->names; kw0 = kw0->next)
			for(kw2 = kw0->kw; kw2; kw2 = kw2->next)
				++m;
	x = x0;
	if (m > KWOUTlen0)
		x = (KWctx**)Malloc(m*sizeof(KWctx*), "kwout");
	m = 0;
	kw0 = kw1;
	if (expand) {
		aoff = kw->nalt;
		roff = kw->nreq;
		for(; kw0 && !kw0->names; kw0 = kw0->next) {
			for(kw2 = kw0->kw; kw2; kw2 = kw2->next) {
				kw2->altoffset = aoff;
				kw2->reqoffset = roff;
				x[m++] = kw2;
				}
			aoff += kw0->nalt;
			roff += kw0->nreq;
			}
		kw->nalt = aoff;
		kw->nreq = roff;
		}
	for(; kw0; kw0 = kw0->next)
		x[m++] = kw0;
	sentinnel = 0;
	if (expand) {
		kw->nkw = m;
		qsort(x, m, sizeof(KWctx*), kwcomp);
		for(i = 1; i < m; ++i)
			if (!strcmp(x[i-1]->names->val, x[i]->names->val)) {
				fprintf(stderr, "%s: duplicate appearance of \"%s\"\n",
					progname, x[i]->names->val);
				++dupnames;
				}
		}
	else if (kw->nkw == 0)
		sentinnel = 1;
	m1 = m + sentinnel;
	printf("\tkw_%d[%d] = {\n" /*}*/, kw->kwkno = ++lastkno, m1);
	needcomma = 1;
	z.nnz = 0;
	skipfuncs = !wantfuncs && !xflag;
	for(i = 0; i < m;) {
		kw1 = x[i++];
		if (!(kwm = kw1->master))
			kwm = kw1;
		if (!(kwp = kw1->pmaster))
			kwp = kw1;
		fputs(xflag ? "\t\t{{" : "\t\t{" /*}}}*/, stdout);
		z.b = z.buf;
		z.cs = CS_Comma;
		if ((nam = kw1->names)) {
			k = kwp->kind | kw1->names->primary | xflag;
			if (kwp->dylib)
				k |= KWKind_Dynlib;
			printf(k < 16
				? "\"%s\",%d"
				: "\"%s\",0x%x", kw1->names->val, k);
			iput(kwm->nkw, &z);
			iput(kwm->alt + kwm->altoffset, &z);
			iput(kwm->req ? kwm->req + kwm->reqoffset : 0, &z);
			}
		else {
			if (xflag)
				printf("0,0x%x", xflag);
			else
				printf("0,0");
			iput(kwm->nkw, &z);
			iput(0, &z);
			iput(0, &z);
			}
		if (guikeywds)
			iput(kw1->agroup, &z);
		if ((n = kwm->kwkno)) {
			zbclear(&z);
			printf(xflag ? ",(KeyWord*)kw_%d": ",kw_%d", n);
			}
		else
			iput(0, &z);
		nput(kwp->lb, &z);
		nput(kwp->ub, &z);
		ffmt = /*{{{*/ !xflag ? "}%s" : "}}%s";
		if (guikeywds) {
			nput(kwm->init, &z);
			sput(kwm->cinit, &z);
			sput(kw1->desc, &z);
			if (!(s = kw1->group) && dfltgroup)
				s = kw1->tag;
			sput(s, &z);
			if (kwm->len)
				sput(kwm->len, &z);
			if (kwm->dylib)
				sput(kwm->dylib, &z);
			goto nofuncs;
			}
		iput(kw1->paoff, &z);
		if (skipfuncs && !kwp->dylib)
			goto nofuncs;
		z.nnz = 0;
		funcs = 0;
		if (kwp->dylib) {
			z.cs = CS_Clbra;
			for(n = 0; n < 3; ++n)
				iput(0,&z);
			zbclear(&z);
			printf(",(void*)\"%s\"", kwp->dylib);
			z.cs = CS_rbr2C;
			z.nnz = 0;
			}
		else if ((S = kwp->funcs)) {
			if (!S->next && !strcmp(S->val,"0"))
				goto xcheck;
			if (!kwp->kw && S->needrev)
				S = reverse(&kwp->funcs);
			funcs = func_cat(&fc, S);
			goto funcput;
			}
		else if ((S = kwm->funcs)) {
			funcs = func_cat(&fc, S);
 funcput:
			z.cs = CS_lrbrC;
			if (funcs && wantfuncs & 1) {
				zbclear(&z);
				printf(",%s", funcs);
				z.cs = CS_rbraC;
				}
			else
				iput(0,&z);
			}
		else if (kwm->names && wantfuncs & 1 && !(kwm->kind & KWKind_Libname)) {
			++missing;
			printf(",/*!!*/");
			for(kws = KWStack; ++kws <= kwtop; ) {
				if (!(kw2 = *kws)->names)
					kw2 = kw2->rparent;
				printf("%s_", kw2->names->val);
				}
			printf("%s_start", kwm->names->val);
			}
 xcheck:
		if (xflag) {
			if (z.nnz) {
				z.cs = z.cs == CS_rbraC ? CS_rbr2C : CS_rbraC;
				ffmt = /*{{{*/ "}}}%s";
				z.nnz = 0;
				}
			else if (z.cs == CS_rbr2C)
				ffmt = /*{{{*/ "}}}%s";
			else {
				z.cs = CS_rbraC;
				ffmt = /*{{*/ "}}%s";
				}
			z.b = z.buf;
			iput(nam ? nam->seqno : 0, &z);
			iput(kw1->agroup, &z);
			sput(funcs, &z);
			sput(kw1->desc, &z);
			sput(kw1->group, &z);
			sput(kw1->defname, &z);
			sput(kw1->len, &z);
			sput(kw1->init, &z);
			sput(kw1->cinit, &z);
			if (z.nnz)
				ffmt = /*{*/ "}%s";
			}
 nofuncs:
		printf(ffmt, i < m1 ? ",\n" : "\n");
		}
	if (sentinnel)
		fputs("\t\t{\"\"}\n", stdout);
	printf(/*{*/"\t\t}");
	if (x != x0)
		free(x);
	if (fc.b != fbuf)
		free(fc.b);
	}

 static void
hashclear(KWctx *kw)
{
	StrList *sl;

	kw->parstate = Dupchecked;
	for(kw = kw->kw; kw; kw = kw->next)
		for(sl = kw->names; sl; sl = sl->next) {
			sl->h.hnext->hprev = sl->h.hprev;
			sl->h.hprev->hnext = sl->h.hnext;
			sl->h.hprev = sl->h.hnext = 0;
			}
	}

 static void
dupfound(HashHead *h0, StrList *sl, StrList *sl1)
{
	StrList *sl2, *sl3;

	++dupnames;
	for(sl2 = sl1; sl2->h.hnext != h0; sl2 = sl3) {
		sl3 = (StrList*)sl2->h.hnext;
		if (strcmp(sl->val, sl3->val))
			break;
		}
	fprintf(stderr, "Warning: \"%s\" on line %d also appears on line%s %d",
		sl->val, sl->Lineno, "s" + (sl2 == sl1), sl1->Lineno);
	while(sl1 != sl2) {
		sl1 = (StrList*)sl1->h.hnext;
		fprintf(stderr, ", %d", sl1->Lineno);
		}
	fputs(".\n", stderr);
	}

 static KWctx **
hashadd(KWctx *kw, KWctx **kwtop)
{
	HashHead *h, *h0;
	KWctx *kw1;
	StrList *sl, *sl1;
	char *s;
	int n;
	unsigned int c, x;

	kw->parstate = Dupchecking;
	for(kw1 = kw->kw; kw1; kw1 = kw1->next)
		for(sl = kw1->names; sl; sl = sl->next) {
			x = 0;
			s = sl->val;
			while((c = *(unsigned char*)s++))
				x = 43*x + c - ' ';
			h0 = &HTab[x % HTlen];
			for(h = h0->hnext; h != h0; h = h->hnext) {
				sl1 = (StrList*)h;
				n = strcmp(sl->val, sl1->val);
				if (n <= 0) {
					if (n == 0)
						dupfound(h0,sl,sl1);
					break;
					}
				}
			sl->h.hnext = h;
			h->hprev = (sl->h.hprev = h->hprev)->hnext = (HashHead *)sl;
			}
	if ((kw1 = kw->kw))
		*++kwtop = kw1;
	return kwtop;
	}

 static void
dup_check(void)
{
	HashHead *h, *he;
	KWctx *kw, *kw1, *kw2, **kwtop, **kwbot, *plist;
	int lev;

	for(h = HTab, he = h + HTlen; h < he; ++h)
		h->hnext = h->hprev = h;
	plist = 0;
	kwtop = kwbot = KWStack;
	*kwtop = &KWbase;
	for(;;) {
		kw = *kwtop;
		lev = kw->level;
		while(plist && plist->level >= lev) {
			hashclear(plist);
			plist = plist->next;
			}
		if (kw->parstate == Dupchecking) {
			hashclear(kw);
			if ((*kwtop = kw->next))
				continue;
			if (kwtop <= kwbot)
				break;
			--kwtop;
			continue;
			}
		if ((kw1 = kw->pparent) && kw1->parstate != Dupchecking) {
			for(kw2 = kw1; plist && kw2; plist = plist->next) {
				while(kw2->level > plist->level && (kw2 = kw2->pparent));
				if (kw2 != plist)
					break;
				hashclear(plist);
				}
			kw2 = 0;
			do {
				if (kw1->kw) {
					kw1->next = kw2;
					kw2 = kw1;
					}
				}
				while((kw1 = kw1->pparent) && kw1->parstate != Dupchecking);
			while(kw2) {
				kwtop = hashadd(kw2, kwtop);
				kw1 = kw2;
				kw2 = kw2->next;
				kw1->next = plist;
				plist = kw1;
				}
			}
		kwtop = hashadd(kw, kwtop);
		}
	}

 static void
prkwname(KWctx *kw, const char *suf, int needrev)
{
	int k;
	StrList *S;
	static const char *kname[] = { "BUG", "INTEGER", "REAL", "STRING" };

	if (needrev)
		reverse(&kw->names);
	S = kw->names;
	printf("%s", S->val);
	while((S = S->next))
		printf(" ALIAS %s", S->val);
	if ((k = kw->kind & KWKind_Mask)) {
		printf(" %s", kname[k]);
		if (kw->kind & KWKind_List)
			printf("LIST");
		}
	if (kw->lb)
		printf(" %s %s", kw->kind & KWKind_strictLb ? ">" : ">=", kw->lb);
	if (kw->ub)
		printf(" %s %s", kw->kind & KWKind_strictUb ? "<" : "<=", kw->ub);
	if (kw->defname && wantfuncs & 1) {
		if (kw->init)
			printf(" := %s %s", kw->defname, kw->init);
		else if (kw->cinit)
			printf(" := %s \"%s\"", kw->defname, kw->cinit);
		}
	else {
		if (kw->init)
			printf(" := %s", kw->init);
		else if (kw->cinit)
			printf(" := \"%s\"", kw->cinit);
		}
	if (!(wantfuncs&1))
		goto sufchk;
	if (kw->tag)
		printf(ignoretags ? " TAG \"%s\"" : " /*TAG \"%s\"*/", kw->tag);
	if (kw->funcs) {
		fputs(" {", stdout);
		reverse(&kw->funcs);
		for(S = kw->funcs;;) {
			fputs(S->val, stdout);
			if (!(S = S->next))
				break;
			putchar(',');
			}
		putchar('}');
		}
	if (kw->len)
		printf(" LEN %s", kw->len);
	if (kw->group)
		printf(" GROUP \"%s\"", kw->group);
	if (kw->desc)
		printf(" DESC \"%s\"", kw->desc);
 sufchk:
	if (kw->dylib)
		printf(" DYNLIB \"%s\"", kw->dylib);
	if (suf)
		fputs(suf, stdout);
	if (needrev == 2)
		reverse(&kw->names);
	}

 static KWctx *
kwrev(KWctx **pkw)
{
	KWctx *kw, *kw0, *kw1;
	kw1 = *pkw;
	for(kw0 = 0; (kw = kw1); kw0 = kw) {
		kw1 = kw->next;
		kw->next = kw0;
		}
	return *pkw = kw0;
	}

 static int
StrListcmp(StrList *sl1, StrList *sl2)
{
	if (sl1 == sl2)
		return 0;
	while(sl1) {
		if (!sl2 || strcmp(sl1->val, sl2->val))
			return 1;
		sl1 = sl1->next;
		sl2 = sl2->next;
		}
	return sl2 != 0;
	}

#if 0	/* Could hash keywords and use the following samekw to find duplicate  */
	/* contained keywords, noting duplicates when found in a new, normally */
	/* NULL pointer in KWctx, but so far this does not appear worthwhile.  */
#define OffsetOf(T,f) ((char*)&((T*)0)->f - (char*)0)
#define NumberOf(x) (sizeof(x)/sizeof(x[0]))
#define ValueOf(T,x,o) (*(T*)((char*)x + o))

 static int
samekw(KWctx *kw1, KWctx *kw2)
{
	char *s1, *s2;
	size_t *sp, *spe;

	static size_t
		Intoff[] = {
			OffsetOf(KWctx,alt),
			OffsetOf(KWctx,kind),
			OffsetOf(KWctx,nalt),
			OffsetOf(KWctx,nkw),
			OffsetOf(KWctx,nreq),
			OffsetOf(KWctx,req)
			},
		Stroff[] = {
			OffsetOf(KWctx,lb),
			OffsetOf(KWctx,ub),
			OffsetOf(KWctx,init),
			OffsetOf(KWctx,cinit),
			OffsetOf(KWctx,defname),
			OffsetOf(KWctx,desc),
			OffsetOf(KWctx,dylib),
			OffsetOf(KWctx,group),
			OffsetOf(KWctx,len)
			};
	while(kw1) {
		if (!kw2)
			return 0;
		if (kw1->pparent != kw2->pparent
		 || kw1->rparent != kw2->rparent)
			return 0;
		if (StrListcmp(kw1->names, kw2->names)
		 || StrListcmp(kw1->funcs, kw2->funcs))
			return 0;
		for(sp = Intoff, spe = sp + NumberOf(Intoff); sp < spe; ++sp)
			if (ValueOf(int,kw1,*sp) != ValueOf(int,kw2,*sp))
				return 0;
		for(sp = Stroff, spe = sp + NumberOf(Stroff); sp < spe; ++sp) {
			s1 = ValueOf(char*,kw1,*sp);
			s2 = ValueOf(char*,kw2,*sp);
			if (s1 != s2 && (!s1 || !s2 || strcmp(s1,s2)))
				return 0;
			}
		kw1 = kw1->next;
		kw2 = kw2->next;
		}
	return kw2 == 0;
	}
#undef ValueOf
#undef NumberOf
#undef OffsetOf
#endif

 static void
paradj(int needrev)
{
	KWctx *kw, *kw1, *kw2, **kwp, **kwtop, **kwbot, *plist;
	int lev;

	kwtop = kwbot = KWStack;
	*kwtop = &KWbase;
	if (needrev)
		kwrev(kwtop);

	/* remove extra parens */

	for(;; --kwtop) {
		for(kw = *kwtop; kw; kw = kw->next) {
 top:
			if (needrev && (kw1 = kw->kw) && kw1->next)
				kwrev(&kw->kw);
			if ((kw1 = kw->pparent) && kw1->parstate != Parchecked) {
				kwp = &kw->pparent;
				plist = kw;
				for(;;) {
					if (kw1->kw) {
						*kwtop++ = kw1;
						kw1->parstate = Parchecked;
						*kwp = kw1;
						kwp = &kw1->pparent;
						kw1->next = plist;
						plist = kw1;
						}
					if (!(*kwp = kw1->pparent)) {
						lev = kw1->level;
						break;
						}
					kw1 = *kwp;
					if (kw1->parstate == Parchecked) {
						lev = kw1->level + 1;
						break;
						}
					}
				for(;; plist = kw1) {
					plist->level = lev++;
					if (plist == kw)
						break;
					kw1 = plist->next;
					plist->next = 0;
					}
				}
			if (kw->kw) {
				if ((*kwtop = kw->next))
					++kwtop;
				kw = kw->kw;
				goto top;
				}
			}
		if (kwtop <= kwbot)
			break;
		}
	}

 static KWctx *
masteradj(KWctx *kw)
{
	/* to help print (nm1 | nm2 foo goo), make nm2 "master" */
	KWctx *kw1, *kw2, *kw3, *kwm, *rv;

	for(kw1 = kw, kw3 = 0; kw1; kw1 = kw1->next) {
		if ((kw2 = kw1->master))
			kw3 = kw2->master = kw2;
		}
	if (!kw3)
		return kw;
	for(rv = kw; kw; kw = kw->next) {
		kwm = kw->master;
		if ((kw1 = kw->next) && kw1->master == kwm && kwm) {
			kw3 = kw->kw ? kw : 0;
			do {
				kw2 = kw1;
				if (kw1->kw) {
					if (kw3)
						botch("Two masters in masteradj!");
					kw3 = kw1;
					}
				}
				while((kw1 = kw1->next) && kw1->master == kwm);
			if (kw3 || kw2->pparent) {
				if (kw2 != kw3) {
					kw2->pparent = kw3->pparent;
					kw2->kw = kw3->kw;
					kw3->kw = 0;
					}
				for(;;kw = kw->next) {
					kw->master = kw2;
					if (kw == kw2)
						break;
					}
				}
			else
				for(;;kw = kw->next) {
					kw->master = 0;
					if (kw == kw2)
						break;
					}
			kw = kw2;
			}
		}
	return rv;
	}

 static void
xmasteradj(KWctx **x, int n)
{
	/* variant of masteradj for pretty_ex1 */
	KWctx *kw, *kw1, *kw2, *kw3, *kwm;
	int i, j;

	for(i = 0, kw3 = 0; i < n; ++i) {
		kw1 = x[i];
		if ((kw2 = kw1->master))
			kw3 = kw2->master = kw2;
		}
	if (!kw3)
		return;
	for(i = 0; i < n; i = j) {
		j = i + 1;
		kw = x[i];
		kwm = kw->master;
		if (j < n && kwm && (kw1 = x[j])->master == kwm) {
			kw3 = kw->kw ? kw : 0;
			do {
				kw2 = kw1;
				if (kw1->kw) {
					if (kw3)
						botch("Two masters in xmasteradj!");
					kw3 = kw1;
					}
				}
				while(++j < n && (kw1 = x[j])->master == kwm);
			if (kw3 || kw2->pparent) {
				if (kw2 != kw3) {
					kw2->pparent = kw3->pparent;
					kw2->kw = kw3->kw;
					kw3->kw = 0;
					}
				for(;;kw = x[++i]) {
					kw->master = kw2;
					if (kw == kw2)
						break;
					}
				}
			else
				for(;; kw = x[++i]) {
					kw->master = 0;
					if (kw == kw2)
						break;
					}
			}
		}
	}

 static void
do_indent(int indent)
{
	for(indent += 8; indent >= 8; indent -= 8)
		putchar('\t');
	for(; indent > 0; --indent)
		putchar(' ');
	}

 static void
do_indent1(int indent)	/* write indent-1 spaces */
{
	while(--indent > 0)
		putchar(' ');
	}

 static const char *KEYKIND[4] = { "KEYWORD ", "KEYWORD01 ", "KEYWORD1 ", "KEYWORD12 " };

 static void
pretty_print(int iinc, int needrev)
{
	KWctx *kw, *kw1, *kw2, *kw3, **kwtop, **kwbot, *plist;
	const char *s;
	int alt, indent, lev, needindent, needspace;

	kwtop = kwbot = KWStack;
	indent = 0;
	plist = 0;
	needindent = 1;
	for(kw1 = KWbase.kw; kw1;) {
		fputs(KEYKIND[(kw1->kind >> KWtopshift) & 3], stdout);
		prkwname(kw1, "\n", needrev);
		*kwtop = masteradj(kw1->kw);
		alt = 0;
		for(;;--kwtop) {
	 top:
			for(kw = *kwtop; kw; kw = kw->next) {
				lev = kw->level;
				if (needindent) {
					do_indent(indent);
					needindent = 0;
					}
				if (plist
				 && (plist->level > lev
				  || (plist->level == lev && plist != kw->pparent))) {
 finishpar:
					*kwtop = kw;
					*++kwtop = plist;
					*++kwtop = masteradj(plist->kw);
					plist->alt = alt;
					alt = 0;
					plist = plist->next;
					goto top;
					}
				if (kw->parstate == Printing) {
					s = !kw->names || /*[(*/kw->req ? ")\n" : "]\n";
					if (kw->kw && !kw->req
					 && kw->next && kw->next->alt == kw->alt) {
						s = /*(*/ ")\n";
						indent -= iinc;
						}
					fputs(s, stdout);
					needindent = 1;
					alt = kw->alt;
					indent -= iinc;
					if (!kw->names) {
						alt = kw->alt;
						goto bot;
						}
					continue;
					}
				if (plist)
					for (kw2 = kw->pparent;
							kw2 && kw2->parstate != Printing;
							kw2 = kw2->pparent) {
						if (plist->level >= kw2->level)
							goto finishpar;
						}
				needindent = 1;
				needspace = 0;
				if (alt == kw->alt) {
					putchar('|');
					needspace = 1;
					}
				if ((kw2 = kw->pparent) && kw2->parstate != Printing) {
					kw3 = 0;
					do {
						if (kw2->parstate == Printing)
							break;
						if (kw2->kw) {
							kw2->next = kw3;
							kw3 = kw2;
							kw2->parstate = Printing;
							if (needspace) {
								putchar('\n');
								needspace = 0;
								do_indent(indent);
								}
							indent += iinc;
							putchar('('); /*)*/
							do_indent1(iinc);
							}
						}
						while((kw2 = kw2->pparent));
					while(kw3) {
						kw2 = kw3->next;
						kw3->next = plist;
						plist = kw3;
						kw3 = kw2;
						}
					}
				if ((kw2 = kw->master) && kw == kw2)
					kw2 = 0;
				if (!kw->kw && (!kw2 || kw2->parstate == Leftput)) {
					if (needspace)
						putchar(' ');
					if (!kw->req) {
						if (alt != kw->alt) {
							fputs("[ ", stdout);
							alt = kw->alt;
							}
						s = /*[*/ " ]\n";
						if ((kw2 = kw->next) && kw2->alt == alt)
							s = "\n";
						prkwname(kw, s, needrev);
						continue;
						}
					alt = kw->alt;
					prkwname(kw, "\n", needrev);
					continue;
					}
				if (kw->parstate == Leftput) {
					if (needspace)
						goto putspace;
					}
				else {
					if (needspace) {
						putchar('\n');
						do_indent(indent);
						}
					putchar(kw->req ? '(' : '['); /*])*/
					if (kw->kw && !kw->req
					 && kw->next && kw->next->alt == kw->alt) {
						do_indent1(iinc);
						putchar('('/*)*/);
						indent += iinc;
						}
					indent += iinc;
 putspace:
					do_indent1(iinc);
					}
				kw->parstate = Printing;
				prkwname(kw, "\n", needrev);
				if (kw2) {
					alt = kw->alt;
					kw2->parstate = Leftput;
					continue;
					}
				alt = 0;
				*kwtop = kw;
				*++kwtop = masteradj(kw->kw);
				goto top;
				}
 bot:
			if (kwtop == kwbot) {
				if (plist) {
					*kwtop++ = plist;
					*kwtop++ = masteradj(plist->kw);
					plist->alt = alt;
					alt = 0;
					plist = plist->next;
					}
				else {
					if (indent) {
						if (needindent)
							do_indent(indent);
						fputs(/*(*/")\n", stdout);
						needindent = 1;
						indent -= iinc;
						}
					break;
					}
				}
			}
		if (!(kw1 = kw1->next))
			break;
		putchar('\n');
		}
	}

 static void
nslist(StrList *h, const char *fmt)
{
	for(; h; h = h->next)
		printf(fmt, h->val);
	}

 static void
Output(StrList *headers, StrList *nsname, int hwant, const char *kname, int libmode, int nwant)
{
	KWctx *kw, *kw1, *kw2, **kwtop, **kwbot;
	StrList *h;
	Uint xflag;
	const char *fmt;
	int kwkno;
	static char
		libfmt1[] = "#ifdef __cplusplus\nextern \"C\" {\n"
			"KeyWord *keyword_add(void);\n}\n#endif\n\n"
			" KeyWord*\nkeyword_add(void) {\n\t",
		libfmt2[] = "kw_%d[0].kw;\n\t}\n",
		libfmt3[] = /*{*/ "kw_%d.kw;\n\t}\n";

	if (hwant && headers) {
		for(h = reverse(&headers); h; h = h->next)
			printf(*h->val == '<' /*>*/
				? "#include %s\n"
				: "#include \"%s\"\n", h->val);
		}
	if (nwant) {
		if (nsname)
			reverse(&nsname);
		else
			nsname = new_StrList("Dakota", nsname);
		putchar('\n');
		nslist(nsname, "namespace %s {\n" /*}*/);
		}
	/* BMA (20160923): Suppress keyword count summary to avoid merge conflicts */
	//printf("\n/** %d distinct keywords (plus %d aliases) **/\n", lastagroup, nalias);
	printf("\nstatic %sKeyWord%s\n", guikeywds ? "Gui" : "", libmode ? "x" : "");
	kwtop = kwbot = KWStack;
	*kwtop = &KWbase;
	if (libmode) {
		xflag = KWKind_Extended;
		Comma = Commasx;
		}
	else {
		xflag = 0;
		Comma = Commas0;
		}
	for(;;--kwtop) {
 top:
		for(kw = *kwtop; kw; kw = kw->next) {
			if ((kw->master && kw->master != kw)
			|| ((!(kw1 = kw->pparent) || kw->parstate == KwStacked) && !kw->kw)
			|| kw->kwkno)
				continue;
			if (kw1 && kw1->parstate != KwStacked) {
				if ((kw1 = kw->pparent) && kw1->parstate != KwStacked) {
					kw1->parstate = KwStacked;
					*kwtop = kw;
					*++kwtop = kw1;
					kw2 = *++kwtop = kwsort(kw1);
					do kw2->parstate = KwStacked;
						while((kw2 = kw2->next));
					while((kw1 = kw1->pparent) && kw1->parstate != KwStacked) {
						kw1->parstate = KwStacked;
						*++kwtop = kw1;
						kw2 = *++kwtop = kwsort(kw1);
						do kw2->parstate = KwStacked;
							while((kw2 = kw2->next));
						}
					goto top;
					}
				if (!kw1 && !kw->kw)
					continue;
				}
			if (!kw->nkw && kw->kw) {
				*kwtop = kw;
				*++kwtop = kwsort(kw);
				goto top;
				}
			kwout(kw, kwtop, xflag);
			}
		if (kwtop == kwbot)
			break;
		}
	printf(";\n\n");
	if (nwant) {
		nslist(reverse(&nsname), /*{*/"} // namespace %s\n");
		putchar('\n');
		reverse(&nsname);
		}

	printf("#ifdef __cplusplus\n");
	printf("extern \"C\" {\n");
	printf("#endif\n");

	kwkno = KWbase.kwkno;
	if (libmode) {
		if (KWbase.nkw > 1 && !kname)
			kname = "KeywordTop";
		printf(libfmt1);
		fmt = libfmt2;
		if (kname) {
			printf("static KeyWordx kw_%d = { \"%s\", 0x%x,%d,%d,%d,(KeyWord*)",
				kwkno + 1, kname, KWKind_Dynmult, KWbase.nkw,
				KWbase.alt, KWbase.req);
			if (nwant) {
				nslist(nsname, "%s::");
				nwant = 0;
				}
			printf("kw_%d };\n\t", kwkno++);
			fmt = libfmt3;
			}
		printf("return &");
		}
	else {
		printf("KeyWord Dakota_Keyword_Top = {\"KeywordTop\",0,%d,%d,%d,",
			KWbase.nkw, KWbase.alt, KWbase.req);
		fmt = "kw_%d};\n";
		}
	if (nwant)
		nslist(nsname, "%s::");
	printf(fmt, kwkno);

	printf("#ifdef __cplusplus\n");
	printf("}\n");
	printf("#endif\n");

	if (startwarn && missing && !guikeywds)
		fprintf(stderr, "%s: %d start routines implicitly named.\n",
			progname, missing);
	}

 typedef struct Valhash Valhash;

 struct
Valhash {
	DefhashHead hh;
	Uint n;
	Valhash *thread;
	};

 static Kwpack *kwpack, *kwpack0, *kwpacke;
 static size_t kwpack_len, strtab_len = 1;
 static DefhashTab FpTab, StrTab;
 static Valhash *FpTabLast, *StrTabLast;

 static Uint
NumOff(char *name)
{
	Valhash *v;
	if (!name || !*name || (*name == '0' && !name[1]))
		return 0;
	v = (Valhash*)DHHlookup(name, &FpTab, 1);
	if (!v->n) {
		v->n = FpTab.nent;
		v->thread = FpTabLast;
		FpTabLast = v;
		}
	return v->n;
	}

 static Uint
StrOff(char *name)
{
	Valhash *v;
	if (!name)
		return 0;
	v = (Valhash*)DHHlookup(name, &StrTab, 1);
	if (!v->n) {
		v->n = strtab_len;
		strtab_len += strlen(name) + 1;
		v->thread = StrTabLast;
		StrTabLast = v;
		}
	return v->n;
	}

 static void
kwjout(KWctx *kw, KWctx **kwtop)
{
	KWctx *kw0, *kw1, *kw2, *kwm, **x, *x0[KWOUTlen0];
	Kwpack *p;
	Uint lastknoj;
	int i, m, n;
	size_t L;

	kw1 = kw->kw;
	n = 0;
	lastknoj = kwpack - kwpack0 + 1;
	if (kw->names) {
		kw0 = kw1;
		if (kw->parstate != KwStacked)
			for(kw2 = kw->pparent; kw2; kw2 = kw2->pparent) {
				kw2->next = kw1;
				kw1 = kw2;
				++n;
				}
		if (!kw0) {
			if (!--n) {
				kw->nkw = kw1->nkw;
				if (kw1->kwkno) {
					kw->kwkno = kw1->kwkno;
					return;
					}
				kw->kwkno = lastknoj;
				kw = kw1;
				kw1 = kw->kw;
				}
			else {
				kw2 = kw->pparent;
				if (kw2->kwknop) {
					kw->kwkno = kw2->kwknop;
					return;
					}
				kw2->kwknop = lastknoj;
				kw1 = kw2;
				while((kw2 = kw2->pparent)) {
					kw2->next = kw1;
					kw1 = kw2;
					}
				}
			}
		}
	for(m = 0, kw0 = kw1; kw0; kw0 = kw0->next)
		++m;
	x = x0;
	if (m > KWOUTlen0)
		x = (KWctx**)Malloc(m*sizeof(KWctx*), "kwjout");
	m = 0;
	for(kw0 = kw1; kw0; kw0 = kw0->next)
		x[m++] = kw0;
	kw->kwkno = lastknoj;
	++lastkno;
	for(i = 0; i < m;) {
		if (++kwpack >= kwpacke) {
			L = 2*kwpack_len;
			kwpack0 = (Kwpack*)Realloc(kwpack0, L*sizeof(Kwpack), "kwjout");
			kwpacke = kwpack0 + L;
			kwpack = kwpack0 + kwpack_len;
			memset(kwpack, 0, kwpack_len*sizeof(Kwpack));
			kwpack_len = L;
			}
		p = kwpack;
		kw1 = x[i++];
		if (!(kwm = kw1->master))
			kwm = kw1;
		p->nkw = kwm->nkw;
		if (kw1->names) {
			p->name = StrOff(kw1->names->val);
			p->kind = kw1->kind | kw1->names->primary;
			p->alt = kwm->alt + kwm->altoffset;
			if (kwm->req)
				p->req = kwm->req + kwm->reqoffset;
			if (kwm->dylib) {
				p->dylib = StrOff(kwm->dylib);
				p->kind |= KWKind_Dynlib;
				}
			}
		p->kw = kwm->kwkno;
		p->Lb = NumOff(kw1->lb);
		p->Ub = NumOff(kw1->ub);
		p->poff = kw1->paoff;
		}
	if (x != x0)
		free(x);
	}

 static Valhash *
Vreverse(Valhash *sh)
{
	Valhash *sh0, *sh1;
	for(sh0 = 0; sh; sh = sh1) {
		sh1 = sh->thread;
		sh->thread = sh0;
		sh0 = sh;
		}
	return sh0;
	}

 static void
kwpack_adjust(Kwpack *kp, Uint nkw)
{
	Kwpack0 *qp, *qpe;

	qp = (Kwpack0*)kp;
	for(qpe = qp + nkw; qp < qpe; ++kp, ++qp) {
		memcpy(qp, kp, sizeof(Kwpack0) - sizeof(int));
		qp->poff = kp->poff;
		}
	}

 static void
jOutput(char *jfname)
{
	FILE *f;
	KWctx *kw, *kw1, *kw2, **kwtop, **kwbot;
	Kwpack *kp;
	Uint nkw0;
	Valhash *vh;
	double d;
	size_t L, Lkw, Ls;
	static KwpHead head = { "DakreorderKwds\n", 1e13, 0x1234abcd };
	static double zbuf[2];

	FpTab.nalloc  = (sizeof(Valhash) + sizeof(void*) - 1)/sizeof(void*);
	StrTab.nalloc = (sizeof(Valhash) + sizeof(void*) - 1)/sizeof(void*);
	kwpack_len = 128;
	kwpack = kwpack0 = (Kwpack*)Malloc(kwpack_len * sizeof(Kwpack), "jOutput");
	memset(kwpack, 0, kwpack_len * sizeof(Kwpack));
	kwpacke = kwpack0 + kwpack_len;
	kwtop = kwbot = KWStack;
	*kwtop = &KWbase;
	for(;;--kwtop) {
 top:
		for(kw = *kwtop; kw; kw = kw->next) {
			if ((kw->master && kw->master != kw)
			|| ((!(kw1 = kw->pparent) || kw->parstate == KwStacked) && !kw->kw)
			|| kw->kwkno)
				continue;
			if (kw1 && kw1->parstate != KwStacked) {
				if ((kw1 = kw->pparent) && kw1->parstate != KwStacked) {
					kw1->parstate = KwStacked;
					*kwtop = kw;
					*++kwtop = kw1;
					kw2 = *++kwtop = kwsort(kw1);
					do kw2->parstate = KwStacked;
						while((kw2 = kw2->next));
					while((kw1 = kw1->pparent) && kw1->parstate != KwStacked) {
						kw1->parstate = KwStacked;
						*++kwtop = kw1;
						kw2 = *++kwtop = kwsort(kw1);
						do kw2->parstate = KwStacked;
							while((kw2 = kw2->next));
						}
					goto top;
					}
				if (!kw1 && !kw->kw)
					continue;
				}
			if (!kw->nkw && kw->kw) {
				*kwtop = kw;
				*++kwtop = kwsort(kw);
				goto top;
				}
			kwjout(kw, kwtop);
			}
		if (kwtop == kwbot)
			break;
		}
	StrTabLast = Vreverse(StrTabLast);
	FpTabLast = Vreverse(FpTabLast);
	nkw0 = kwpack - kwpack0;
	head.nkw = nkw0 + 1;
	kp = kwpack0;
	kp->kw = nkw0 - KWbase.nkw + 1;
	kp->nkw = KWbase.nkw;
	kp->alt = KWbase.alt;
	kp->req = KWbase.req;
	head.strtab_offset = (FpTab.nent+1)*sizeof(double);
	Ls = (strtab_len + 7) & ~7;
	head.kw_offset = head.strtab_offset + Ls;
	if ((head.pad = dynlibseen))
		Lkw = sizeof(Kwpack);
	else {
		Lkw = sizeof(Kwpack0);
		kwpack_adjust(kp, head.nkw);
		}
	head.end_offset = head.kw_offset + head.nkw * Lkw;
	f = fopen(jfname, "wb");
	if (!f) {
		fprintf(stderr, "Could not open \"%s\"\n", jfname);
		exit(1);
		}
	fwrite(&head, sizeof(head), 1, f);
	for(vh = FpTabLast; vh; vh = vh->thread) {
		d = strtod(vh->hh.name,0);
		fwrite(&d, sizeof(double), 1, f);
		}
	fwrite(&zbuf,9,1,f);
	for(vh = StrTabLast; vh; vh = vh->thread)
		fwrite(vh->hh.name, strlen(vh->hh.name) + 1, 1, f);
	if ((L = Ls - strtab_len))
		fwrite(&zbuf, L, 1, f);	/* align */
	fwrite(kp, Lkw, head.nkw, f);
	fclose(f);
	}

 typedef int (*KWcomp)(const void*, const void*);

 static const char **zalt;

 static int
kwcomp1(const void *a, const void *b)
{
	const KWctx *ka, *kb;

	ka = *(const KWctx**)a;
	kb = *(const KWctx**)b;
	if (ka->req) {
		if (!kb->req)
			return -1;
		}
	else if (kb->req)
		return 1;
	if (ka->alt == kb->alt) {
		if (ka->master != kb->master) {
			if (!ka->master)
				return -1;
			if (!kb->master)
				return 1;
			ka = ka->master;
			kb = kb->master;
			}
		return strcmp(ka->names->val, kb->names->val);
		}
	return strcmp(zalt[ka->alt], zalt[kb->alt]);
	}

 static int
kwcomp2(const void *a, const void *b)
{
	const KWctx *ka, *kb, *pa, *pb;

	ka = *(const KWctx**)a;
	kb = *(const KWctx**)b;
	if (ka->req) {
		if (!kb->req)
			return -1;
		}
	else if (kb->req)
		return 1;
	if (ka->alt == kb->alt) {
		pa = ka->master ? ka->master->pparent : ka->pparent;
		pb = kb->master ? kb->master->pparent : kb->pparent;
		if (pa != pb)
			return ka->objno - kb->objno;
		if (ka->master != kb->master) {
			if (!ka->master)
				return -1;
			if (!kb->master)
				return 1;
			ka = ka->master;
			kb = kb->master;
			}
		return strcmp(ka->names->val, kb->names->val);
		}
	return strcmp(zalt[ka->alt], zalt[kb->alt]);
	}

 static void
do_sort(KWctx *kw, KWcomp KWcmp)
{
	KWctx *kw1, *kwm, **x, *x0[KWPlen];
	const char *z0[KWPlen];
	int i, n, na;

	for(n = 0, kw1 = kw->kw; kw1; kw1 = kw1->next) {
		++n;
		reverse(&kw1->names);
		}
	if (n <= 1)
		return;
	x = n > KWPlen ? (KWctx**)Malloc(n*sizeof(KWctx*), "do_sort") : x0;
	na = kw->nalt + 1;
	zalt = na <= KWPlen ? z0 : (const char**)Malloc(na*sizeof(const char*), "do_sort");
	memset(zalt, 0, na*sizeof(const char*));
	for(n = 0, kw1 = kw->kw; kw1; kw1 = kw1->next) {
		x[n++] = kw1;
		if ((kwm = kw1->master))
			kwm->master = kwm;
		if (!zalt[i = kw1->alt]
		 || strcmp(kw1->names->val, zalt[i]) < 0)
			zalt[i] = kw1->names->val;
		}
	qsort(x, n--, sizeof(KWctx*), KWcmp);
	for(i = 0; i < n; ++i)
		x[i]->next = x[i+1];
	x[n]->next = 0;
	kw->kw = x[0];
	if (zalt != z0)
		free(zalt);
	if (x != x0)
		free(x);
	}

 static void
ksort(KWcomp KWcmp)
{
	KWctx *kw, *kw1, **kwtop, **kwtop0, **kwbot;

	kwtop = kwbot = KWStack;

	for(*kwtop = &KWbase;; --kwtop) {
 top:
		kw = *kwtop;
 top1:
		while(kw) {
			if ((kw1 = kw->pparent) && kw1->parstate != Sorted) {
				*kwtop = kw;
				kwtop0 = kwtop;
				do {
					kw1->parstate = Sorted;
					*++kwtop = kw1;
					}
					while((kw1 = kw1->pparent) && kw1->parstate != Sorted);
				if (kwtop > kwtop0)
					goto top;
				}
			if (kw->kw) {
				do_sort(kw, KWcmp);
				if ((*kwtop = kw->next))
					++kwtop;
				kw = *kwtop = kw->kw;
				goto top1;
				}
			kw = kw->next;
			}
		if (kwtop <= kwbot)
			break;
		}
	}

 static void
pretty_ex1(KWctx *kw, KWctx *prevkw, KWctx *nextkw, int indent, int iinc, int needrev)
{
	KWctx *kw0, *kw1, *kw2, *kwm, *kwp, *kwp0, **x, *y, *y1;
	const char *suf, *z0[KWPlen];
	int ep, i, i0, ia, j, n, na, np;

	n = 0;
	x = 0;
	y = 0;
	kw1 = kw->kw;
	if (needrev == 1)
		needrev = 2;
	if ((kwp0 = kwp = kw->pparent) && (!kw->master || kw->master == kw)) {
		for(;;) {
			if (kwp->parstate != Printing) {
				kwp->parstate = Printing;
				for(kw1 = kwp->kw; kw1; kw1 = kw1->next)
					kw1->pparent = 0;	/* avoid infinite loops */
				}
			for(kw1 = kwp->kw; kw1; kw1 = kw1->next)
				++n;
			if (!(kw1 = kwp->pparent))
				break;
			kwp = kwp->pparent = kw1;
			}
		np = n;
		for(kw1 = kw->kw; kw1; kw1 = kw1->next)
			++n;
		y = y1 = (KWctx*)Malloc(np*sizeof(KWctx) + n*sizeof(KWctx*), "pretty_ex1");
		x = (KWctx**)(y + np);
		for(i = 0, kw1 = kw->kw; kw1; kw1 = kw1->next)
			x[i++] = kw1;
		i0 = i;
		for(kwp = kwp0; kwp; kwp = kwp->pparent) {
			j = i;
			na = i + kwp->nalt;
			for(kw1 = kwp->kw; kw1; kw1 = kw1->next, ++y1) {
				memcpy(y1, kw1, sizeof(KWctx));
				x[i++] = y1;
				y1->alt += j;
				}
			}
		if (!needrev) {
			for(j = i0; j < n; ++j) {
				kw1 = x[j];
				kw1->nreq = -(kw1->nreq + 1);
				}
			zalt = ++na <= KWPlen ? z0
				: (const char**)Malloc(na*sizeof(const char*), "pretty_ex1");
			memset(zalt, 0, na*sizeof(const char*));
			for(j = 0; j < n; ++j) {
				kw1 = x[j];
				if ((kwm = kw1->master))
					kwm->master = kwm;
				if (!zalt[i = kw1->alt]
				 || strcmp(kw1->names->val, zalt[i]) < 0)
					zalt[i] = kw1->names->val;
				}
			qsort(x, n, sizeof(KWctx*), kwcomp1);
			if (zalt != z0)
				free(zalt);
			for(j = 0; j < n; ++j) {
				kw1 = x[j];
				if (kw1->nreq < 0) {
					kw1->nreq = -(kw1->nreq + 1);
					kw1->pparent = 0;
					}
				}
			}
		kw1 = x[0];
		xmasteradj(x, n);
		}
	if (kw1) {
		if (!x)
			kw1 = masteradj(kw->kw);
		ep = kw->req;
		if (prevkw && kw->alt == prevkw->alt) {
			ia = kw->master && kw->master == prevkw->master ? iinc : 0;
			do_indent(indent + ia);
			if (!prevkw->master || prevkw->master != kw->master) {
				fputs("|\n", stdout);
				do_indent(indent + ia);
				putchar('(');
				}
			else
				putchar('|');
			}
		else {
			do_indent(indent);
			if (kw->req)
				putchar('('/*)*/);
			else {
				putchar('[');
				if (nextkw && !nextkw->req && kw->alt == nextkw->alt) {
					ep = 1;
					do_indent1(iinc);
					putchar('(');
					indent += iinc;
					}
				}
			}
		do_indent1(iinc);
		indent += iinc;
		prkwname(kw, "\n", needrev);
		i = 0;
		kw0 = 0;
		while(kw1) {
			if (x)
				kw2 = ++i >= n ? 0 : x[i];
			else
				kw2 = kw1->next;
			pretty_ex1(kw1, kw0, kw2, indent, iinc, needrev);
			kw0 = kw1;
			kw1 = kw2;
			}
		if (x)
			free(y);
		do_indent(indent);
		fputs(ep ? ")\n" : "]\n", stdout);
		}
	else {
		if (prevkw && kw->alt == prevkw->alt
		 && kw->master && kw->master == prevkw->master)
			indent += iinc;
		do_indent(indent);
		suf = "\n";
		if (prevkw && kw->alt == prevkw->alt) {
			if (kw->master && kw->master != prevkw->master) {
				fputs("|\n", stdout);
				do_indent(indent);
				fputs("( ", stdout);
				}
			else {
				fputs("| ", stdout);
				if (!kw->req)
					goto rbcheck;
				}
			}
		else if (kw->req) {
			if (prevkw) {
				if (prevkw->master != kw->master)
					goto putlpar;
				}
			else if (kw != kw->master && kw->master) {
 putlpar:
				fputs("( ", stdout);
				}
			}
		else {
			fputs("[ ", stdout);
 rbcheck:
			if (!nextkw || nextkw->alt != kw->alt)
				suf = " ]\n";
			}
		prkwname(kw, suf, needrev);
		}
	}

 static void
pretty_expand(int iinc, int needrev)
{
	KWctx *kw, *kw0, *kw1, *kw2;

	if ((kw = KWbase.kw)) for(;;) {
		masteradj(kw->kw);
		fputs(KEYKIND[(kw->kind >> KWtopshift) & 3], stdout);
		prkwname(kw, "\n", needrev);
		kw0 = 0;
		for(kw1 = kw->kw; kw1; kw1 = kw2) {
			kw2 = kw1->next;
			pretty_ex1(kw1, kw0, kw2, 0, iinc, needrev);
			kw0 = kw1;
			}
		if (!(kw = kw->next))
			break;
		putchar('\n');
		}
	}

 static int
Str_option(char *s, char ***pav, int *want, StrList **ps)
{
	if (*++s || (s = *++*pav)) {
		if (*s == '-' && ! s[1])
			*want = 0;
		else
			*ps = new_StrList(s, *ps);
		return 0;
		}
	return 1;
	}

 static int
version(void)
{
	printf("nidrgen version 20100528\n");
	return 0;
	}

 static int
tagcomp(const void *a, const void *b)
{
	return strcmp(*(const char**)a, *(const char**)b);
	}

 static int
unused_tags()
{
	Taghash *th;
	char **x;
	size_t i, nbad;

	nbad = 0;
	for(th = (Taghash*)TagTab.lastdh; th; th = (Taghash*)th->hh.thread)
		if (!th->uses)
			++nbad;
	if (!nbad)
		return 0;
	x = (char**)Malloc(nbad*sizeof(char*), "unsed_tags");
	for(i = 0, th = (Taghash*)TagTab.lastdh; th; th = (Taghash*)th->hh.thread)
		if (!th->uses)
			x[i++] = th->hh.name;
	fprintf(stderr, "%lu unused tag definition%s in %s:\n",
		(unsigned long)nbad, "s" + (nbad == 1), tagfname);
	if (nbad > 1)
		qsort(x, nbad, sizeof(char*), tagcomp);
	for(i = 0; i < nbad; ++i)
		fprintf(stderr, "\t\"%s\"\n", x[i]);
	free(x);
	return 1;
	}

 static int
defcomp(const void *a, const void *b)
{
	return strcmp((*(const Defhash**)a)->hh.name, (*(const Defhash**)b)->hh.name);
	}

 static void
write_defs(FILE *f)
{
	Defhash *d, **x, **x1;
	DefhashHead *dh;
	size_t i, n;

	n = DHTab.nent;
	x = x1 = (Defhash**)Malloc(n*sizeof(Defhash*), "write_defs");
	for(dh = DHTab.lastdh, i = 0; dh; dh = dh->thread) {
		if (++i > n)
			botch("lastdef count botch");
		*x1++ = (Defhash*)dh;
		}
	qsort(x, n, sizeof(Defhash*), defcomp);
	for(i = 0; i < n; ++i) {
		d = x[i];
		fprintf(f, d->qkeep ? "#define %s \"%s\"\n"
			: "#define %s %s\n",
			d->hh.name, d->value);
		}
	}

 int
main(int argc, char **argv)
{
	FILE *f;
	StrList *headers, *mtagname, *nsname, *outtagname;
	char *defsout, *jfname, *ofname, *s, *se;
	const char *kname, *libex, *libname;
	int defwrite, dupcheck, hwant, libmode, needrev, nwant, pprint, sort, x;

	KWStack[0] = &KWbase;
	DHTab.nalloc  = (sizeof(Defhash) + sizeof(void*) - 1)/sizeof(void*);
	TagTab.nalloc = (sizeof(Taghash) + sizeof(void*) - 1)/sizeof(void*);
	KWbase.names = new_StrList("$top",0);
	progname = argv[0];
	if (argc <= 1)
 usage1:	return usage(1);
	headers = mtagname = nsname = outtagname = 0;
	kname = libex = libname = defsout = jfname = ofname = 0;
	dupcheck = hwant = needrev = nwant = 1;
	defwrite = libmode = pprint = sort = 0;
 nextopt:
	while((s = *++argv) && *s == '-' && s[1]) {
		while(*++s)
		switch(*s) {
		 case 'D':
			defwrite = libexpand = 1;
			libex = "-D";
			break;
		 case 'T':
			if (Str_option(s, &argv, &hwant, &outtagname))
				goto usage1;
			goto nextopt;
		 case 'd':
			dupcheck = 0;
			break;
		 case 'e':
			expand = 1;
			break;
		 case 'f':
			wantfuncs = 0;
			break;
		 case 'g':
			guikeywds = libexpand = 1;
			libex = "-g";
			break;
		 case 'G':
			dfltgroup = 0;
			break;
		 case 'h':
			if (Str_option(s, &argv, &hwant, &headers))
				goto usage1;
			goto nextopt;
		 case 'j':
			if (!*++s && !(s = *++argv))
				goto usage1;
			jfname = s;
			libexpand = 1;
			libex = "-j";
			goto nextopt;
		 case 'k':
			if (!*++s && !(s = *++argv))
				goto usage1;
			kname = s;
			libmode = 1;
			goto nextopt;
		 case 'l':
			libmode = 1;
			break;
		 case 'L':
 have_L:
			if (*++s)
				libname = s;
			else if (!(libname = *++argv))
				return usage(1);
			goto nextopt;
		 case 'm':
			if (Str_option(s, &argv, &hwant, &mtagname))
				goto usage1;
			goto nextopt;
		 case 'n':
			if (Str_option(s, &argv, &nwant, &nsname))
				goto usage1;
			goto nextopt;
		 case 't':
			ignoretags = 1;
			break;
		 case 's':
			sort = 1;
			needrev = 0;
			libex = "-s";
			for(;;++s) {
				switch(s[1]) {
				  case 'L': ++s; goto have_L;
				  case 'e': expand = 1;
					    continue;
				  case 'f': wantfuncs = 0;
				  case 'p': continue;
				  }
				break;
				}
			/* no break */
		 case 'p':
			libexpand = 1;
			pprint = 2;
			if (!libex)
				libex = "-p";
			for(;;++s) {
				switch(s[1]) {
				  case 'L':
					++s;
					goto have_L;
				  case 'e':
					expand = 1;
					continue;
				  case 'f':
					wantfuncs = 0;
					continue;
				  case 's':
					sort = 1;
					needrev = 0;
					continue;
				  case 't':
					ignoretags = 1;
					continue;
				  }
				break;
				}
			if (*++s) {
				pprint = (int)strtol(s,&se,10);
				if (pprint < 1 || se)
					goto usage1;
				}
			goto nextopt;
		 case 'v':
			return version();
		 case 'w':
			startwarn = 0;
			break;
		 case '?':
			return usage(s[1] != 0);
		 case '-':
			if (!strcmp(++s,"help"))
				return usage(0);
			if (!strcmp(s,"version"))
				return version();
			if (!*s) {
				s = *++argv;
				goto argsdone;
				}
			/* no break */
		  default:
			goto usage1;
		  }
		}
 argsdone:
	if (libname) {
		if (s && argv[1] && argv[2])
			return usage(1);
		if ((defsout = s))
			ofname = argv[1];
		if (mtagname) {
			fprintf(stderr, "\n%s: ignoreing -m due to -L\n", progname);
			mtagname = 0;
			}
		goto lcheck;
		}
	if (!s || (argv[1] && argv[2] && argv[3] && argv[4]))
		return usage(1);
	if ((tagfname = argv[1])) {
		specfile = s;
		s = tagfname;
		lastseen = Saw_kwname;
		lastkw1 = &KWbase;
		if ((defsout = argv[2]))
			ofname = argv[3];
		}
	else if (mtagout) {
		fprintf(stderr, "%s: -m specified but no tagfile given.\n", progname);
		return 1;
		}
	if (*s == '-' && !s[1])
		yyin = stdin;
	else if (!(yyin = fopen(infname = s, "r"))) {
		cantopen(s);
		return 1;
		}
	if (specfile)
		tagin = yyin;
	if ((mtagname && (!(mtagout = fopen(s = mtagname->val, "w"))))
	 || (outtagname && !(tagout = fopen(s = outtagname->val, "w")))) {
		cantopen(s);
		return 1;
		}
 lcheck:
	if (libmode) {
		wantfuncs |= 2;
		if (libexpand && kname) {
			libmode = 0;
			fprintf(stderr, "\n%s: ignoring -k due to %s\n",
				progname, libex);
			}
		}
	else
		libexpand = 0;
	if (libname)
#ifndef NO_NIDR_DYNLIB
		libread1(0, libname, 0, 2);
#else
		{
		fprintf(stderr, "\ndlopen for \"%s\" is NOT SUPPORTED\n", libname);
		return usage(1);
#endif
		}
	else if ((x = yylex())) {
		fprintf(stderr, "\n%s: Surprise return %d from yylex()\n", progname, x);
		return x;
		}
	if (mtagout)
		fclose(mtagout);
	if (tagout)
		fclose(tagout);
	if (tagout || mtagout)
		goto done;
	if (Br_top > Br_stack) {
		fputs("Missing ", stderr);
		do putc(/*[*/ *Br_top-- == '(' ? ')' : ']', stderr);
			while(Br_top >Br_stack);
		fputs(" at end of file.\n", stderr);
		return 1;
		}
	if (nsquawk > liberror)
		return 1;
	if (tagfname && unused_tags())
		return 1;
	if (KWbase.kw) {
		paradj(needrev);
		if (dupcheck)
			dup_check();
		}
	if (jfname)
		ofname = 0;
	else if (ofname && !freopen(ofname, "w", stdout)) {
		fprintf(stderr, "%s: Cannot open output file \"%s\"\n", progname, ofname);
		return 1;
		}
	if (KWbase.kw) {
		if (pprint) {
			if (expand) {
				if (sort)
					ksort(kwcomp1);
				pretty_expand(pprint, needrev);
				}
			else {
				if (sort)
					ksort(kwcomp2);
				pretty_print(pprint, needrev);
				}
			}
		else if (!defwrite) {
			if (jfname)
				jOutput(jfname);
			else
				Output(headers, nsname, hwant, kname,
					guikeywds ? 0 : libmode, nwant);
			}
		else if (DHTab.lastdh)
			write_defs(stdout);
		}
	if (!defwrite && defsout) {
		f = fopen(defsout,"w");
		if (!f) {
			fprintf(stderr, "%s:  cannot open defs file \"%s\"\n",
				progname, defsout);
			return 1;
			}
		write_defs(f);
		fclose(f);
		}
 done:
	return dupnames > 0 ? 1 : liberror ? 2 : 0;
	}
