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

/* nidr.h */

#ifndef NIDR_H
#define NIDR_H

#include <stddef.h>	/* for size_t;	    for pure C++, change to <cstddef> */
typedef struct KeyWord KeyWord;	/* for C compilation */
typedef struct Values Values;	/* for C compilation */

typedef unsigned int	Uint;

#ifdef NO_DOUBLE
typedef float		Real;
#else
typedef double		Real;
#endif

struct Values {
	size_t	n;	/* for n = 0, pass (Values*)0 to the func. */
	Real	*r;
	int	*i;
	const char **s;
	int rstate;	/* for handling n*value and L:U and L:step:U */
	};

typedef void (*Kwfunc)(const char *keyname, Values *val, void **g, void *v);

enum KeyWordKind {
	KWKind_Void = 0,	/* no associated value */
	KWKind_Int  = 1,	/* INTEGER */
	KWKind_Real = 2,	/* REAL */
	KWKind_Str  = 3,	/* STRING */
	KWKind_Mask = 3,	/* for use with (kind & KWKind_Mask) */
	/* The following may be or-ed into one of the above */
	KWKind_List = 4, /* INTEGERLIST, REALLIST, or STRINGLIST (not for KWKind_Void) */
	KWKind_primary	= 0x08,	 /* primary in its group of aliases */
	KWKind_strictLb	= 0x10,	 /* > value specified */
	KWKind_caneqLb	= 0x20,	 /* >= value specified */
	KWKind_Lb	= 0x30,	 /* Lb mask: > or >= specified */
	KWKind_strictUb	= 0x40, /* < value specified */
	KWKind_caneqUb	= 0x80, /* <= value specified */
	KWKind_Ub	= 0xc0, /* Ub mask: < or <= specified */
	KWKind_01	= 0x100, /* a top-level keyword that can appear at most once */
	KWKind_1	= 0x200, /* a top-level keyword that must appear exactly once */
	KWKind_12	= 0x300, /* a top-level keyword that must appear at least once */
	KWtopshift	= 8,
	KWKind_init	= 0x400, /* has := value */
	KWKind_Len1OK	= 0x800, /* list of values can be one long */
	KWKind_Stacked	= 0x1000, /* internal use (by nidrgen) */
	KWKind_Hashed	= 0x2000, /* internal use */
	KWKind_Dynlib	= 0x4000, /* The final value is NULL and the vf field names a */
				  /* shared library that provides contained keywords */
				  /* and any needed start and final routines. */
	KWKind_Loaded	= 0x8000, /* We've loaded the library; vs points to the */
				  /* substituted KeyWord. */
	KWKind_Dynmult	= 0x10000, /* Library contains several keywords */
	KWKind_Extended = 0x20000, /* Library version of KeyWord:  KeyWord */
				   /* followed by extra GuiKeyWord fields */
	KWKind_Libname	= 0x40000  /* LIBNAME specified: a string value in the */
				   /* input specifies the name of a library to */
				   /* load for contained keywords. */
	};

typedef struct Comment Comment;

typedef struct KWfuncs {
	Kwfunc start;
	void *vs;	/* v arg for start (e.g., pointer to member object) */
	Kwfunc final;
	void *vf;	/* v arg for final */
	} KWfuncs;

struct KeyWord {
	const char *name; /* name of this keyword */
	Uint kind;	/* see enum KeyWordKind */
	Uint nkw;	/* number of keywords nested within this one */
	Uint alt;	/* > 0 ==> this is in alternate group alt */
	Uint req;	/* > 0 ==> this is required item number req */
	KeyWord *kw;	/* the nkw keywords nested within this one */
	Real Lb, Ub;	/* lower and upper bounds, strict or not according to kind */
	int paoff;	/* offset to preferred alias */
	KWfuncs f;
	Comment *comment;
	KeyWord *kwnext;
	};

typedef struct GuiKeyWord GuiKeyWord;
struct GuiKeyWord {
	const char *name; /* name of this keyword */
	Uint kind;	/* see enum KeyWordKind */
	Uint nkw;	/* number of keywords nested within this one */
	Uint alt;	/* > 0 ==> this is in alternate group alt */
	Uint req;	/* > 0 ==> this is required item number req */
	Uint agroup;	/* alias group */
	GuiKeyWord *kw;	/* the nkw keywords nested within this one */
	Real Lb, Ub;	/* lower and upper bounds, strict or not according to kind */
	Real init;
	const char *cinit;
	const char *desc;
	const char *group; /* display group */
	const char *alen;  /* for array types, the keyword specifying the array length */
	const char *dylib; /* "library" containing more of this keyword */
	};

 typedef struct
KeyWordx {
	KeyWord kw;
	Uint seqno, agroup;
	const char *funcs;
	const char *desc;
	const char *group;
	const char *defname;
	const char *alen;
	const char *init;
	const char *cinit;
	} KeyWordx;

 typedef struct KWinfo KWinfo;
 struct
KWinfo {
	KeyWord *kw, *kw1;
	void *g;
	KeyWord **alt, **req;
	int *altct;
	const char *name;
	Uint needstart;
	int nalt, nreq;
	};

/* stuff for dakreorder and nidrgen */

 typedef struct
KwpHead {
	char id[16];
	double fpkind;
	Uint bkind;
	Uint nkw;
	Uint strtab_offset;
	Uint kw_offset;
	Uint end_offset;
	Uint pad; /* so we have an even number of Uint values */
	} KwpHead;

 typedef struct
Kwpack {
	Uint name;
	Uint kind;
	Uint nkw;
	Uint alt;
	Uint req;
	Uint kw;
	Uint Lb;
	Uint Ub;
	Uint dylib;
	int poff;
	} Kwpack;

 typedef struct
Kwpack0 {
	Uint name;
	Uint kind;
	Uint nkw;
	Uint alt;
	Uint req;
	Uint kw;
	Uint Lb;
	Uint Ub;
	int poff;
	} Kwpack0;	/* version to use when pad field of KwpHead is zero */

 typedef struct
NIDR_KWlib {
	struct NIDR_KWlib *next;
	char *libname;
	void *h;	/* library handle */
	KeyWord *kw0;	/* replaced keyword */
	KeyWord *oldtop;
	KeyWord kw;	/* replacement */
	} NIDR_KWlib;

#ifdef __cplusplus
#define externC extern "C"
#else
#define externC extern
#endif

externC NIDR_KWlib* nidr_lib_record(void *h, const char *name);
externC NIDR_KWlib *NIDR_Libs;
#endif
