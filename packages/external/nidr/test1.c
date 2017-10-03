#include <stdio.h>
#include "nidr.h"

char nidr_please_refer[] =
	"Please refer to the DAKOTA Reference and/or User's Manual or the\n" 
	"dakota.input.summary file distributed with this executable.";

 static void
showvals(Values *v)
{
	Real *r;
	const char **s;
	int *i;
	size_t j, n;

	if (v && (n = v->n) > 0) {
		printf("%u values:\n", n);
		if (i = v->i)
			for(j = 0; j < n; j++)
				printf("\ti[%d] = %d\n", j, i[j]);
		if (r = v->r)
			for(j = 0; j < n; j++)
				printf("\tr[%d] = %g\n", j, r[j]);
		if (s = v->s)
			for(j = 0; j < n; j++)
				printf("s[%d] = \"%s\"\n", j, s[j]);
		}
	}

 static void
kwstart(const char *name, Values *vals, void **gp, void *v)
{
	printf("Start of \"%s\": *gp = %d, v = #%x\n", name, (*(int*)gp)++, v);
	showvals(vals);
	}

 static void
kwend(const char *name, Values *vals, void **gp, void *v)
{
	printf("End of \"%s\": *gp = %d, v = #%x\n", name, (*(int*)gp)++, v);
	showvals(vals);
	}

 static void
dummystart(const char *name, Values *vals, void **gp, void *v)
{}

#ifdef INCLUDENAME
#include INCLUDENAME
#else
#define V(x) (void*)x

 static KeyWord
prog_thresh = { "progress_threshold", KWKind_Real, 0,0,1,0, kwstart, V(V(0x1)), kwend, V(V(0x2))};

 static KeyWord
unc_kw[3] = {
	{ "adaptive_hybrid", KWKind_Void, 1,0,0,&prog_thresh,kwstart,V(0x3),kwend,V(0x4)},
	{ "method_list", KWKind_Str|KWKind_List, 0,0,1,0,kwstart,V(0x7),kwend,V(0x8)},
	{ "num_solutions_transferred", KWKind_Int,0,0,0,0,kwstart,V(0x5),kwend,V(0x6)}
	},
cpl_kw[3] = {
	{ "global_method_pointer", KWKind_Str, 0,0,1,0,kwstart,V(0xb),kwend,V(0xc)},
	{ "local_method_pointer", KWKind_Str, 0,0,2,0,kwstart,V(0xd),kwend,V(0xe)},
	{ "local_search_probability", KWKind_Real,0,0,0,0,kwstart,V(0x20),kwend,V(0x21)}
	};

 static KeyWord
ml_kw[2] = {
	{ "coupled", KWKind_Void, 3,1,1,cpl_kw,kwstart,V(0x22),kwend,V(0x23)},
	{ "uncoupled", KWKind_Void, 3,1,1,unc_kw,kwstart,V(0x9),kwend,V(0xa)}
	},
sbokw[2] = {
	{ "max_iterations", KWKind_Int, 0,0,0,0, kwstart,V(0x28),kwend,V(0x29)},
	{ "opt_method_pointer", KWKind_Str, 0,0,1,0, kwstart,V(0x26),kwend,V(0x27)}
	},
tab_gr_file = { "tabular_graphics_file", KWKind_Str, 0,0,0,0, kwstart, V(0x123), kwend, V(0x456) },
stkw[4] = {
	{ "graphics", KWKind_Void, 0, 0, 0, 0, kwstart, V(0x33), kwend, V(0x34) },
	{ "multi_level", KWKind_Void, 2,1,1,ml_kw,kwstart,V(0x24),kwend,V(0x25)},
	{ "surrogate_based_opt", KWKind_Void, 2,1,1,sbokw,kwstart,V(0x28),kwend,V(0x29)},
	{ "tabular_graphics_data", KWKind_Void, 1,0,0,&tab_gr_file, kwstart, V(0x23), kwend,V(0x45)}
	},
methkw[2] = {
	{ "id_method", KWKind_Str, 0,0,0,0,kwstart,V(0x2c),kwend,V(0x2d)},
	{ "model_pointer", KWKind_Str, 0,0,0,0,kwstart,V(0x2e),kwend,V(0x2f)}
	},
top[2] = {
	{ "method", KWKind_Void, 2,0,0,methkw,kwstart,V(0x30),kwend,V(0x31)},
	{ "strategy", KWKind_Void, 4,0,0,stkw,kwstart,V(0x2a),kwend,V(0x2b)}
	};

KeyWord Dakota_Keyword_Top = {"KeywordTop",KWKind_Void,2,0,0,top,kwstart,V(0x987),kwend,V(0x9871)};
#endif

extern int nidrparse(void);
extern void nidr_setup(const char*);
extern FILE *nidrin;

 int
main(int argc, char **argv)
{
	char *progname, *s;

	progname = *argv;
	if ((s = argv[1]) && *s == '-') {
		nidr_setup(s+1);
		++argv;
		--argc;
		}
	if (argc > 1 && !(nidrin = fopen(argv[1],"r"))) {
		printf("%s: Could not open \"%s\"\n", progname, argv[1]);
		return 1;
		}
	printf("\n\nnidrparse returned %d\n", nidrparse());
	return 0;
	}

 void
abort_handler(int n) {}
