# Users building code depending on NIDR only need libnidr.a and its
# (generated) dependencies are checked into the repository.  We'll
# maintain this convention for now.
#
# Developers of NIDR clients (e.g., DAKOTA), will routinely need
# nidrgen (and possibly Dt) when their input specification changes.
# nidrgen has a dependency on lex/sed when nidrgen.l changes
#
# Developers of NIDR proper will need to use lex, yacc, sed to
# regenerate source files and check them into the repository under
# nidr/.  (make src-scanner)

# To build nidrgen for Jaguar (static if possible)
#   configure with --enable-spec-maint
#   rm nidrgen.o nidrgen
#     linux64: make AM_CFLAGS=-m64 AM_LDFLAGS="-static -s" nidrgen
#     linux32: make AM_CFLAGS=-m32 AM_LDFLAGS="-static -s" nidrgen
#     darwin:  make AM_LDFLAGS=-s nidrgen
#     mingw:   need to make w/o -ldl


# TODO:
# * Only need lex, yacc, sed if regenerating from new sources (test
#   conditionally)
# * Review all rules in makefile.orig
# * Migrate rules for static builds to here
# * Migrate dylib tests here
# * Offer option to only build nidr-scanner.c nidr-parser.c on demand
#   (and therefore have to check into repo)?  Seems could eliminate,
#   but should check with Dave
#
# * Put headers in nidr subdir at least when installing


# Always build libnidr, possibly from sources checked into the repo
lib_LIBRARIES = libnidr.a
libnidr_a_SOURCES = nidr-scanner.c nidr-parser.c avltree.c nidr.c

## The sense of NIDR's DYNLIB flag is opposite: if no DL_SOLVER, use the 
## pre-processor to turn off NIDR's dynamic loading capability
## For now, only libnidr respects this define, not nidrgen
if !WITH_DL_SOLVER
##CPPFLAGS += -DNO_NIDR_DYNLIB
libnidr_a_CPPFLAGS = -DNO_NIDR_DYNLIB $(AM_CPPFLAGS)
endif

# For now we assume specification maintainers live on a Linux or
# similar platform, which has lex, yacc, requires libdl
# TODO: conditionally probe for these
if ENABLE_SPEC_MAINT
  noinst_PROGRAMS = nidrgen Dt
  # NOTE: Dave's rule was _conditionally_ adding -ldl only if needed
  # (otherwise will break on Windows)
  nidrgen_LDADD = -lm -ldl


### Avoid regenerating files needlessly, which can be counterproductive on
### systems with older versions of [f]lex and bison:  require explicit
### "make src-scanner" invocation to update $S/nidr-scanner and $S/nidr-parser.c
### and explicit "make nidrgen.c" to recreate nidrgen.c from nidrgen.l.
src-scanner: nidr-scanner.c nidr-parser.h nidr-parser.c


# Using custom rules due to sed
# Using separate rules for .c and .h to manage deps, despite running YACC twice
nidr-parser.h: nidrparse.y
	$(YACC) -d nidrparse.y
	$(SED) 's/yy/nidr/' y.tab.h > nidr-parser.h
	$(RM) rm y.tab.c y.tab.h

nidr-parser.c: nidrparse.y
	$(YACC) -d nidrparse.y
	$(SED) 's/y.tab.c/nidr-parser.c/;s!#line.*!/* & */!;s/yy/nidr/g;s/reset_lex_state/reset_nidrlex_state/' y.tab.c > nidr-parser.c
	$(RM) rm y.tab.c y.tab.h

nidr-scanner.c: nidrscan.l
	$(LEX) nidrscan.l
	$(SED) 's/lex.yy.c/nidr-scanner.c/;s!#line.*!/* & */!;s/yy/nidr/g;s/reset_lex_state/reset_nidrlex_state/' lex.yy.c > nidr-scanner.c
	$(RM) lex.yy.c

nidrgen.c: nidrgen.l
	$(LEX) nidrgen.l
	$(SED) 's/lex.yy.c/nidrgen.c/;s!#line.*!/* & */!;s/yyparse/nidrparse/g' lex.yy.c > nidrgen.c
	$(RM) lex.yy.c

endif


# Headers to install
include_HEADERS = nidr.h
