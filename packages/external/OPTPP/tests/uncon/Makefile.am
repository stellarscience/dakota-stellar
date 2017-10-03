#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set list of of tests to be built run by 'make check' and provide the
# relevant source files.

TESTS = tstqnewton tstgnewton tstnewton tstfdnlf1 tstcg tstLBFGS
check_PROGRAMS = $(TESTS)

tstqnewton_SOURCES = tstqnewton.C rosen.C tstfcn.h
tstgnewton_SOURCES = tstgnewton.C rosen.C tstfcn.h
tstnewton_SOURCES = tstnewton.C rosen.C tstfcn.h
tstfdnlf1_SOURCES = tstfdnlf1.C rosen.C tstfcn.h
tstcg_SOURCES = tstcg.C rosen.C tstfcn.h
tstLBFGS_SOURCES = tstLBFGS.C rosen.C tstfcn.h

TEUCHOS_INCL	= @TEUCHOS_CPPFLAGS@
TEUCHOS_ARCHIVE	= @TEUCHOS_LDFLAGS@ -lteuchos

# Provide location of additional include files.

INCLUDES = $(TEUCHOS_INCL) -I$(top_srcdir)/include

# Provide libraries to be linked in.

OPTPP_DEPS = $(top_builddir)/lib/liboptpp.la \
	     $(TEUCHOS_ARCHIVE) \
	     $(NPSOL_LIB) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

tstqnewton_LDADD = $(OPTPP_DEPS)
tstgnewton_LDADD = $(OPTPP_DEPS)
tstnewton_LDADD = $(OPTPP_DEPS)
tstfdnlf1_LDADD = $(OPTPP_DEPS)
tstcg_LDADD = $(OPTPP_DEPS)
tstLBFGS_LDADD = $(OPTPP_DEPS)

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL

# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *.ti *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -rf ii_files SunWS_cache so_locations
