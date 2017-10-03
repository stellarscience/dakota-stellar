#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set list of of tests to be built run by 'make check' and provide the
# relevant source files.

TESTS = tstillum tstrosen tstlinear tstnonlinear
check_PROGRAMS = $(TESTS)

tstillum_SOURCES = tstillum.C illum.C MAIN.c tstfcn.h
tstrosen_SOURCES = tstrosen.C rosen.C MAIN.c tstfcn.h
tstlinear_SOURCES = tstlinear.C rosen.C MAIN.c tstfcn.h
tstnonlinear_SOURCES = tstnonlinear.C sixtyfive.C MAIN.c tstfcn.h

TEUCHOS_INCL	= @TEUCHOS_CPPFLAGS@
TEUCHOS_ARCHIVE	= @TEUCHOS_LDFLAGS@ -lteuchos

# Provide location of additional include files.

INCLUDES = $(TEUCHOS_INCL) -I$(top_srcdir)/include

# Provide libraries to be linked in.

OPTPP_DEPS = $(top_builddir)/lib/liboptpp.la \
	     $(TEUCHOS_ARCHIVE) \
	     $(NPSOL_LIB) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

tstillum_LDADD = $(OPTPP_DEPS)
tstrosen_LDADD = $(OPTPP_DEPS)
tstlinear_LDADD = $(OPTPP_DEPS)
tstnonlinear_LDADD = $(OPTPP_DEPS)

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL

# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *.ti *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -rf ii_files SunWS_cache so_locations
