#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set list of of tests to be built run by 'make check' and provide the
# relevant source files.

TESTS = tsthock1 tsthock2 tsthock5 tsthock6 tsthock7 tsthock10	    \
	tsthock13 tsthock14 tsthock26 tsthock28 tsthock35 tsthock65 \
	tsthock77 tsthock78
check_PROGRAMS = $(TESTS)

tsthock1_SOURCES = tsthock1.C hockfcns.C hockfcns.h
tsthock2_SOURCES = tsthock2.C hockfcns.C hockfcns.h
tsthock5_SOURCES = tsthock5.C hockfcns.C hockfcns.h
tsthock6_SOURCES = tsthock6.C hockfcns.C hockfcns.h
tsthock7_SOURCES = tsthock7.C hockfcns.C hockfcns.h
tsthock10_SOURCES = tsthock10.C hockfcns.C hockfcns.h
tsthock13_SOURCES = tsthock13.C hockfcns.C hockfcns.h
tsthock14_SOURCES = tsthock14.C hockfcns.C hockfcns.h
tsthock26_SOURCES = tsthock26.C hockfcns.C hockfcns.h
tsthock28_SOURCES = tsthock28.C hockfcns.C hockfcns.h
tsthock35_SOURCES = tsthock35.C hockfcns.C hockfcns.h
tsthock65_SOURCES = tsthock65.C hockfcns.C hockfcns.h
tsthock77_SOURCES = tsthock77.C hockfcns.C hockfcns.h
tsthock78_SOURCES = tsthock78.C hockfcns.C hockfcns.h

TEUCHOS_INCL	= @TEUCHOS_CPPFLAGS@
TEUCHOS_ARCHIVE	= @TEUCHOS_LDFLAGS@ -lteuchos

# Provide location of additional include files.

INCLUDES = $(TEUCHOS_INCL) -I$(top_srcdir)/include

# Provide libraries to be linked in.

OPTPP_DEPS = $(top_builddir)/lib/liboptpp.la \
	     $(TEUCHOS_ARCHIVE) \
	     $(NPSOL_LIB) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

tsthock1_LDADD = $(OPTPP_DEPS)
tsthock2_LDADD = $(OPTPP_DEPS)
tsthock5_LDADD = $(OPTPP_DEPS)
tsthock6_LDADD = $(OPTPP_DEPS)
tsthock7_LDADD = $(OPTPP_DEPS)
tsthock10_LDADD = $(OPTPP_DEPS)
tsthock13_LDADD = $(OPTPP_DEPS)
tsthock14_LDADD = $(OPTPP_DEPS)
tsthock26_LDADD = $(OPTPP_DEPS)
tsthock28_LDADD = $(OPTPP_DEPS)
tsthock35_LDADD = $(OPTPP_DEPS)
tsthock65_LDADD = $(OPTPP_DEPS)
tsthock77_LDADD = $(OPTPP_DEPS)
tsthock78_LDADD = $(OPTPP_DEPS)

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL

# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *.ti *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -rf ii_files SunWS_cache so_locations
