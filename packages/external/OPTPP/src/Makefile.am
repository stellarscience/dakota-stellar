#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Source subdirectories to be included during the build and
# convenience libraries to be built in each subdirectory.

SUBDIRS = Base Constraints GSS Newton PDS Utils
SUBLIBS = Base/libbase.la	Constraints/libconstraints.la	\
	  GSS/libgss.la		Newton/libnewton.la		\
	  PDS/libpds.la		Utils/libutils.la
DIST_SUBDIRS = Base Constraints GSS Newton PDS Utils

# Set main library name, include source files in this directory,
# and combine convenience libraries into main library.

lib_LTLIBRARIES = @top_builddir@/lib/liboptpp.la
@top_builddir@_lib_liboptpp_la_SOURCES = globals.C
@top_builddir@_lib_liboptpp_la_LIBADD = $(SUBLIBS)

# Provide location of additional include files.

TEUCHOS_INCL	= @TEUCHOS_CPPFLAGS@
INCLUDES = $(TEUCHOS_INCL) -I$(top_srcdir)/include

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL

# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *.ti *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -rf ii_files SunWS_cache so_locations
