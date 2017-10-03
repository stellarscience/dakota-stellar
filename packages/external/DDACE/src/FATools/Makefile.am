## Process this file with automake to produce 'Makefile.in'.

## Convenience library for the FATools subdirectory.
noinst_LTLIBRARIES = libfatools.la

## SOURCE dependencies for this convenience library.
libfatools_la_SOURCES = FuncApprox.cpp \
			FuncApproxBase.cpp \
			Mars.cpp \
			mars36_fort.f

## Include directory for DDACE headers.
INCLUDES = -I$(top_srcdir)/include

## Remove the Makefile.in when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = COPYRIGHT GNU_LGPL LICENSE_MARS
