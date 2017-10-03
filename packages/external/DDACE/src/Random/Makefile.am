## Process this file with automake to produce 'Makefile.in'.

## Convenience library for the Random subdirectory.
noinst_LTLIBRARIES = librandom.la

## SOURCE dependencies for this convenience library.
librandom_la_SOURCES = dcdflib.c \
		       Distribution.cpp \
		       dlaran.f \
		       ipmpar.c \
		       NormalDistribution.cpp \
		       PseudoRandomTestsOnly.cpp \
		       UniformDistribution.cpp

## Include directory for DDACE headers.
INCLUDES = -I$(top_srcdir)/include

## Remove the Makefile.in when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = COPYRIGHT GNU_LGPL
