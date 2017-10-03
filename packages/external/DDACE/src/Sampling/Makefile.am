## Process this file with automake to produce 'Makefile.in'.

## Convenience library for the Sampling subdirectory.
noinst_LTLIBRARIES = libsampling.la

## SOURCE dependencies for this convenience library.
libsampling_la_SOURCES = DDaceArraySampler.cpp \
			 DDaceBoxBehnkenSampler.cpp \
			 DDaceCentralCompositeSampler.cpp \
			 DDaceFactorialSampler.cpp \
			 DDaceLHSampler.cpp \
			 DDaceOALHSampler.cpp \
			 DDaceOASampler.cpp \
			 DDaceRandomSampler.cpp \
			 DDaceSamplePoint.cpp \
			 DDaceSampler.cpp \
			 DDaceUserInputSampler.cpp

## Include directory for DDACE headers.
INCLUDES = -I$(top_srcdir)/include

## Remove the Makefile.in when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = COPYRIGHT GNU_LGPL
