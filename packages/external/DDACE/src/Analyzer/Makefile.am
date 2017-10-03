## Process this file with automake to produce 'Makefile.in'.

## NOTE: There are source files in this directory which aren't currently
##       being compiled that depend on MARS in ../FATools; they are:
##
##       ddaceAnalysis.cpp
##       ResponseSurfaceAnalyzer.cpp
##
##       It will be necessary to fence the MARS-dependent sections of these
##       sources with the the preprocessor conditional HAVE_MARS if they are
##       ever re-enabled.

## Convenience library for the Analyzer subdirectory.
noinst_LTLIBRARIES = libanalyzer.la

## SOURCE dependencies for this convenience library.
libanalyzer_la_SOURCES = ColumnHeader.cpp \
                         DataValue.cpp \
                         Factor.cpp \
                         MainEffectsAnalyzer3.cpp \
			 MainEffectsConverter.cpp \
			 MainEffectsExcelOutput.cpp \
			 OneWayANOVA.cpp \
                         Response.cpp \
                         Statistics.cpp \
                         ValueAndRowIndexAndColumnIndex.cpp \
                         VectorCountingNumbersAndCount.cpp

## Include directory for DDACE headers.
INCLUDES = -I$(top_srcdir)/include
if WITH_GSL
INCLUDES += $(GSL_CPPFLAGS)
endif

## Remove the Makefile.in when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = COPYRIGHT GNU_LGPL
