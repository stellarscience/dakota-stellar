## Process this file with automake to produce 'Makefile.in'.

## Need the m4 directory for non-native macro tests.
ACLOCAL_AMFLAGS = -I ../../config

## Set the FSUDace library name
lib_LTLIBRARIES = libfsudace.la

## FSUDace library sources.
libfsudace_la_SOURCES = fsu_halton.C fsu_hammersley.C fsu_cvt.C \
                        fsu_latinize.C fsu_quality.C fsu_sub.C

## FSUDace header which needs to be distributed.
include_HEADERS = fsu.H

## FSUDace standalone binaries.
bin_PROGRAMS = fsu_halton_standalone            \
               fsu_hammersley_standalone        \
               fsu_cvt_standalone               \
               fsu_cvt_problems                 \
               fsu_latinize_standalone          \
               fsu_latinize_problems            \
               fsu_quality_standalone           \
               fsu_quality_problems

## Sources and library link for same.
fsu_halton_standalone_SOURCES = fsu_halton_standalone.C
fsu_halton_standalone_LDADD = $(top_builddir)/libfsudace.la
fsu_hammersley_standalone_SOURCES = fsu_hammersley_standalone.C
fsu_hammersley_standalone_LDADD = $(top_builddir)/libfsudace.la
fsu_cvt_standalone_SOURCES = fsu_cvt_standalone.C
fsu_cvt_standalone_LDADD = $(top_builddir)/libfsudace.la
fsu_cvt_problems_SOURCES = fsu_cvt_problems.C
fsu_cvt_problems_LDADD = $(top_builddir)/libfsudace.la
fsu_latinize_standalone_SOURCES = fsu_latinize_standalone.C
fsu_latinize_standalone_LDADD = $(top_builddir)/libfsudace.la
fsu_latinize_problems_SOURCES = fsu_latinize_problems.C
fsu_latinize_problems_LDADD = $(top_builddir)/libfsudace.la
fsu_quality_standalone_SOURCES = fsu_quality_standalone.C
fsu_quality_standalone_LDADD = $(top_builddir)/libfsudace.la
fsu_quality_problems_SOURCES = fsu_quality_problems.C
fsu_quality_problems_LDADD = $(top_builddir)/libfsudace.la

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 config.guess config.sub \
		       configure depcomp fsudace_config.h.in install-sh \
		       ltmain.sh missing

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = GNU_LGPL                           \
             examples/cvt_07_00010.latin.txt    \
             examples/cvt_07_00010_qual.txt     \
             examples/cvt_07_00010.txt          \
             examples/fsu_cvt_problems.out      \
             examples/fsu_latinize_problems.out \
             examples/fsu_quality_problems.out  \
             examples/halton_02_00100.latin.txt \
             examples/halton_02_00100_qual.txt  \
             examples/halton_02_00100.txt
