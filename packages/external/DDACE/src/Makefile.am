## Process this file with automake to produce 'Makefile.in'.

## Source subdirectories to be recursed into during the build.
## Convenience libraries to be built in each of the SUBDIRS.
## Add FATools iff ENABLE_STANDALONE has been set!
SUBDIRS = Analyzer Bose
SUBLIBS = Analyzer/libanalyzer.la Bose/libbose.la
if ENABLE_STANDALONE
SUBDIRS += FATools
SUBLIBS += FATools/libfatools.la
endif
SUBDIRS += Random Sampling
SUBLIBS += Random/librandom.la Sampling/libsampling.la

## Set the DDACE library name.
lib_LTLIBRARIES = libddace.la

## No sources in this directory!
libddace_la_SOURCES = 

## Combine the convenience libraries into the final product.
libddace_la_LIBADD = $(SUBLIBS)

## Remove the Makefile.in when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = COPYRIGHT GNU_LGPL
