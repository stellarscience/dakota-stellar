## Process this file with automake to produce 'Makefile.in'.

## Convenience library for the Bose subdirectory.
noinst_LTLIBRARIES = libbose.la

## SOURCE dependencies for this convenience library.
libbose_la_SOURCES = akconst.c \
		     boselink.c \
		     construct.c \
		     galois.c \
		     gfields.c \
		     memory.c \
		     oa.c \
		     primes.c

## Include directory for DDACE headers.
INCLUDES = -I$(top_srcdir)/include

## Remove the Makefile.in when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = COPYRIGHT GNU_LGPL LICENSE_BOSE
