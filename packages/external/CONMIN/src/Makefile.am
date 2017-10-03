## Process this file with automake to produce 'Makefile.in'.

## Set the CONMIN library name.
lib_LTLIBRARIES = libconmin.la

## CONMIN library source.
libconmin_la_SOURCES = conmin.f

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations
