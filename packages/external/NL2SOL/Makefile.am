## Process this file with automake to produce 'Makefile.in'.

## Set the NL2SOL library name
lib_LTLIBRARIES = libcport.la

## NL2SOL library sources.
libcport_la_SOURCES = dn2f.c dn2g.c drn2g.c dg7lit.c df7hes.c dl7svx.c \
		      dn2lrd.c do7prd.c dc7vfn.c dl7nvr.c dl7tsq.c dn2cvp.c \
		      dn2rdp.c dn2fb.c dn2gb.c drn2gb.c dg7itb.c dparck.c \
		      divset.c i7mdcn.c i1mach.c dv7dfl.c dl7sqr.c i7pnvr.c \
		      ditsum.c stopx.c dl7msb.c dl7mst.c dd7mlp.c dl7tvm.c \
		      i7copy.c ds7dmp.c dg7qsb.c dv7vmp.c dg7qts.c dl7srt.c \
		      dl7svn.c ds7bqn.c dq7rsh.c dh2rfg.c dh2rfa.c i7shft.c \
		      dv7shf.c dl7ivm.c dl7itv.c dv7ipr.c ds7ipr.c drldst.c \
		      da7sst.c dl7vml.c ds7lup.c ds7lvm.c df7dhb.c dv7cpy.c \
		      dq7apl.c dr7tvm.c dd7upd.c dv7scp.c dq7rad.c dv2nrm.c \
		      dv7scl.c dd7tpr.c dr7mdc.c d1mach.c dv2axy.c

## Header file which needs to be distributed.
include_HEADERS = f2c.h

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = COPYING INSTALL Makefile.in aclocal.m4 config.sub \
		       config.guess configure depcomp nl2sol_config.h.in \
		       install-sh ltmain.sh missing

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Extra file which must be included in the distribution.
EXTRA_DIST = port.readme
