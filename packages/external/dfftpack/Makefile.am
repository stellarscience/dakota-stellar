## Process this file with automake to produce 'Makefile.in'.

## Set the dfftpack library name.
noinst_LTLIBRARIES = libdfftpack.la

## dfftpack library sources.
libdfftpack_la_SOURCES = \
zfftb.f		\
cfftb1.f	\
zfftf.f		\
cfftf1.f	\
zffti.f		\
cffti1.f	\
dcosqb.f	\
cosqb1.f	\
dcosqf.f	\
cosqf1.f	\
dcosqi.f	\
dcost.f		\
dcosti.f	\
ezfft1.f	\
dzfftb.f	\
dzfftf.f	\
dzffti.f	\
passb.f		\
passb2.f	\
passb3.f	\
passb4.f	\
passb5.f	\
passf.f		\
passf2.f	\
passf3.f	\
passf4.f	\
passf5.f	\
radb2.f		\
radb3.f		\
radb4.f		\
radb5.f		\
radbg.f		\
radf2.f		\
radf3.f		\
radf4.f		\
radf5.f		\
radfg.f		\
dfftb.f		\
rfftb1.f	\
dfftf.f		\
rfftf1.f	\
dffti.f		\
rffti1.f	\
dsinqb.f	\
dsinqf.f	\
dsinqi.f	\
dsint.f		\
sint1.f		\
dsinti.f

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations
