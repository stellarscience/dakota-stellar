#AM_CPPFLAGS = -I$(top_srcdir)/kernel -I$(top_srcdir)/dft	\
#-I$(top_srcdir)/rdft -I$(top_srcdir)/api
#
#if THREADS
#if COMBINED_THREADS
#noinst_LTLIBRARIES = libfftw3@PREC_SUFFIX@_threads.la
#else
#lib_LTLIBRARIES = libfftw3@PREC_SUFFIX@_threads.la
#endif
#endif

include_directories(${FFTW_SOURCE_DIR}/api
                    ${FFTW_SOURCE_DIR}/kernel
                    ${FFTW_SOURCE_DIR}/dft
                    ${FFTW_SOURCE_DIR}/rdft
                    ${FFTW_BINARY_DIR})

add_library(fftw3_threads 
  api.c conf.c threads.c threads.h openmp.c dft-vrank-geq1.c ct.c
  rdft-vrank-geq1.c hc2hc.c vrank-geq1-rdft2.c f77api.c f77funcs.h
)

