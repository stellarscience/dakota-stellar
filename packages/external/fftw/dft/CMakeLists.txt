INCLUDE_DIRECTORIES(${FFTW_SOURCE_DIR}/kernel
                    ${FFTW_BINARY_DIR})

IF(HAVE_SIMD)
  ADD_SUBDIRECTORY(simd)
ENDIF(HAVE_SIMD)

ADD_SUBDIRECTORY(scalar)

SET(DFT_SRCs 
  bluestein.c buffered.c conf.c ct.c dftw-direct.c	
  dftw-directsq.c dftw-generic.c dftw-genericbuf.c direct.c generic.c	
  indirect.c indirect-transpose.c kdft-dif.c kdft-difsq.c kdft-dit.c	
  kdft.c nop.c plan.c problem.c rader.c rank-geq2.c solve.c vrank-geq1.c	
  zero.c codelet-dft.h ct.h dft.h
)
ADD_LIBRARY(dft ${DFT_SRCs})
TARGET_LINK_LIBRARIES(dft kernel dft_scalar_codelets)

#TARGET_LINK_LIBRARIES(dft kernel)# dft_standard)
