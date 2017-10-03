INCLUDE_DIRECTORIES(
  ${FFTW_SOURCE_DIR}/kernel
  ${FFTW_SOURCE_DIR}/dft
  ${FFTW_SOURCE_DIR}/simd
  ${FFTW_BINARY_DIR})

#if(HAVE_CELL) #not sure when this should be defined...
#  ADD_SUBDIRECTORY(spu)
#endif(HAVE_CELL)

if(SINGLE)
  SET(PLANS plans-single.c)
else(SINGLE)
  SET(PLANS plans-double.c)
endif(SINGLE)

ADD_LIBRARY(cell 
  cell.c conf.c copy.c dft-direct-cell.c
  transpose.c ${PLANS} fftw-cell.h)
