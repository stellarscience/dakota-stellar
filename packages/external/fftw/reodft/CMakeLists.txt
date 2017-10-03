INCLUDE_DIRECTORIES(${FFTW_SOURCE_DIR}/kernel
                    ${FFTW_SOURCE_DIR}/rdft
                    ${FFTW_BINARY_DIR}
)

ADD_LIBRARY(reodft 
  conf.c reodft.h reodft010e-r2hc.c reodft11e-radix2.c reodft11e-r2hc-odd.c
  redft00e-r2hc-pad.c rodft00e-r2hc-pad.c reodft00e-splitradix.c
)

