INCLUDE_DIRECTORIES(
  ${FFTW_SOURCE_DIR}/kernel
  ${FFTW_SOURCE_DIR}/dft
  ${FFTW_BINARY_DIR}
)

ADD_SUBDIRECTORY(codelets)

ADD_LIBRARY(dft_scalar n.c t.c f.h n.h q.h t.h)

