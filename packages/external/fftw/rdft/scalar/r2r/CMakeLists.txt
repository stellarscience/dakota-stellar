INCLUDE_DIRECTORIES(
  ${FFTW_SOURCE_DIR}/kernel
  ${FFTW_SOURCE_DIR}/rdft
  ${FFTW_SOURCE_DIR}/rdft/scalar
)

ADD_LIBRARY(rdft_scalar_r2r e01_8.c e10_8.c codlist.c)
TARGET_LINK_LIBRARIES(rdft_scalar_r2r rdft_scalar)

