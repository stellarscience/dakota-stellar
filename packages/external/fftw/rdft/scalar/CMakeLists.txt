INCLUDE_DIRECTORIES(${FFTW_SOURCE_DIR} ${FFTW_SOURCE_DIR}/rdft)

ADD_SUBDIRECTORY(r2cf)
ADD_SUBDIRECTORY(r2cb)
ADD_SUBDIRECTORY(r2r)

ADD_LIBRARY(rdft_scalar 
 hb.h r2cb.h r2cbIII.h hf.h hfb.c r2c.c r2cf.h r2cfII.h r2r.c r2r.h hc2c.c
 hc2cf.h hc2cb.h
) 
TARGET_LINK_LIBRARIES(rdft_scalar 
  rdft_scalar_r2cb
  rdft_scalar_r2cf 
  rdft_scalar_r2r
)

