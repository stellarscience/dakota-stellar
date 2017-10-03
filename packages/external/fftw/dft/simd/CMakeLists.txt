INCLUDE_DIRECTORIES(
 ${FFTW_SOURCE_DIR}/kernel
 ${FFTW_SOURCE_DIR}/dft
 ${FFTW_SOURCE_DIR}/simd
)

ADD_SUBDIRECTORY(codelets)

ADD_LIBRARY(dft_simd 
  n1b.c n1f.c n2b.c n2f.c n2s.c q1b.c q1f.c t.c ts.c n1b.h n1f.h n2b.h n2f.h
  n2s.h q1b.h q1f.h t1b.h t1f.h t2b.h t2f.h t3b.h t3f.h ts.h t1fu.h t1bu.h
)

