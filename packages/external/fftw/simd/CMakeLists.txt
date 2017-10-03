ADD_SUBDIRECTORY(nonportable)
INCLUDE_DIRECTORIES(${FFTW_SOURCE_DIR}/kernel)

ADD_LIBRARY(simd
  altivec.c sse.c sse2.c simd-altivec.h simd-sse.h simd-sse2.h simd.h taint.c
  x86-cpuid.h mips_ps.c mips_ps.h simd-mips_ps.h
)

