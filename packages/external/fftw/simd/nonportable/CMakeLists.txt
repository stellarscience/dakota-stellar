# this directory contains files that can only be compiled with
# SIMD_CFLAGS

include_directories(${PROJECT_SOURCE_DIR}/packages/fftw/kernel
                    ${PROJECT_SOURCE_DIR}/packages/fftw/simd
                    ${PROJECT_BINARY_DIR}/packages/fftw
)
add_library(simd_nonportable sse.c sse2.c)

