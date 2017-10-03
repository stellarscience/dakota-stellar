include_directories(${FFTW_SOURCE_DIR}/libbench2
                    ${FFTW_SOURCE_DIR}/api
                    ${FFTW_SOURCE_DIR}/dft
                    ${FFTW_SOURCE_DIR}/rdft
                    ${FFTW_SOURCE_DIR}/kernel
                    ${FFTW_BINARY_DIR})

add_executable(bench bench.c hook.c fftw-bench.c fftw-bench.h)
#IF(THREADS and COMBINED_THREADS)
target_link_libraries(bench
  bench2 fftw3_threads fftw3 ${CMAKE_THREAD_LIBS_INIT})

#CMake doesn't support linking against object files
#adding this new library so that fftw-wisdom can link properly
#ADD_LIBRARY(bench1 bench.c)
add_library(bench1 bench.c hook.c fftw-bench.c fftw-bench.h)
target_link_libraries(bench1
  fftw3_threads fftw3 bench2 ${CMAKE_THREAD_LIBS_INIT})

