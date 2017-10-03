include_directories(${FFTW_SOURCE_DIR}/libbench2
                    ${FFTW_SOURCE_DIR}/api
                    ${FFTW_BINARY_DIR})
add_executable(fftw_wisdom fftw-wisdom.c)
target_link_libraries(fftw_wisdom bench1 bench2)

