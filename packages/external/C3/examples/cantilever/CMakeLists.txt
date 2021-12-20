cmake_minimum_required (VERSION 2.6)


add_executable(c3_cantilever cantilever.c)
target_link_libraries(c3_cantilever -g
  c3
  ${blas_lapack}
  -lm)



