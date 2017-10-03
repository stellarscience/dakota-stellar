set(Utils_sources
  BoolVector.C file_cutils.c ioformat.C mcholesky.C OptppExceptions.C
  OptppFatalError.C print.C timers.c
)

if(NOT BLAS_LIBS)
  #list(APPEND libutils_srcs linalg.c)
  message(FATAL_ERROR "  CMake cannot build Opt++ without
    BLAS. Please provide full path to libblas or equivalent.")
endif()
