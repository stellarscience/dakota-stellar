cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

project(NCSUOpt Fortran)

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget NCSUOptTargets CACHE STRING 
  "Export set name for ${CMAKE_PROJECT_NAME}")

add_library(ncsuopt DIRect.f DIRbatch.f DIRsubrout.f)

install(
  TARGETS ncsuopt
  EXPORT  ${ExportTarget}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  )

# Install package configuration and import files.
configure_file(${NCSUOpt_SOURCE_DIR}/cmake/NCSUOptConfig.cmake.in
               ${NCSUOpt_BINARY_DIR}/cmake/NCSUOptConfig.cmake @ONLY)
install(FILES ${NCSUOpt_BINARY_DIR}/cmake/NCSUOptConfig.cmake
        DESTINATION lib/cmake/NCSUOpt)

# Workaround until we fix Dakota's overall use of exports
if(ExportTarget STREQUAL NCSUOptTargets)
  message(STATUS 
    "Installing NCSUOpt exported targets to lib/cmake/NCSUOpt")
  install(EXPORT ${ExportTarget} DESTINATION lib/cmake/NCSUOpt)
endif()
