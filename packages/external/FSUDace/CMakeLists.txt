cmake_minimum_required(VERSION 2.8)
project(FSUDace)

add_library(fsudace
  fsu_halton.C fsu_hammersley.C fsu_cvt.C
  fsu_latinize.C fsu_quality.C fsu_sub.C
)

# Install standalone binaries, but not testers
set(fsu_standalone
  fsu_halton_standalone
  fsu_hammersley_standalone
  fsu_cvt_standalone
  fsu_latinize_standalone
  fsu_quality_standalone
  )

set(fsu_test_problems
  fsu_cvt_problems
  fsu_latinize_problems
  fsu_quality_problems
  )

foreach(exe ${fsu_standalone} ${fsu_test_problems})
  add_executable(${exe} ${exe}.C)
  target_link_libraries(${exe} fsudace)
endforeach()

# Install headers.
install(
  FILES fsu.H
  DESTINATION include
  )

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget FSUDaceTargets CACHE STRING 
  "Export set name for ${CMAKE_PROJECT_NAME}")

# Install binaries.
install(
  TARGETS fsudace ${fsu_standalone}
  EXPORT  ${ExportTarget}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  )

# Install package configuration and import files.
configure_file(${FSUDace_SOURCE_DIR}/cmake/FSUDaceConfig.cmake.in
               ${FSUDace_BINARY_DIR}/cmake/FSUDaceConfig.cmake @ONLY)
install(FILES ${FSUDace_BINARY_DIR}/cmake/FSUDaceConfig.cmake
        DESTINATION lib/cmake/FSUDace)

# Workaround until we fix Dakota's overall use of exports
if(ExportTarget STREQUAL FSUDaceTargets)
  message(STATUS 
    "Installing FSUDace exported targets to lib/cmake/FSUDace")
  install(EXPORT ${ExportTarget} DESTINATION lib/cmake/FSUDace)
endif()
