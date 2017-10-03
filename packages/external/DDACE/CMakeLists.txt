cmake_minimum_required(VERSION 2.8)
project(ddace C CXX Fortran)

include(CheckFunctionExists)

remove_definitions("-DHAVE_CONFIG_H")

include(FortranCInterface)
FortranCInterface_HEADER(ddace_fortran.h MACRO_NAMESPACE DDACE_FORTRAN_)

include_directories(${ddace_SOURCE_DIR}/include ${ddace_BINARY_DIR})

option(DDACE_ENABLE_STANDALONE "Enable DDACE standalone?" OFF)

if(Boost_DIR)
  include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
  add_definitions(-DHAVE_BOOST)
else()
  find_package(Boost)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
    add_definitions(-DHAVE_BOOST)
  endif()
endif()

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ddaceTargets CACHE STRING 
  "Export set name for ${CMAKE_PROJECT_NAME}")

macro(ddace_install_library lib)
  install(TARGETS ${lib} EXPORT ${ExportTarget}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    )
endmacro()

add_subdirectory(src)

option(DDACE_ENABLE_TESTS "Enable DDACE tests?" ON)
if(DDACE_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

# Install package configuration and import files.
configure_file(${ddace_SOURCE_DIR}/cmake/ddace-config.cmake.in
               ${ddace_BINARY_DIR}/config/cmake/ddace-config.cmake @ONLY)
install(FILES ${ddace_BINARY_DIR}/config/cmake/ddace-config.cmake
        DESTINATION lib/cmake/ddace)

# Workaround until we fix Dakota's overall use of exports
if(ExportTarget STREQUAL ddaceTargets)
  message(STATUS 
    "Installing ddace exported targets to lib/cmake/ddace")
  install(EXPORT ${ExportTarget} DESTINATION lib/cmake/ddace)
endif()

# TODO: move to directory-protected headers
#install(DIRECTORY include/ DESTINATION include/ddace PATTERN .svn EXCLUDE)
#install(FILES ${ddace_BINARY_DIR}/ddace_fortran.h
#        DESTINATION include/ddace)
file (GLOB ddace_headers include/*.h)
install(FILES ${ddace_headers} DESTINATION include)
install(FILES ${ddace_BINARY_DIR}/ddace_fortran.h DESTINATION include)
