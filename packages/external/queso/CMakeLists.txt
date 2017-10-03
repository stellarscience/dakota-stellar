# Experimental QUESO CMake that builds only the sources (no tests,
# examples, etc.); for testing purposes only.  Does not properly
# manage various QUESO build/conditional compile settings.
#
cmake_minimum_required(VERSION 2.8.9)
project(Queso CXX)


# -----
# Required packages
# -----

# MPI is optional and Dakota manages MPI for now
#find_package(MPI REQUIRED)

# This protection shouldn't be necessary, but having trouble with Acro
# and non-system-provided Boost.  Probably need to change all
# invocations to Boost_FOUND instead of Boost_DIR.
if(NOT Boost_FOUND)
  find_package(Boost 1.44 REQUIRED COMPONENTS program_options)
endif()
find_package(GSL REQUIRED) # Version 1.15 or newer?


# -----
# Generated headers
# -----

# Fake out includes with empty file
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/queso)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config_queso.h.cmake.in 
  ${CMAKE_CURRENT_BINARY_DIR}/include/queso/config_queso.h COPYONLY)

# Set the defines we think are used from config_queso.h
#add_definitions("-DQUESO_HAVE_BOOST")...
# TODO: consider managing these at top-level to keep in sync with dakota/src
add_definitions("-DQUESO_HAVE_BOOST_SCOPED_PTR_HPP")
if(DAKOTA_HAVE_MPI)
  add_definitions("-DQUESO_HAVE_MPI")
endif()

set(GENERIC_MAJOR_VERSION 1)
set(GENERIC_MINOR_VERSION 50)
set(GENERIC_MICRO_VERSION 0)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/core/inc/queso.h.in 
  ${CMAKE_CURRENT_BINARY_DIR}/include/queso/queso.h)


# -----
# Properties for compilation/link
# -----

include_directories(${CMAKE_CURRENT_BINARY_DIR}/include 
  ${MPI_CXX_INCLUDE_PATH}
  ${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  )

# Don't yet need since not linking anything here...
#link_directories(${Boost_LIBRARY_DIRS} ${GSL_LIBRARY_DIRS})

# Don't force, so parent can override
set(BUILD_SHARED_LIBS TRUE CACHE BOOL "QUESO building shared libraries")


# -----
# Create libraries
# -----

set(queso_sources)
foreach(queso_subdir basic core misc stats gp surrogates)

  # Mimic autotools build process and put headers in build tree
  # Likely can eliminate...
  file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/src/${queso_subdir}/inc/*.h)
  foreach(h ${headers})
    configure_file(${h} ${CMAKE_CURRENT_BINARY_DIR}/include/queso/ COPYONLY)
  endforeach()

  file(GLOB dir_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    src/${queso_subdir}/src/*.C)
  list(APPEND queso_sources ${dir_sources})

  # Add one library per directory
  #add_library(queso_${queso_subdir} ${dir_sources})
  #target_link_libraries(queso_${queso_subdir} ${Boost_LIBRARIES})

endforeach()

# Add one library in total
add_library(queso ${queso_sources})
target_link_libraries(queso 
  ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES})
# Apply MPI compiler/link options from Dakota MPI
DakotaApplyMPISettings(queso)

install(TARGETS queso EXPORT ${ExportTarget} DESTINATION lib)
