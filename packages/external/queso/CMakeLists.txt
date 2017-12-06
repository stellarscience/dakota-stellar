# Experimental QUESO CMake that builds only the sources (no tests,
# examples, etc.).  Does not properly manage all QUESO
# build/conditional compile settings.

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
  # We no longer use Boost program options, instead Getpot
  ##find_package(Boost 1.44 REQUIRED COMPONENTS program_options)
  find_package(Boost 1.44 REQUIRED)
endif()

find_package(GSL REQUIRED) # Version 1.15 or newer?
# Required by queso.h.in (cannot prefix with QUESO_)
set(GSL_LIBS ${GSL_LIBRARY_DIRS})


# -----
# CMake defines for generated headers
# -----

if(DAKOTA_HAVE_MPI)
  set(QUESO_HAVE_MPI 1)
endif()

# Boost options
# One CMake option QUESO_HAVE_BOOST to control all the following
#  QUESO_HAVE_BOOST
#  QUESO_HAVE_BOOST_MATH_SPECIAL_FUNCTIONS_HPP
#  QUESO_HAVE_BOOST_PROGRAM_OPTIONS_HPP
#  QUESO_HAVE_BOOST_SCOPED_PTR_HPP
#  QUESO_HAVE_BOOST_SHARED_PTR_HPP
set(QUESO_HAVE_BOOST 1)


# TODO: Check for C++11 support once using CMake 3.x (CXX_STANDARD)

# For now, a user can force use of C++11 features in QUESO by setting
# these variables: CMAKE_CXX_STANDARD=11, CMAKE_STANDARD_REQUIRED. And
# setting any necessary flags in CMAKE_CXX_FLAGS, e.g., -std=c++11
#
# We assume:
#  * If the standard is set and isn't C++98, that it's C++11 or newer
#  * That CMake didn't fall back to C++98 in its configuration
#    (CXX_STANDARD_REQUIRED).
#  * That all needed features are working if the standard is 11
if(CMAKE_CXX_STANDARD AND NOT CMAKE_CXX_STANDARD EQUAL 98)
  set(QUESO_HAVE_CXX11 1)
endif()

# Check for include files
include(CheckIncludeFile)
check_include_file("dlfcn.h" QUESO_HAVE_DLFCN_H)
check_include_file("inttypes.h" QUESO_HAVE_INTTYPES_H)
check_include_file("memory.h" QUESO_HAVE_MEMORY_H)
check_include_file("stdint.h" QUESO_HAVE_STDINT_H)
check_include_file("stdlib.h" QUESO_HAVE_STDLIB_H)
check_include_file("strings.h" QUESO_HAVE_STRINGS_H)
check_include_file("string.h" QUESO_HAVE_STRING_H)
check_include_file("sys/stat.h" QUESO_HAVE_SYS_STAT_H)
check_include_file("sys/types.h" QUESO_HAVE_SYS_TYPES_H)
check_include_file("unistd.h" QUESO_HAVE_UNISTD_H)

# TODO: Add other include files from config_queso.h.cmake.in
# There is a stub for QUESO_STDC_HEADERS but not sure what specific header(s)
# to key off of.  Should it be a bool AND of some of the above along
# with additional ones? The macro does not appear in the pruned Queso
# source.  RWH


# -----
# Build Information
# Definitions used by queso.h.in and config_queso.h.cmake.in
# -----

# TODO: see if we can get upstream to consolidate the two configured header files

set(GENERIC_MAJOR_VERSION 0)
set(GENERIC_MINOR_VERSION 56)
set(GENERIC_MICRO_VERSION 0)

set(BUILD_USER "anonymous")
set(BUILD_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
# set BUILD_HOST
cmake_host_system_information(RESULT BUILD_HOST QUERY HOSTNAME)
# queso.h is going away; for now set BUILD_DATE:
string(TIMESTAMP BUILD_DATE "%Y-%m-%d %H:%M:%S")

if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION )
  file( STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION Queso_GIT_REV )
  set(BUILD_VERSION ${Queso_GIT_REV})
else()
  find_package(Git)
  if(GIT_FOUND)
    get_filename_component(abs_source_dir ${Dakota_SOURCE_DIR} REALPATH)
    # Get the abbreviated SHA1 of the most recent commit
    execute_process(COMMAND ${GIT_EXECUTABLE} log --pretty=format:%h -1 
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE Queso_GIT_ABBREV_SHA1)
    # Get the date and time of the most recent commit
    execute_process(COMMAND ${GIT_EXECUTABLE} log --pretty=format:%ci -1 
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE Queso_GIT_DATETIME)
    # Extract the date
    string(SUBSTRING "${Queso_GIT_DATETIME}" 0 10 Queso_GIT_DATE)
    # Build revision string
    string(CONFIGURE "Dakota packages/external @Queso_GIT_ABBREV_SHA1@ (@Queso_GIT_DATE@)" Queso_GIT_REV)
  endif(GIT_FOUND) # GIT_FOUND

  set(BUILD_VERSION ${Queso_GIT_REV})
endif( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION )

set(VERSION 
  "${GENERIC_MAJOR_VERSION}.${GENERIC_MINOR_VERSION}.${GENERIC_MICRO_VERSION}")
set(BUILD_DEVSTATUS "Development Build")
set(CXX "${CMAKE_CXX_COMPILER}")
set(CXXFLAGS "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")


# -----
# Generate include files
# -----

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/queso)

#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/core/inc/queso.h.in
#  ${CMAKE_CURRENT_BINARY_DIR}/include/queso/queso.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config_queso.h.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/queso/config_queso.h)


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
# Create libraries and copy headers
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

# Vagaries of file globbing: version.C has a main()
list(REMOVE_ITEM queso_sources "src/core/src/version.C")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/contrib/getpot/getpot.h 
  ${CMAKE_CURRENT_BINARY_DIR}/include/queso/ COPYONLY)

# Add one library in total
add_library(queso ${queso_sources})
target_link_libraries(queso 
  ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES})
# Apply MPI compiler/link options from Dakota MPI
DakotaApplyMPISettings(queso)

install(TARGETS queso EXPORT ${ExportTarget} DESTINATION lib)
