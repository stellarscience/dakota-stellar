cmake_minimum_required(VERSION 2.8.4)
project(optpp C CXX Fortran)

set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_MODULE_PATH}
)

include(FortranCInterface)

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

option(DAKOTA_OPTPP "Enable DAKOTA enhancements to OPT++" TRUE)
if (DAKOTA_OPTPP)
  add_definitions(-DDAKOTA_OPTPP)
else()
  find_library(npsol_lib libnpsol.a /usr/local/lib/libnpsol.a)
endif (DAKOTA_OPTPP)

# first check for a system blas and lapack
if(NOT DEFINED BLAS_LIBS OR NOT DEFINED LAPACK_LIBS)
  if(NOT DEFINED BLAS_LIBS)
    find_library(BLAS_LIBS blas)
  endif()
  if(NOT DEFINED LAPACK_LIBS)
    find_library(LAPACK_LIBS lapack)
  endif()
  if(NOT BLAS_LIBS OR NOT LAPACK_LIBS)
    # if not a system blas and lapack, then look for a cmake built LAPACK
    # with find_package
    find_package(LAPACK REQUIRED NO_MODULE)
    set(BLAS_LIBS blas)
    set(LAPACK_LIBS lapack)
  endif()
endif()

FortranCInterface_HEADER(include/cblas_config.h MACRO_NAMESPACE OPTPP_)


# OPT++ Teuchos management; 3 cases to consider:
# * Installed Trilinos
# * Teuchos provided elsewhere, e.g., DAKOTA
# * local Teuchos (separately checked out into teuchos/)

# when building inside Trilinos, the path to Teuchos will already be set
if (NOT BUILD_IN_TRILINOS)
  # first probe for system-installed Trilinos
  # this will respect Trilinos_DIR if already set
  find_package(Trilinos QUIET)

  if ( NOT Trilinos_DIR )

    # if no one has configured a src Teuchos via Teuchos_DIR,
    # see if a local one is available
    if(Teuchos_DIR)
      message( "in ${CMAKE_CURRENT_BINARY_DIR} Teuchos_DIR already set to ${Teuchos_DIR}" )
    else()
      # teuchos checked out under OPT++
      if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/teuchos/cmake/tribits )
	set( Teuchos_DIR ${CMAKE_CURRENT_BINARY_DIR}/teuchos/packages/teuchos )
	set( Trilinos_ENABLE_Teuchos ON CACHE BOOL "Enable teuchos only" FORCE )
        set( TPL_BLAS_LIBRARIES ${BLAS_LIBS} )
        set( TPL_LAPACK_LIBRARIES ${LAPACK_LIBS} )
	add_subdirectory(teuchos)
      elseif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/teuchos)
        # "Historical" directory hierarchy, before TriBITS
	set( Teuchos_DIR ${CMAKE_CURRENT_BINARY_DIR}/teuchos )
	message( "Setting Teuchos_DIR to ${Teuchos_DIR}" )
	add_subdirectory(teuchos)
	set(Teuchos_LIBRARY_DIRS "")
      else()
	message(WARNING 
	  "OPT++ requires teuchos, but could not find local nor parent copy!")
      endif()
    endif() # Teuchos_DIR

    if( EXISTS ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake )
      # Move the "TriBITS" config file to teuchos subdir (rather than the BuildRoot)
      file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/teuchos/TrilinosConfig.cmake)

      # Additional setting to prevent multiple targets with the same name
      set(Trilinos_TARGETS_IMPORTED 1)
    endif()

    find_package( Teuchos NO_MODULE REQUIRED )
    #message( STATUS "Teuchos_LIBRARIES for opt++ link: ${Teuchos_LIBRARIES}" )
  else()
    message( "Using system trilinos in ${Trilinos_DIR}" )
  endif() # NOT Trilinos_DIR
endif() # NOT BUILD_IN_TRILINOS
link_directories(${Teuchos_LIBRARY_DIRS})


#the following configured file doesn't do anything except generate an empty
#source file that interfaces still wants to include...
configure_file(
  ${optpp_SOURCE_DIR}/cmake/OPT++_config_cmake.h.in
  ${optpp_BINARY_DIR}/include/OPT++_config.h
  COPYONLY
)


# WJB - ToDo: work with Zack to ensure consistency across Dakota/Opt++/Acro
#if(NOT BUILT_BY_ACRO)
#endif()

option(OPTPP_HAVE_MPI "Should OPTPP enable MPI?" OFF)
if(OPTPP_HAVE_MPI)
  find_package(MPI)

  if(MPI_FOUND)
    add_definitions(-DOPTPP_HAVE_MPI)
    include_directories(${MPI_INCLUDE_PATH})
  else()
    message(FATAL_ERROR " Opt++ has MPI bindings enabled but cannot find an
            installed MPI. Please set MPIHOME environment var to MPI location.")
  endif(MPI_FOUND)
endif(OPTPP_HAVE_MPI)

#option(OPTPP_SHARED "Turn ON if you have a shared file system" OFF)
#if(OPTPP_SHARED)
#  set(SHARED TRUE)
#  if(NOT MPI_FOUND)
#    message(FATAL_ERROR " optpp cannot be built on a shared file system
#    without MPI.  Please set OPTPP_HAVE_MPI to ON")
#  endif()
#endif()

include(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(times HAVE_TIMES)
if(HAVE_TIMES)
  add_definitions(-DHAVE_TIMES)
endif(HAVE_TIMES)

include(CheckIncludeFiles)
CHECK_INCLUDE_FILES(sys/times.h HAVE_SYS_TIMES_H)
if(HAVE_SYS_TIMES_H)
  add_definitions(-DHAVE_SYS_TIMES_H)
endif(HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILES(sys/time.h HAVE_SYS_TIME_H)
if(HAVE_SYS_TIME_H)
  add_definitions(-DHAVE_SYS_TIME_H)
endif(HAVE_SYS_TIME_H)
add_definitions(-DHAVE_STD)
    
if(DAKOTA_OPTPP)
  set(OPTPP_INCLUDE_DIRS
    ${optpp_BINARY_DIR}/include
    ${optpp_SOURCE_DIR}/include
    ${Teuchos_INCLUDE_DIRS}
    PARENT_SCOPE
  )
endif()

include_directories(
  ${optpp_BINARY_DIR}/include
  ${optpp_SOURCE_DIR}/include
  ${Teuchos_INCLUDE_DIRS}
  ${AMPL_INCLUDE_DIRS}
)

remove_definitions(-DHAVE_CONFIG_H)

add_subdirectory(src)

option(OPTPP_ENABLE_TESTS "Enable OPT++ tests." ON)
if(OPTPP_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif(OPTPP_ENABLE_TESTS)

file(GLOB optpp_include_files include/*.h)
install(FILES ${optpp_BINARY_DIR}/include/OPT++_config.h
  ${optpp_include_files} DESTINATION include)
