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

# TODO: Update to use CMake BLAS/LAPACK probes and properly use their output
if(NOT DEFINED BLAS_LIBS OR NOT DEFINED LAPACK_LIBS)
  # Historically on MSVC, tried to use CMake config file approach first.
  # Could probably just use the Unix logic below instead...
  if(MSVC)
    find_package(LAPACK REQUIRED NO_MODULE)
    set(BLAS_LIBS blas)
    set(LAPACK_LIBS lapack)
  else()
    # first check for a system blas and lapack
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
endif()

FortranCInterface_HEADER(include/optpp_blas_config.h MACRO_NAMESPACE OPTPP_)


# OPT++ Teuchos management; 3 cases to consider:
#  1. Installed Trilinos
#  2. Teuchos provided elsewhere in build, e.g., Dakota
#  3. Local Teuchos (separately checked out into trilinos/)

if (NOT BUILD_IN_TRILINOS)


  # Workaround to skip finding system Trilinos until this probe is
  # simplified and follows find_package semantics. Double negative to
  # preserve historical behavior without overcomplicating things.
  if(NOT DAKOTA_NO_FIND_TRILINOS)
    # First probe for system-installed Trilinos, respecting Trilinos_DIR if set
    find_package(Trilinos QUIET)
  endif()

  if(Trilinos_FOUND)
    if(Trilinos_DIR)
      message(STATUS
	"OPT++ using previously specified Trilinos in ${Trilinos_DIR}")
    else()
      message(STATUS "OPT++ using external Trilinos")
    endif()
  else()

    # If no parent project configured Teuchos, do so, using Teuchos_DIR if set
    if(Teuchos_DIR)
      message(STATUS
	"OPT++ using previously specified Teuchos in ${Teuchos_DIR}")
    elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/trilinos)

      set(Teuchos_DIR
        ${CMAKE_CURRENT_BINARY_DIR}/trilinos/packages/teuchos)
      set(Trilinos_ENABLE_Teuchos ON CACHE BOOL
        "OPT++ enabling Trilinos Teuchos" FORCE)

      # Map key OPT++ variables to TriBITS variables
      set( TPL_BLAS_LIBRARIES ${BLAS_LIBS} )
      set( TPL_LAPACK_LIBRARIES ${LAPACK_LIBS} )
      # OPT++ doesn't use any Teuchos MPI features; may want to force off
      #set( TPL_ENABLE_MPI ${DAKOTA_HAVE_MPI} )
# This doesn't do as name implies; setting OFF doesn't generate Config.cmake
# at all; doesn't just control whether installed!  Want Config.cmake in build
#        set(Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES OFF CACHE BOOL
#	  "OPT++ is the top-level package; don't write Trilinos config files")

      message(STATUS "OPT++ setting Teuchos_DIR to ${Teuchos_DIR}")
      add_subdirectory(trilinos)

    else()
      message(WARNING	"OPT++ requires teuchos, but could not find local nor parent copy!")
    endif() # Teuchos_DIR

    # Additional setting to prevent multiple targets with the same name
    set(Trilinos_TARGETS_IMPORTED 1)

    find_package( Teuchos NO_MODULE REQUIRED )

  endif() # Trilinos_DIR

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
