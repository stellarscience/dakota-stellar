#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

cmake_minimum_required(VERSION 3.12)

set( CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_MODULE_PATH}
)

# Pecos requires C++11
include(DakotaCxxOptions)
dakota_cxx_standard()

project("Pecos" C CXX Fortran)

# --- CMake modules ---
include(FortranCInterface)
include(CTest)


# Libraries from Pecos itself
set(Pecos_LIBS)
# Libraries from CMake-built packages required by Pecos
set(Pecos_PKG_LIBS)
# External third-party libraries installed on the system
set(Pecos_TPL_LIBS)

# --- Global options ---

option(BUILD_SHARED_LIBS "Build Pecos with shared libraries?" ON)

remove_definitions("-DHAVE_CONFIG_H")

if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS
      "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")


# --- External TPLs ---

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

find_package(Boost 1.58 REQUIRED)

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
	"Pecos using previously specified Trilinos in ${Trilinos_DIR}")
    else()
      message(STATUS "Pecos using external Trilinos")
    endif()
  else()

    # If no parent project configured Teuchos, do so, using Teuchos_DIR if set
    if(Teuchos_DIR)
      message(STATUS
	"Pecos using previously specified Teuchos in ${Teuchos_DIR}")
    elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/packages/trilinos)

      set(Teuchos_DIR
        ${CMAKE_CURRENT_BINARY_DIR}/packages/trilinos/packages/teuchos)
      set(Trilinos_ENABLE_Teuchos ON CACHE BOOL
        "Pecos enabling Trilinos Teuchos" FORCE)

      # Map key Pecos variables to TriBITS variables
      set( TPL_BLAS_LIBRARIES ${BLAS_LIBS} )
      set( TPL_LAPACK_LIBRARIES ${LAPACK_LIBS} )
      # Pecos doesn't use any Teuchos MPI features; may want to force off
      #set( TPL_ENABLE_MPI ${DAKOTA_HAVE_MPI} )
      # This doesn't do as name implies; setting OFF doesn't generate Config.cmake
      # at all; doesn't just control whether installed!  Want Config.cmake in build
      #        set(Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES OFF CACHE BOOL
      #	  "Pecos is the top-level package; don't write Trilinos config files")

      message(STATUS "Pecos setting Teuchos_DIR to ${Teuchos_DIR}")

      # This mirrors the Trilinos release process setting, so Pecos
      # can work with Trilinos master, when needed.  Specifically it
      # will suppress errors about missing packages.
      set(Trilinos_ENABLE_DEVELOPMENT_MODE OFF CACHE BOOL
	"Pecos disabling Trilinos development mode")
      # This would be lighter weight, but might result in strong warnings
      ##set(Trilinos_ASSERT_MISSING_PACKAGES OFF CACHE BOOL "Dakota being lazy")

      set(Trilinos_GENERATE_REPO_VERSION_FILE OFF CACHE BOOL
        "Pecos disabling generation of TrilinosRepoVersion.txt")

      # BMA: Since Pecos surrogates not using complex yet, disabled this
      ##set(Trilinos_ENABLE_COMPLEX_DOUBLE TRUE CACHE BOOL
      ##  "Pecos requires some complex blas wrappers")

      add_subdirectory(packages/trilinos)

    else()
      message(WARNING
	"Pecos requires teuchos, but could not find local nor parent copy!")
    endif() # Teuchos_DIR

    # Additional setting to prevent multiple targets with the same name
    set(Trilinos_TARGETS_IMPORTED 1)

    find_package( Teuchos NO_MODULE REQUIRED )

  endif() # Trilinos_DIR

endif() # NOT BUILD_IN_TRILINOS

link_directories(${Teuchos_LIBRARY_DIRS})

# Unconditional dependencies
# (Teuchos is in TPL even if we're building it to avoid export issues
# with Trilinos and Dakota)
list(APPEND Pecos_TPL_LIBS
  ${Teuchos_LIBRARIES} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS})


# --- Potential locally included packages ---

option(HAVE_FFT "Use either the fftw or dfftpack packages" ON)
option(HAVE_DFFTPACK "Build the dfftpack package.  This OPTION has no effect if
  HAVE_FFT is OFF" ON)
option(HAVE_FFTW "Find and use an installed fftw package. This OPTION has no
  effect if HAVE_FFT is OFF" OFF)

option(HAVE_LHS "Build the LHS package." ON)
option(HAVE_SPARSE_GRID "Build the VPISparseGrid package." ON)

add_subdirectory(packages)

# --- Dependency configuration based on previously defined packages ---

function(CheckPackage package)
  if(HAVE_${package})
    add_definitions("-DHAVE_${package}")
    string(TOLOWER ${package} package_lower)
    set(Pecos_PKG_LIBS ${Pecos_PKG_LIBS} ${package_lower} PARENT_SCOPE)
  endif()
endfunction(CheckPackage)

CheckPackage(LHS)
if(HAVE_LHS)
  list(APPEND Pecos_PKG_LIBS lhs_mods lhs_mod)
endif(HAVE_LHS)

if(HAVE_FFT)
  CheckPackage(DFFTPACK)

  if(HAVE_FFTW)
    # BMA: I couldn't get this to work with the ExternalProject; commenting:
    #find_package(Fftw3)
    list(APPEND Pecos_PKG_LIBS fftw3)
    link_directories("${FFTW_BINARY_DIR}/.libs")
  endif(HAVE_FFTW)

endif(HAVE_FFT)

if(HAVE_SPARSE_GRID)
  add_definitions("-DHAVE_SPARSE_GRID")
  list(APPEND Pecos_PKG_LIBS sparsegrid)
endif(HAVE_SPARSE_GRID)


# --- Options for Pecos components

##option(PECOS_ENABLE_SAMPLING "Enable sampling module" ON)
option(PECOS_ENABLE_ORTHOGPOLY "Enable core pecos module" ON)
option(PECOS_ENABLE_SURROGATES "Enable surrogates module" ON)
option(PECOS_ENABLE_PYTHONAPI "Enable Python wrappers" OFF)
option(PECOS_ENABLE_TESTS "Enable Pecos tests?" ON)

# Foundational utilities
add_subdirectory(util)
list(APPEND Pecos_LIBS pecos_util)

# Core Pecos
if(PECOS_ENABLE_ORTHOGPOLY)
  add_subdirectory(src)
  list(APPEND Pecos_LIBS pecos_src)
endif()

# Tests for core Pecos
if(BUILD_TESTING AND PECOS_ENABLE_TESTS AND PECOS_ENABLE_ORTHOGPOLY)
  add_subdirectory(test)
endif()

# Modular surrogate model library
# TODO: optionally depends on Pecos core (need tests with and without)
if(PECOS_ENABLE_SURROGATES)
  add_subdirectory(surrogates)

##  list(APPEND Pecos_LIBS models)
##  if(PECOS_ENABLE_ORTHOGPOLY)
##    list(APPEND Pecos_LIBS pecos_wrapper)
##  endif()
endif()

# Python wrappers are optional; this approach is having CMake invoke
# distutils setup.py. Must add this subdirectory after all other core
# libraries. TODO: check for required dependencies.
if(PECOS_ENABLE_PYTHONAPI)
  add_subdirectory(python)
endif()

set(Pecos_LINK_LIBRARIES "${Pecos_LIBS};${Pecos_PKG_LIBS}" CACHE INTERNAL
  "Pecos link libraries")
