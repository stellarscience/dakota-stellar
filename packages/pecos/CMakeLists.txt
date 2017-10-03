cmake_minimum_required(VERSION 2.8)

project("Pecos" C CXX Fortran)

# --- CMake modules ---
include(FortranCInterface)
include(CTest)

set( CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/package_arch
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config_tests
  ${CMAKE_MODULE_PATH}
)


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

# Mandate a system or user-provided Boost
if(Boost_DIR)
  include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
else()
  find_package(Boost REQUIRED)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
  endif()
endif(Boost_DIR)

# when building inside Trilinos, the path to Teuchos will already be set
if (NOT BUILD_IN_TRILINOS)

  # first probe for system-installed Trilinos
  # this will respect Trilinos_DIR if already set
  find_package(Trilinos QUIET)
  
  if ( NOT Trilinos_DIR )

    # if no one has configured a local src Teuchos, do so
    # this will respect Teuchos_DIR if already set
    if( NOT Teuchos_DIR )
      if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/packages/teuchos/cmake/tribits )
	set( Teuchos_DIR
          ${CMAKE_CURRENT_BINARY_DIR}/packages/teuchos/packages/teuchos )
	set( Trilinos_ENABLE_Teuchos ON CACHE BOOL
          "ENABLE TriBITS build of teuchos" FORCE )

	# Map key Dakota variables to TriBITS variables
	set( TPL_ENABLE_MPI ${DAKOTA_HAVE_MPI} )
	set( TPL_BLAS_LIBRARIES ${BLAS_LIBS} )
	set( TPL_LAPACK_LIBRARIES ${LAPACK_LIBS} )

      else()
	# Use "historical" directory hierarchy, before TriBITS
	set( Teuchos_DIR ${CMAKE_CURRENT_BINARY_DIR}/packages/teuchos )
      endif() # tribits source directory EXISTS

      message( "Setting Teuchos_DIR to ${Teuchos_DIR}" )

      # Additional setting to prevent multiple targets with the same name
      set(Trilinos_TARGETS_IMPORTED 1)

      add_subdirectory(packages/teuchos)
      #set(Teuchos_LIBRARY_DIRS "")

    else()
      message( "in ${CMAKE_CURRENT_BINARY_DIR} Teuchos_DIR already set to ${Teuchos_DIR}" )
    endif() # NOT Teuchos_DIR

    if( EXISTS ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake )
      # Move the "TriBITS" config file to teuchos subdir (not the BuildRoot)
      file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake
	${CMAKE_CURRENT_BINARY_DIR}/packages/teuchos/TrilinosConfig.cmake)
    endif()

    #set(Teuchos_LIBRARIES teuchos)
    find_package( Teuchos NO_MODULE REQUIRED )
    message( STATUS "Teuchos_LIBRARIES for pecos link: ${Teuchos_LIBRARIES}" )
    
  else()
    message( "Using system trilinos in ${Trilinos_DIR}" )
  endif() # NOT Trilinos_DIR
endif() # NOT BUILD_IN_TRILINOS

link_directories(${Teuchos_LIBRARY_DIRS})


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

# Keep track of what libraries should be included in libpecos
set(SUBLIBS pecos_src)

function(CheckPackage package)
  if(HAVE_${package})
    add_definitions("-DHAVE_${package}")
    string(TOLOWER ${package} package_lower)
    set(SUBLIBS ${SUBLIBS} ${package_lower} PARENT_SCOPE)
  endif()
endfunction(CheckPackage)

CheckPackage(LHS)
if(HAVE_LHS)
  list(APPEND SUBLIBS lhs_mods lhs_mod)
endif(HAVE_LHS)

if(HAVE_FFT)
  CheckPackage(DFFTPACK)

  if(HAVE_FFTW)
    # BMA: I couldn't get this to work with the ExternalProject; commenting:
    #find_package(Fftw3)
    list(APPEND SUBLIBS fftw3)
    link_directories("${FFTW_BINARY_DIR}/.libs")
  endif(HAVE_FFTW)

endif(HAVE_FFT)

if(HAVE_SPARSE_GRID)
  add_definitions("-DHAVE_SPARSE_GRID")
  set(SUBLIBS ${SUBLIBS} sparsegrid)
endif(HAVE_SPARSE_GRID)

add_subdirectory(src)

option(PECOS_ENABLE_TESTS "Enable Pecos-specific tests?" ON)
if(BUILD_TESTING AND PECOS_ENABLE_TESTS)
  add_subdirectory(test)
endif()

# libpecos.la must be built PRIOR to the linkage of test executables
# WJB- ToDo: prefer no empty C++ files, so consult with Kitware for a better way
configure_file("${Pecos_SOURCE_DIR}/cmake/empty.cxx.in"
  "${Pecos_BINARY_DIR}/empty.cxx"
  @ONLY
)

add_library(pecos empty.cxx)
target_link_libraries(pecos ${SUBLIBS})
install(TARGETS pecos EXPORT ${ExportTarget} DESTINATION lib)

set(Pecos_LINK_LIBRARIES ${SUBLIBS} CACHE INTERNAL "Pecos link libraries")

