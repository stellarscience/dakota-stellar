# $Id: CMakeLists.txt 218 2013-11-27 16:37:43Z briadam $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/src/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2013 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Copy the LICENSE file to the target location.
CONFIGURE_FILE (LICENSE_HOPSPACK
                ${HOPSPACK_BINARY_DIR}/LICENSE_HOPSPACK
                COPYONLY)


IF (lapack)
  #---- Probe for Fortran, without enabling, for help with linking
  #---- against some BLAS/LAPACK distributions.  Generate a header with
  #---- macro HOPSPACK_F90_GLOBAL to manage underscore convention.
  include(CheckLanguage)
  check_language(Fortran)
  set(HAVE_FORTRAN_COMPILER FALSE)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    include(FortranCInterface)
    FortranCInterface_HEADER(HOPSPACK_f90_config.h MACRO_NAMESPACE HOPSPACK_F90_)
    set(HAVE_FORTRAN_COMPILER TRUE)
  endif()
ENDIF (lapack)


# Build needs to find generated header file:
include_directories(${CMAKE_CURRENT_BINARY_DIR})


# Process all subdirectory sources into a single library (hopspack)
set(local_library hopspack)

## Top-level library sources (none).
set(lib${local_library}_sources "")

# Build a single library from all subdirectories
set(local_src_subdirs 
  src-shared src-evaluator src-executor src-citizens src-framework src-main
)

# Aggregate all the sources into a single library for installation
foreach(sd_it ${local_src_subdirs})

  add_subdirectory(${sd_it})

  # grab the trailing part of any paths
  string(REGEX MATCH "[^/]+$" last_path_name ${sd_it})

  # extract the variables ${last_path_name}_sources from the subdir's env
  get_directory_property(sd_sources DIRECTORY ${sd_it} 
    DEFINITION ${last_path_name}_sources)

  # add the sources to the library's list and includes
  foreach ( src_file ${sd_sources} )
    list ( APPEND lib${local_library}_sources "${sd_it}/${src_file}" )
  endforeach()

  # add the subdir itself to includes
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${sd_it})

  # extract the variables ${last_path_name}_includes from the subdir's env
  get_directory_property(sd_includes DIRECTORY ${sd_it} 
    DEFINITION ${last_path_name}_includes)
  foreach (include_dir ${sd_includes} )
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/${sd_it}/${include_dir}")
  endforeach()

endforeach()

# Additional include directory needed for generated file
INCLUDE_DIRECTORIES (${HOPSPACK_BINARY_DIR_SRC}/src-shared)

add_library(${local_library} ${lib${local_library}_sources})

target_link_libraries(${local_library} ${OPSYS_LIBRARIES})


#---- Installation instructions.
IF (INSTALL_LIB)
  INSTALL (TARGETS ${local_library} EXPORT ${ExportTarget} DESTINATION lib)
ENDIF (INSTALL_LIB)
IF (INSTALL_HEADERS)
  FILE (GLOB_RECURSE hopspack_headers *.hpp *.h)
  INSTALL (FILES ${CMAKE_CURRENT_BINARY_DIR}/src-shared/HOPSPACK_common.hpp
    ${hopspack_headers} DESTINATION include)
ENDIF (INSTALL_HEADERS)

INSTALL (FILES LICENSE_HOPSPACK
         DESTINATION .
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
