# $Id: CMakeLists.txt 222 2014-01-21 21:30:55Z briadam $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2014 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Specify the CMake version required.
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.9)

#---- Declaring the project name will define ${HOPSPACK_SOURCE_DIR} to be
#---- the path to this source file, and ${HOPSPACK_BINARY_DIR} to be the path
#---- to the directory where cmake has been invoked.
PROJECT (HOPSPACK CXX C)

SET (HOPSPACK_SRC_DIR ${HOPSPACK_SOURCE_DIR}/src)

#---- Source subdirectory of the binary tree.
SET (HOPSPACK_BINARY_DIR_SRC ${HOPSPACK_BINARY_DIR}/src)

#---- Set the export name for install targets; parent packages likely want
#---- to override this to the name of their target.
SET (ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
     "Name for the export target for ${CMAKE_PROJECT_NAME}")

#---- Configure for MPI.
INCLUDE (ConfigureMPI.cmake)

#---- Configure for multi-threading.
INCLUDE (ConfigureMT.cmake)

#---- Cannot request both MPI and multi-threading as the executable.
IF (MPI_FOUND AND MT_FOUND)
    MESSAGE (FATAL_ERROR "HOPSPACK: Cannot request both -Dmpi and -Dmt.")
ENDIF (MPI_FOUND AND MT_FOUND)

#---- Configure for LAPACK.
INCLUDE (ConfigureLapack.cmake)

#---- Configure the build type (debug vs release, etc.)
INCLUDE (ConfigureBuildType.cmake)

#---- Configure system libraries.
INCLUDE (ConfigureSysLibs.cmake)

#---- Set the default install directory to be under the build directory.
IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET (CMAKE_INSTALL_PREFIX installed_HOPSPACK)
ENDIF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#---- Set the default executable install location to CMAKE_INSTALL_PREFIX/.
IF (NOT HOPSPACK_EXEC_DEST)
    SET (HOPSPACK_EXEC_DEST ".")
ENDIF()

#---- Option to enable/disable GPL licensed CDDLIB.
OPTION (HAVE_CDDLIB
        "Enable GPL-licensed package CDDLIB for degenerate constraints" TRUE)
IF (HAVE_CDDLIB)
    MESSAGE (STATUS "HOPSPACK: Including CDDLIB for degenerate constraints.")
ELSE (HAVE_CDDLIB)
    MESSAGE (STATUS "HOPSPACK: Excluding CDDLIB for degenerate constraints.")
ENDIF (HAVE_CDDLIB)


#---- List subdirectories in the order they are to be processed.
ADD_SUBDIRECTORY (src)
ADD_SUBDIRECTORY (examples)


#---- Option to allow configure-time disabling of unit tests.
OPTION (HOPSPACK_ENABLE_TESTS "Enable HOPSPACK tests?" ON)
IF (HOPSPACK_ENABLE_TESTS)
    ENABLE_TESTING()

    IF (lapack)
        ADD_SUBDIRECTORY (test/lapack-tests)
    ENDIF (lapack)
    ADD_SUBDIRECTORY (test/penalty-tests)
    ADD_SUBDIRECTORY (test/startpoint-tests)
ENDIF (HOPSPACK_ENABLE_TESTS)


#---- Options to conditionally install headers and libraries
#---- (following Trilinos naming convention).
OPTION (INSTALL_HEADERS "Install headers?" OFF)
OPTION (INSTALL_LIB "Install libraries?" OFF)
