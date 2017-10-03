# $Id: CMakeLists.txt 217 2013-11-25 21:59:49Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/test/startpoint-tests/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#               Copyright 2009-2013 Sandia Corporation
#
#   Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
#   license for use of this work by or on behalf of the U.S. Government.
# ************************************************************************


#---- Copy the README file to the target location.
SET (PARAMS_FILE_NAME test/startpoint-tests/README.txt)
CONFIGURE_FILE (${HOPSPACK_SOURCE_DIR}/${PARAMS_FILE_NAME}
                ${HOPSPACK_BINARY_DIR}/${PARAMS_FILE_NAME}
                COPYONLY)

#---- Specify directories with required include files.
INCLUDE_DIRECTORIES (
  ${HOPSPACK_BINARY_DIR_SRC}/src-shared
  ${HOPSPACK_SRC_DIR}/src-shared
)

#---- Define how to build the target executable.
ADD_EXECUTABLE (HOPSPACK_main_startpoint_tests
                  HOPSPACK_main_startpoint_tests.cpp)

#---- List libraries needed to link the target.
#---- LAPACK_LIB was defined in the top-level CMakeLists.txt.  It satisfies
#---- dependencies in "shared".
TARGET_LINK_LIBRARIES (
  HOPSPACK_main_startpoint_tests
    hopspack
    ${LAPACK_LIBS}
    ${LAPACK_ADD_LIBS}
)

#---- Add this as an automatic test.
ADD_TEST (StartpointTests HOPSPACK_main_startpoint_tests)

#---- Installation instructions.
INSTALL (TARGETS HOPSPACK_main_startpoint_tests
         DESTINATION test)
INSTALL (FILES README.txt
         DESTINATION test
         RENAME README_startpoint_tests.txt
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
