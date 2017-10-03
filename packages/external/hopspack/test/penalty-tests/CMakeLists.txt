# $Id: CMakeLists.txt 217 2013-11-25 21:59:49Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/test/penalty-tests/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2013 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Copy the README file to the target location.
CONFIGURE_FILE (${HOPSPACK_SOURCE_DIR}/test/penalty-tests/README.txt
                ${HOPSPACK_BINARY_DIR}/test/penalty-tests/README.txt
                COPYONLY)

#---- Specify directories with required include files.
INCLUDE_DIRECTORIES (
  ${HOPSPACK_BINARY_DIR_SRC}/src-shared
  ${HOPSPACK_SRC_DIR}/src-shared
)

#---- Define how to build the target executable.
ADD_EXECUTABLE (
  HOPSPACK_main_penalty_tests
    HOPSPACK_main_penalty_tests.cpp
)

#---- List libraries needed to link the target.
#---- LAPACK_LIB was defined in the top-level CMakeLists.txt.  It satisfies
#---- dependencies in "shared".
TARGET_LINK_LIBRARIES (
  HOPSPACK_main_penalty_tests
    hopspack
    ${LAPACK_LIBS}
    ${LAPACK_ADD_LIBS}
)

#---- Add this as an automatic test.
ADD_TEST (PenaltyTests HOPSPACK_main_penalty_tests)

#---- Installation instructions.
INSTALL (TARGETS HOPSPACK_main_penalty_tests
         DESTINATION test)
INSTALL (FILES README.txt
         DESTINATION test
         RENAME README_penalty_tests.txt
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
