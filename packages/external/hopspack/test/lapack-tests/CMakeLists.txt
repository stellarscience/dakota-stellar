# $Id: CMakeLists.txt 217 2013-11-25 21:59:49Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/test/lapack-tests/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2013 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Copy the README file to the target location.
CONFIGURE_FILE (${HOPSPACK_SOURCE_DIR}/test/lapack-tests/README.txt
                ${HOPSPACK_BINARY_DIR}/test/lapack-tests/README.txt
                COPYONLY)

#---- Specify directories with required include files.
INCLUDE_DIRECTORIES (
  ${HOPSPACK_BINARY_DIR_SRC}/src-shared
  ${HOPSPACK_SRC_DIR}/src-shared
)

#---- Define how to build the target executable.
ADD_EXECUTABLE (
  HOPSPACK_main_lapack_tests
    HOPSPACK_main_lapack_tests.cpp
)

#---- List libraries needed to link the target.
#---- Variables like LAPACK_LIB were defined in the top-level CMakeLists.txt.
#---- The order of libraries is important because this is a static build.
TARGET_LINK_LIBRARIES (
  HOPSPACK_main_lapack_tests
    hopspack
    ${LAPACK_LIBS}
    ${LAPACK_ADD_LIBS}
)

#---- Add this as an automatic test.
ADD_TEST (LapackTests HOPSPACK_main_lapack_tests)

#---- Installation instructions.
INSTALL (TARGETS HOPSPACK_main_lapack_tests
         DESTINATION test)
INSTALL (FILES README.txt
         DESTINATION test
         RENAME README_lapack_tests.txt
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
