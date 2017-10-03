# $Id: CMakeLists.txt 222 2014-01-21 21:30:55Z briadam $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-main/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2013 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Specify directories with required include files.
INCLUDE_DIRECTORIES (
  ${HOPSPACK_BINARY_DIR_SRC}/src-shared
  ${HOPSPACK_SRC_DIR}/src-shared
  ${HOPSPACK_SRC_DIR}/src-citizens
  ${HOPSPACK_SRC_DIR}/src-evaluator
  ${HOPSPACK_SRC_DIR}/src-executor
  ${HOPSPACK_SRC_DIR}/src-framework
)

#---- src-main source files for libhopspack
SET (src-main_sources HOPSPACK_Hopspack.cpp)


#---- Define which target executable to build.
IF (MPI_FOUND)
    #---- Build MPI-enabled executable (multiple processes).
    SET (HOPSPACK_EXEC_NAME HOPSPACK_main_mpi)

ELSE (MPI_FOUND)
    IF (MT_FOUND)
        #---- Build multi-threaded executable (single process, multiple threads).
        SET (HOPSPACK_EXEC_NAME HOPSPACK_main_threaded)

    ELSE (MT_FOUND)
        #---- Build serial executable (single process, single thread).
        SET (HOPSPACK_EXEC_NAME HOPSPACK_main_serial)

    ENDIF (MT_FOUND)
ENDIF (MPI_FOUND)

#---- List the source files to be built.
ADD_EXECUTABLE (
  ${HOPSPACK_EXEC_NAME}
    ${HOPSPACK_EXEC_NAME}.cpp
)

#---- List libraries needed to link the target.
#---- Variables like LAPACK_LIB were defined in the top-level CMakeLists.txt.
#---- The order of libraries is important because this is a static build.
TARGET_LINK_LIBRARIES (
  ${HOPSPACK_EXEC_NAME}
    hopspack
    ${LAPACK_LIBS}
    ${LAPACK_ADD_LIBS}
    ${OPSYS_LIBRARIES}
)

#---- Installation instructions.
INSTALL (TARGETS ${HOPSPACK_EXEC_NAME}
         DESTINATION ${HOPSPACK_EXEC_DEST})
