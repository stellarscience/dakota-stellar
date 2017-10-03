# $Id: CMakeLists.txt 194 2011-07-26 21:38:55Z briadam $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2010 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Create HOPSPACK_common.hpp by modifying the "HAVE_varname" define
#---- statements for this particular build.
CONFIGURE_FILE (
    ${HOPSPACK_SRC_DIR}/src-shared/HOPSPACK_common.hpp.cmake
    ${HOPSPACK_BINARY_DIR_SRC}/src-shared/HOPSPACK_common.hpp
)

#---- List the source files to be built.
SET (BASIC_LIST
    HOPSPACK_CallbackToMediator.cpp
    HOPSPACK_DataPoint.cpp
    HOPSPACK_float.cpp
    HOPSPACK_LapackWrappers.cpp
    HOPSPACK_LinConstr.cpp
    HOPSPACK_Matrix.cpp
    HOPSPACK_NonlConstrPenalty.cpp
    HOPSPACK_ParameterEntry.cpp
    HOPSPACK_ParameterList.cpp
    HOPSPACK_Print.cpp
    HOPSPACK_ProblemDef.cpp
    HOPSPACK_SolveLinConstrProj.cpp
    HOPSPACK_SystemTimer.cpp
    HOPSPACK_utils.cpp
    HOPSPACK_Vector.cpp
)
IF (MPI_FOUND)
    SET (SOURCE_LIST ${BASIC_LIST} HOPSPACK_GenProcComm.cpp)
ELSE (MPI_FOUND)
    IF (MT_FOUND)
        SET (SOURCE_LIST ${BASIC_LIST} HOPSPACK_ThreadRunnable.cpp
                                       HOPSPACK_ThreadSimpleLock.cpp
                                       HOPSPACK_ThreadSynchObject.cpp)
    ELSE (MT_FOUND)
        SET (SOURCE_LIST ${BASIC_LIST})
    ENDIF (MT_FOUND)
ENDIF (MPI_FOUND)

#---- src-shared directory sources.
set (src-shared_sources ${SOURCE_LIST})
