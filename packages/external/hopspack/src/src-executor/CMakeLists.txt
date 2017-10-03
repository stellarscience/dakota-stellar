# $Id: CMakeLists.txt 217 2013-11-25 21:59:49Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-executor/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2013 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- List the source files to be built.
IF (MPI_FOUND)
    SET (src-executor_sources HOPSPACK_ExecutorMpi.cpp)

ELSE (MPI_FOUND)
    IF (MT_FOUND)
        SET (src-executor_sources
             HOPSPACK_ExecutorMultiThreaded.cpp
             HOPSPACK_ThreadedEvalWorker.cpp
            )
    ELSE (MT_FOUND)
        SET (src-executor_sources HOPSPACK_ExecutorSerial.cpp)

    ENDIF (MT_FOUND)
ENDIF (MPI_FOUND)
