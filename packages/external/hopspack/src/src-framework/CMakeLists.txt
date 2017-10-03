# $Id: CMakeLists.txt 194 2011-07-26 21:38:55Z briadam $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- List the source files to be built.
SET (src-framework_sources
  HOPSPACK_CacheManager.cpp
  HOPSPACK_CachePoint.cpp
  HOPSPACK_Conveyor.cpp
  HOPSPACK_ConveyorList.cpp
  HOPSPACK_Mediator.cpp
  HOPSPACK_EvalCounter.cpp
  HOPSPACK_ScaledComparison.cpp
)
