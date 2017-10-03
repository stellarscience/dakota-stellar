# $Id: CMakeLists.txt 220 2014-01-02 21:24:59Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/examples/linked-evaluator-example/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2014 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- List the source files to be built.
SET (src-evaluator_sources
  HOPSPACK_EvaluatorFactory.cpp
  ExampleLinkedEvaluator.cpp
)
