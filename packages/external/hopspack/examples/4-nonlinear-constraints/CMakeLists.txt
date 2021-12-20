# $Id: CMakeLists.txt 217 2013-11-25 21:59:49Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/examples/4-nonlinear-constraints/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2013 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************

include(ModifyTargetRpath)

#---- Copy the parameters file to the target location.
SET (PARAMS_FILE_NAME examples/4-nonlinear-constraints/example4_params.txt)
CONFIGURE_FILE (${HOPSPACK_SOURCE_DIR}/${PARAMS_FILE_NAME}
                ${HOPSPACK_BINARY_DIR}/${PARAMS_FILE_NAME}
                COPYONLY)

#---- Define how to build the target executable.
ADD_EXECUTABLE (
  nonlinear_constraints
    nonlinear_constraints.cpp
)

#---- Installation instructions.
INSTALL (TARGETS nonlinear_constraints
         DESTINATION ${DAKOTA_EXAMPLES_INSTALL}/examples/hopspack/4-nonlinear-constraints)
ModifyTargetRpath( nonlinear_constraints 5)
INSTALL (FILES ${HOPSPACK_SOURCE_DIR}/${PARAMS_FILE_NAME}
         DESTINATION ${DAKOTA_EXAMPLES_INSTALL}/examples/hopspack/4-nonlinear-constraints)
