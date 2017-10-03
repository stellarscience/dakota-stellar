# $Id: CMakeLists.txt 217 2013-11-25 21:59:49Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/examples/2-linear-constraints/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2013 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Copy the parameters file to the target location.
SET (PARAMS_FILE_NAME examples/2-linear-constraints/example2_params.txt)
CONFIGURE_FILE (${HOPSPACK_SOURCE_DIR}/${PARAMS_FILE_NAME}
                ${HOPSPACK_BINARY_DIR}/${PARAMS_FILE_NAME}
                COPYONLY)

#---- Define how to build the target executable.
ADD_EXECUTABLE (
  linear_constraints
    linear_constraints.cpp
)

#---- Installation instructions.
INSTALL (TARGETS linear_constraints
         DESTINATION examples/2-linear-constraints)
INSTALL (FILES ${HOPSPACK_SOURCE_DIR}/${PARAMS_FILE_NAME}
         DESTINATION examples/2-linear-constraints)
