# $Id: CMakeLists.txt 217 2013-11-25 21:59:49Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/examples/5-multi-start/CMakeLists.txt $
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
SET (PARAMS_FILE_NAME examples/5-multi-start/example5_params.txt)
CONFIGURE_FILE (${HOPSPACK_SOURCE_DIR}/${PARAMS_FILE_NAME}
                ${HOPSPACK_BINARY_DIR}/${PARAMS_FILE_NAME}
                COPYONLY)

#---- Define how to build the target executable.
ADD_EXECUTABLE (
  multi_start
    multi_start.cpp
)

#---- Installation instructions.
INSTALL (TARGETS multi_start
         DESTINATION ${DAKOTA_EXAMPLES_INSTALL}/examples/hopspack/5-multi-start)
ModifyTargetRpath( multi_start 5)
INSTALL (FILES ${HOPSPACK_SOURCE_DIR}/${PARAMS_FILE_NAME}
         DESTINATION ${DAKOTA_EXAMPLES_INSTALL}/examples/hopspack/5-multi-start)
