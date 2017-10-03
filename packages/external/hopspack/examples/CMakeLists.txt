# $Id: CMakeLists.txt 217 2013-11-25 21:59:49Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/examples/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2010 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Copy the README file to the target location.
CONFIGURE_FILE (${HOPSPACK_SOURCE_DIR}/examples/README.txt
                ${HOPSPACK_BINARY_DIR}/examples/README.txt
                COPYONLY)

#---- List subdirectories in the order they are to be processed.
ADD_SUBDIRECTORY (1-var-bnds-only)
ADD_SUBDIRECTORY (2-linear-constraints)
ADD_SUBDIRECTORY (3-degen-linear-constraints)
ADD_SUBDIRECTORY (4-nonlinear-constraints)
ADD_SUBDIRECTORY (5-multi-start)

#---- Installation instructions.
INSTALL (FILES README.txt
         DESTINATION examples
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
