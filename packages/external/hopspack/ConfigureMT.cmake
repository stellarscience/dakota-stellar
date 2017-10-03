# $Id: ConfigureMT.cmake 149 2009-11-12 02:40:41Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/ConfigureMT.cmake $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Define a boolean user option "mt" that can be set from the
#---- command line ("-Dmt=true" or "=on" or "=yes"),
#---- or from the CMake cache.

OPTION (mt "TRUE means compile and link with multi-thread capability" FALSE)

#---- If FALSE then do not link MT.
#----
#---- If TRUE, then on return sets the variables:
#----   MT_FOUND        - for linking
#----   OPSYS_LIBRARIES - for linking
#----   HAVE_MT         - for header files


#-------------------------------------------------------------------------
IF (mt)
#-------------------------------------------------------------------------

   SET (MT_FOUND TRUE)

   INCLUDE (FindThreads)
   SET (OPSYS_LIBRARIES ${OPSYS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

   #---- The HOPSPACK header file needs to know that MPI is used.
   SET (HAVE_MT TRUE)

   MESSAGE (STATUS "HOPSPACK: Building with multi-thread capability.")


#-------------------------------------------------------------------------
ENDIF (mt)
#-------------------------------------------------------------------------
