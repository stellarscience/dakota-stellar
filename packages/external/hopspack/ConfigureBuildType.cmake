# $Id: ConfigureBuildType.cmake 221 2014-01-06 21:48:13Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/ConfigureBuildType.cmake $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2014 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Define a boolean user option "HOPSPACK_DEBUG" that can be set from the
#---- command line ("-DHOPSPACK_DEBUG<:BOOL>=true" or "=on" or "=yes"),
#---- or from the CMake cache.  Default to false, so it doesn't interfere with
#---- CMAKE_BUILD_TYPE if set in another context.

OPTION (HOPSPACK_DEBUG "TRUE means compile as debug objects" FALSE)

IF (HOPSPACK_DEBUG)
    SET (CMAKE_BUILD_TYPE DEBUG)
    MESSAGE (STATUS "HOPSPACK: Forcing CMAKE_BUILD_TYPE=DEBUG.")
    MESSAGE (STATUS "HOPSPACK: Makefiles will create debug objects.")
ENDIF (HOPSPACK_DEBUG)


#-------------------------------------------------------------------------


#---- Uncomment the first line to make the native compiler calls visible.
#---- For Windows MSVC also uncomment the the next 2 lines; however, even
#---- then not all make operations are visible.
#SET (CMAKE_VERBOSE_MAKEFILE ON)
#SET (CMAKE_START_TEMP_FILE "")
#SET (CMAKE_END_TEMP_FILE "")

IF (CMAKE_VERBOSE_MAKEFILE)
    MESSAGE (STATUS "HOPSPACK: Makefiles requested to display compile commands.")
ENDIF (CMAKE_VERBOSE_MAKEFILE)
