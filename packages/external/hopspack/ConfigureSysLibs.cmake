# $Id: ConfigureSysLibs.cmake 149 2009-11-12 02:40:41Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/ConfigureSysLibs.cmake $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Define system libraries required by HOPSPACK, and pass the information
#---- on to C++ header files as appropriate.
#----
#---- On return sets the variables:
#----   OPSYS_LIBRARIES


#---- Real time timers are enabled if system libraries can be found.
#---- The HOPSPACK header file is informed thru HAVE_REALTIME_CLOCK.
SET (HAVE_REALTIME_CLOCK TRUE)
IF (WIN32)
    FIND_LIBRARY (WINMM_LIBRARY winmm ${WIN32_LIBRARY_SEARCHPATHS})
    IF (NOT WINMM_LIBRARY)
        MESSAGE (STATUS "HOPSPACK: Could not find WinMM.lib, timers disabled.")
        SET (HAVE_REALTIME_CLOCK FALSE)
    ELSE (NOT WINMM_LIBRARY)
        SET (OPSYS_LIBRARIES ${OPSYS_LIBRARIES} ${WINMM_LIBRARY})
    ENDIF (NOT WINMM_LIBRARY)
ENDIF (WIN32)
