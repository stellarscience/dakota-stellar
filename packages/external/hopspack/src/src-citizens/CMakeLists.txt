# $Id: CMakeLists.txt 217 2013-11-25 21:59:49Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/CMakeLists.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2013 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Specify subdirectories with required include files.
SET (src-citizens_includes
  citizen-gss-common
  citizen-gss
  citizen-gss/cddlib
  citizen-gss-ms
  citizen-gss-nlc
)

#---- List the source files to be built.
SET (src-citizens_sources
  HOPSPACK_Citizen.cpp
  citizen-gss-common/HOPSPACK_gssChildReturnCodes.cpp
  citizen-gss/HOPSPACK_CitizenGSS.cpp
  citizen-gss/HOPSPACK_GssDirections.cpp
  citizen-gss/HOPSPACK_GssIterator.cpp
  citizen-gss/HOPSPACK_GssList.cpp
  citizen-gss/HOPSPACK_GssPoint.cpp
  citizen-gss-ms/HOPSPACK_CitizenGssMS.cpp
  citizen-gss-ms/HOPSPACK_GeneratorTBD.cpp
  citizen-gss-ms/HOPSPACK_MultiStartRepository.cpp
  citizen-gss-ms/HOPSPACK_PointGeneratorInterface.cpp
  citizen-gss-nlc/HOPSPACK_CitizenGssNlc.cpp
)

#---- Conditionally include the GPL CDDLIB package
IF (HAVE_CDDLIB)
  LIST (APPEND src-citizens_sources 
    citizen-gss/cddlib/cddcore.c
    citizen-gss/cddlib/cddio.c
    citizen-gss/cddlib/cddlib.c
    citizen-gss/cddlib/cddlp.c
    citizen-gss/cddlib/cddmp.c
    citizen-gss/cddlib/cddproj.c
    citizen-gss/cddlib/setoper.c
    citizen-gss/HOPSPACK_CddLibWrapper.c
  )
  ADD_DEFINITIONS ("-DHAVE_CDDLIB")

  CONFIGURE_FILE (citizen-gss/cddlib/LICENSE_CDDLIB
                  ${HOPSPACK_BINARY_DIR}/LICENSE_CDDLIB
                  COPYONLY)
  INSTALL (FILES citizen-gss/cddlib/LICENSE_CDDLIB
           DESTINATION .
           PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
ENDIF()
