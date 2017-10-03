# $Id: ConfigureMPI.cmake 149 2009-11-12 02:40:41Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/ConfigureMPI.cmake $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Define a boolean user option "mpi" that can be set from the
#---- command line ("-Dmpi=true" or "=on" or "=yes"),
#---- or from the CMake cache.

OPTION (mpi "TRUE means compile and link with MPI" FALSE)

#---- If FALSE then do not link MPI.
#----
#---- If TRUE then find an MPI compiler and library for linking.
#---- CMake includes a FindMPI module that searches in default locations.
#---- If it fails, or to override the choice with a preferred version,
#---- set using command line options.
#----
#---- Command line options:
#----   -DMPI_COMPILER= full path to special MPI compiler
#----   -DMPI_INCLUDE_PATH= location of mpi.h
#----   -DMPI_LIBRARY= name of special MPI library
#---- It may not be necessary to supply all three options.
#----
#---- If TRUE, then on return sets the variables:
#----   MPI_FOUND - for linking
#----   HAVE_MPI  - for header files


#-------------------------------------------------------------------------
IF (NOT mpi)
#-------------------------------------------------------------------------
    MESSAGE (STATUS "HOPSPACK: Using the standard C++ compiler (no MPI):")
    MESSAGE (STATUS "  " ${CMAKE_CXX_COMPILER})

#-------------------------------------------------------------------------
ELSE (NOT mpi)
#-------------------------------------------------------------------------

    #---- Enable the following options to force a particular version of MPI.
    #---- Enable by uncommenting, or passing on the command line.
    #---- It may not be necessary to supply all three.
#    SET (MPI_COMPILER /tmp/mpich1.2.7p1/bin/mpicxx)
#    SET (MPI_INCLUDE_PATH location of mpi.h)
#    SET (MPI_LIBRARY any mpi library to link)

    IF (MPI_COMPILER OR MPI_INCLUDE_PATH OR MPI_LIBRARY)
        #---- The user defined MPI.
        SET (MPI_FOUND TRUE)

        #---- The HOPSPACK header file needs to know that MPI is used.
        SET (HAVE_MPI TRUE)

        INCLUDE_DIRECTORIES (${MPI_INCLUDE_PATH})

    ELSE (MPI_COMPILER OR MPI_INCLUDE_PATH OR MPI_LIBRARY)
        #---- Search standard locations for an MPI installation.
        #---- FIND_PACKAGE invokes cmake/share/cmake-2.6/Modules/FindMPI.cmake
        #---- to look for MPI, setting MPI_FOUND as TRUE or FALSE.
        #---- See comments in FindMPI.cmake for more information.
        #---- If successful, FIND_PACKAGE will also define the following:
        #----   MPI_INCLUDE_PATH
        #----   MPI_LIBRARY
        #----   MPI_COMPILE_FLAGS
        #----   MPI_LINK_FLAGS
        #----   MPI_EXECUTABLE
        #----   MPIEXEC
        #----   MPIEXEC_MAX_NUMPROCS, MPIEXEC_NUMPROC_FLAG
        #----   MPIEXEC_PREFLAGS, MPIEXEC_POSTFLAGS
        FIND_PACKAGE (MPI)

        IF (MPI_FOUND)
            #---- The HOPSPACK header file needs to know that MPI is used.
            SET (HAVE_MPI TRUE)
        ELSE (MPI_FOUND)
            MESSAGE (FATAL_ERROR "HOPSPACK: MPI requested but not found.")
        ENDIF (MPI_FOUND)

    ENDIF (MPI_COMPILER OR MPI_INCLUDE_PATH OR MPI_LIBRARY)

    #---- If an MPI compiler was specified, then use it.
    IF (MPI_COMPILER)
        INCLUDE (CMakeForceCompiler)
        CMAKE_FORCE_C_COMPILER (${MPI_COMPILER} UndefinedCompilerId)
        CMAKE_FORCE_CXX_COMPILER (${MPI_COMPILER} UndefinedCompilerId)
    ENDIF (MPI_COMPILER)

    MESSAGE (STATUS "HOPSPACK: Building with MPI, compiler is:")
    MESSAGE (STATUS "  " ${CMAKE_CXX_COMPILER})


#-------------------------------------------------------------------------
ENDIF (NOT mpi)
#-------------------------------------------------------------------------
