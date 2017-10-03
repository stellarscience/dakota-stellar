# $Id: ConfigureLapack.cmake 217 2013-11-25 21:59:49Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/ConfigureLapack.cmake $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009-2013 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


#---- Define a boolean user option "lapack" that can be set from the
#---- command line ("-Dlapack=true" or "=on" or "=yes"),
#---- or from the CMake cache.

OPTION (lapack "TRUE means link with LAPACK/BLAS" TRUE)

#---- If FALSE then do not link LAPACK, with the consequence that HOPSPACK
#---- cannot accept problems with linear constraints.
#----
#---- If TRUE then find an LAPACK/BLAS library for linking.
#---- To link a preferred library, set using command line options or edit
#---- this file as described below in the comments after "ELSE (NOT lapack)".
#----
#---- Command line options:
#----   -DLAPACK_LIBS=
#----   -DLAPACK_ADD_LIBS=
#---- The second option is for any additional libraries required to resolve
#---- missing symbols in the LAPACK libraries.
#---- Remember that arguments are usually enclosed in double quotes, and
#---- multiple libraries are separated by a semicolon; eg:
#----   -DLAPACK_LIBS="c:\Temp\lapack.a;c:\Temp\blas.a"
#----
#---- If TRUE, then on return sets the variables:
#----   LAPACK_LIBS            - for linking
#----   LAPACK_ADD_LIBS        - for linking
#----   HAVE_LAPACK            - for header files
#----   HAVE_BLAS_F2C_WRAPPERS - for header files

#---- Set HOPSPACK_TEST_LAPACK_FUNCS=FALSE before including this file
#---- to bypass the calls to check_for_lapack_function below.  This
#---- macro fails to detect the presence of LAPACK functions within
#---- Fortran libraries on Windows.
IF (NOT DEFINED HOPSPACK_TEST_LAPACK_FUNCS) 
    SET(HOPSPACK_TEST_LAPACK_FUNCS TRUE) 
ENDIF () 


##########################################################################
#---- Macro to check for an LAPACK/BLAS function in each of the libraries
#---- specified by variable LAPACK_LIBS.  The macro tries appending
#---- an underscore and prefixing with "f2c_".
#---- On return sets the variables:
#----   HAVE_BLAS_F2C_WRAPPERS - TRUE or undefined
#----   LAPACK_ADD_LIBS        - may append another library
#----   HAVE_FN_${FNNAME}      - TRUE or FALSE

MACRO (check_for_lapack_function FNNAME)

    IF (NOT DEFINED LAPACK_LIBS)
        MESSAGE (FATAL_ERROR "HOPSPACK: Need to define an LAPACK/BLAS library.")
    ENDIF (NOT DEFINED LAPACK_LIBS)

    #---- First, try with the built-in utility CHECK_FUNCTION_EXISTS.
    #---- Unfortunately, this seems able to test only one library in isolation.
    INCLUDE (CheckFunctionExists)

    #---- Try appending an underscore.
    FOREACH (loopVar  ${LAPACK_LIBS})
        IF (NOT HAVEFN_${FNNAME})
            SET (CMAKE_REQUIRED_LIBRARIES ${loopVar})
            CHECK_FUNCTION_EXISTS (${FNNAME}_ HAVEFN_${FNNAME})
        ENDIF (NOT HAVEFN_${FNNAME})
    ENDFOREACH (loopVar)

    #---- Try prefixing f2c_ in case the library uses F2C wrappers.
    IF (NOT HAVEFN_${FNNAME})
        FOREACH (loopVar  ${LAPACK_LIBS})
            IF (NOT HAVEFN_${FNNAME})
                CHECK_FUNCTION_EXISTS (f2c_${FNNAME} HAVEFN_f2c_${FNNAME})

                IF (HAVEFN_f2c_${FNNAME})
                    SET (HAVEFN_${FNNAME} TRUE)
                    #---- The HOPSPACK header file needs to know about wrappers.
                    SET (HAVE_BLAS_F2C_WRAPPERS TRUE)
                ENDIF (HAVEFN_f2c_${FNNAME})
            ENDIF (NOT HAVEFN_${FNNAME})
        ENDFOREACH (loopVar)
    ENDIF (NOT HAVEFN_${FNNAME})

    #---- The CHECK_FUNCTION_EXISTS utility failed to find the function.

    #---- Try compiling directly, with all the libraries at once.
    #---- (Note that if LAPACK_LIBS is not a full path, then need to add
    #----    CMAKE_FLAGS -DLINK_DIRECTORIES:STRING=path )
    IF (NOT HAVEFN_${FNNAME})
        SET (ALL_LIB_LIST ${LAPACK_LIBS} ${LAPACK_ADD_LIBS})
        TRY_COMPILE (HAVEFN_${FNNAME}
                     ${HOPSPACK_BINARY_DIR}/test
                     ${HOPSPACK_SOURCE_DIR}/test/cmake_test_link_lapack.cpp
                     COMPILE_DEFINITIONS "-DTEST_${FNNAME}"
                     CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${ALL_LIB_LIST}"
                     OUTPUT_VARIABLE OUTPUT
                    )
        #---- To see output from the compiler, uncomment the following line:
        #MESSAGE (STATUS "Compiler Output " ${OUTPUT})

        IF (HAVEFN_${FNNAME})
            SET (detected_lapack_fortran TRUE)
            MESSAGE (STATUS "Looking for " ${FNNAME}_
                            " - found with try_compile")
        ENDIF (HAVEFN_${FNNAME})
    ENDIF (NOT HAVEFN_${FNNAME})

    #---- Try compiling directly, guessing LAPACK was built from Fortran
    #---- and needs Fortran symbols.
    IF (NOT HAVEFN_${FNNAME})
        SET (ALL_LIB_LIST ${LAPACK_LIBS} ${LAPACK_ADD_LIBS} g2c)
        TRY_COMPILE (HAVEFN_${FNNAME}
                     ${HOPSPACK_BINARY_DIR}/test
                     ${HOPSPACK_SOURCE_DIR}/test/cmake_test_link_lapack.cpp
                     COMPILE_DEFINITIONS "-DTEST_${FNNAME}"
                     CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${ALL_LIB_LIST}"
                     OUTPUT_VARIABLE OUTPUT
                    )
        #---- To see output from the compiler, uncomment the following line:
        #MESSAGE (STATUS "Compiler Output " ${OUTPUT})

        IF (HAVEFN_${FNNAME})
            SET (detected_lapack_fortran TRUE)
            MESSAGE (STATUS "Looking for " ${FNNAME}_
                            " - found with try_compile and g2c")
            SET (LAPACK_ADD_LIBS ${LAPACK_ADD_LIBS} g2c)
        ENDIF (HAVEFN_${FNNAME})
    ENDIF (NOT HAVEFN_${FNNAME})

    IF (NOT HAVEFN_${FNNAME})
        MESSAGE (STATUS "HOPSPACK: Cannot find '" ${FNNAME} "', looked in")
        FOREACH (loopVar ${LAPACK_LIBS})
            MESSAGE (STATUS "          " ${loopVar})
        ENDFOREACH (loopVar)
        FOREACH (loopVar ${LAPACK_ADD_LIBS})
            MESSAGE (STATUS "    using " ${loopVar})
        ENDFOREACH (loopVar)
        MESSAGE (FATAL_ERROR "HOPSPACK: Need an LAPACK/BLAS library.")
    ENDIF (NOT HAVEFN_${FNNAME})

ENDMACRO (check_for_lapack_function)
##########################################################################


#-------------------------------------------------------------------------
IF (NOT lapack)
#-------------------------------------------------------------------------
    MESSAGE (STATUS "HOPSPACK: No LAPACK, linear constraints disabled.")

#-------------------------------------------------------------------------
ELSE (NOT lapack)
#-------------------------------------------------------------------------
    MESSAGE (STATUS "HOPSPACK: LAPACK requested, linear constraints allowed.")

    #---- Enable EXACTLY ONE of the following options to link LAPACK/BLAS
    #---- libraries.  Enable by uncommenting, or passing on the command line.
    #---- See the HOPSPACK User Manual for details.

    #---- 1. Search the system library paths for an LAPACK library:
    FIND_LIBRARY (LAPACK_LIBS NAMES lapack
                  DOC "Location of LAPACK library")

    #---- 2. Let CMake search for a standard Fortran-based LAPACK library:
#    FIND_PACKAGE (LAPACK)

    #---- 3. Provide a specific location and name of LAPACK and/or BLAS:
    #----    The example below locates libraries installed and built in /tmp
    #----    with gcc on a Linux system.  The actual library files built
    #----    by the LAPACK distribution are named "libblas_LINUX.a"
    #----    and "liblapack_LINUX.a"
#    FIND_LIBRARY (LIB1 NAMES lapack_LINUX
#                  DOC "Location of LAPACK library"
#                  PATHS /tmp/lapack-3.1.1)
#    FIND_LIBRARY (LIB2 NAMES blas_LINUX
#                  DOC "Location of BLAS library"
#                  PATHS /tmp/lapack-3.1.1)
#    SET (LAPACK_LIBS ${LIB1} ${LIB2})

    #---- 4. Provide a specific location of a CLAPACK library on Windows.
    #----    The example below locates the library installed in c:\Temp.
#    FIND_LIBRARY (LAPACK_LIBS NAMES clapack
#                  DOC "Location of MSVC CLAPACK library"
#                  PATHS c:\\Temp\\CLAPACK3-Windows\\CLAPACK\\Release)


    #---- Now test for the presence of LAPACK functions used by HOPSPACK.
    IF (HOPSPACK_TEST_LAPACK_FUNCS)
        check_for_lapack_function (ddot)
        check_for_lapack_function (dgemv)
        check_for_lapack_function (dgemm)
#        THESE ARE SOMETIMES PROBLEMATIC; ASSUME BLAS CHECK IS GOOD ENOUGH.
#        check_for_lapack_function (dgesvd)
#        check_for_lapack_function (dgglse)
#        check_for_lapack_function (dgelss)
#        check_for_lapack_function (dgelqf)
    ENDIF ()

    MESSAGE (STATUS "HOPSPACK: Using LAPACK libraries:")
    FOREACH (loopVar ${LAPACK_LIBS})
        MESSAGE (STATUS "  " ${loopVar})
    ENDFOREACH (loopVar)
    FOREACH (loopVar ${LAPACK_ADD_LIBS})
        MESSAGE (STATUS "  " ${loopVar})
    ENDFOREACH (loopVar)

    SET (HAVE_LAPACK TRUE)

#-------------------------------------------------------------------------
ENDIF (NOT lapack)
#-------------------------------------------------------------------------
