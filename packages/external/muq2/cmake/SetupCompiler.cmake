

include(CheckCXXCompilerFlag)

if(MUQ_USE_MPI)

  find_package(MPI REQUIRED)

  list(APPEND MUQ_LINK_LIBS ${MPI_CXX_LIBRARIES})
  list(APPEND MUQ_EXTERNAL_INCLUDES ${MPI_CXX_INCLUDE_DIRS})

  include_directories(${MPI_CXX_INCLUDE_DIRS})
  link_directories(${MPI_CXX_LIBRARIES})

  set(MUQ_HAS_MPI 1)

else(MUQ_USE_MPI)
  set(MUQ_HAS_MPI 0)
endif(MUQ_USE_MPI)

set(CMAKE_CXX_FLAGS_DEBUG  "-O0")
set(CMAKE_CXX_FLAGS_RELEASE  "-O3")

# default to a release build
message(STATUS "User defined build type = " ${CMAKE_BUILD_TYPE})
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()
message(STATUS "Final build type = " ${CMAKE_BUILD_TYPE})

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check to see if const& and by value need to be treated separately in AnyConstCast
INCLUDE(AnyCastCheck)

IF(MUQ_USE_OPENMP)
	CHECK_CXX_COMPILER_FLAG("-fopenmp" HAS_FOPENMP)
	CHECK_CXX_COMPILER_FLAG("-pthread" HAS_PTHREAD)

	if(HAS_FOPENMP AND HAS_PTHREAD)
  	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -pthread -ldl")
   else()
	   if(HAS_FOPENMP)
	     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -ldl")
	  else()
		  message(WARNING "The flag MUQ_USE_OPENMP is ON, but the compiler does not seem to support the -fopenmp flag.  OPENMP will not be used.")
	  endif()
  endif()
ENDIF(MUQ_USE_OPENMP)


# this is required for cmake version 3.0.0 and later
if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
endif()

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")
