#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# Dakota's modular surrogate model package

cmake_minimum_required(VERSION 2.8.9)

project(DakotaSurrogates CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ---
# Third-party software
# ---

# Teuchos is required
#find_package(teuchos)

include_directories(
  #${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS}
  ${Teuchos_INCLUDE_DIRS}
  ${Pecos_SOURCE_DIR}/util/src
  )


# ---
# Core libraries
# ---

# TODO: convenience function to build into a single library

# each depends only on previous directory

# core utilities for data structures and IO
#add_subdirectory(utils)

# math-specific utilities (Teuchos-related, boost numerics wrappers, etc.)

# disabled as moving to top-level
#add_subdirectory(math_utils)
#add_subdirectory(linear_solvers)

add_subdirectory(models)

add_subdirectory(pecos_wrapper)

#message("ind dir1 ${INCLUDE_DIRECTORIES}")
#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()


option(DAKOTA_ENABLE_DOCS "Enable documentation (requires doxygen)" OFF)
if( DAKOTA_ENABLE_DOCS )
    find_package( Doxygen )
    find_package( LATEX )
    if( DOXYGEN_FOUND )
        if ( LATEX_COMPILER )
	    message( STATUS "Latex was found" )
	else()
	    message( STATUS "Latex was not found. Doxygen will not be able to produce PDF documents" )
        endif()
        add_subdirectory( doc )
    else()
        message( "PyDakota documents cannot be built because doxygen could not be found.")
    endif()
endif()
