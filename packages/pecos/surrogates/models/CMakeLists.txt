#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

add_subdirectory(src)
if(BUILD_TESTING AND PECOS_ENABLE_TESTS)
  add_subdirectory(unit)
endif()
