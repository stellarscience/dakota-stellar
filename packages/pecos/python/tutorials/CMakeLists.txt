#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Multivariate Function Approximation\ With Dakota.ipynb.in" "${CMAKE_CURRENT_BINARY_DIR}/Multivariate Function Approximation\ With Dakota.ipynb" COPY_ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_notebooks.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/test_notebooks.py COPY_ONLY)
