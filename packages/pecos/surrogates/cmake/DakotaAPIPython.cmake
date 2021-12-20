#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

macro(dakota_api_python)

  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})

  find_package(PythonLibs)
  include_directories(${PYTHON_INCLUDE_PATH})

endmacro()