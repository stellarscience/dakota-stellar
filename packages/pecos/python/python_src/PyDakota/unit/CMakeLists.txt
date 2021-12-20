#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# Simple test for SWIG interface (example from SWIG documentation)

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#SET(CMAKE_SWIG_FLAGS "")
# Helpful variables for later: CMAKE_SWIG_FLAGS CMAKE_SWIG_OUTDIR

set_source_files_properties(swig_example.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(swig_example.i PROPERTIES SWIG_FLAGS "-includeall")
swig_add_module(swig_example python swig_example.i swig_example.cpp)
# PYTHON_LIBRARIES?!?
swig_link_libraries(swig_example ${Python_LIBRARIES})

# TODO: make a build time copy with a dependency
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_swig_example.py
  ${CMAKE_CURRENT_BINARY_DIR}/test_swig_example.py COPY_ONLY)
add_test(NAME test_swig_example
  COMMAND "${Python_EXECUTABLE}" "test_swig_example.py"
          "$<TARGET_FILE_DIR:_swig_example>"
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_monomial.py
  ${CMAKE_CURRENT_BINARY_DIR}/test_monomial.py COPY_ONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_cross_validation.py
  ${CMAKE_CURRENT_BINARY_DIR}/test_cross_validation.py COPY_ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_compressed_sensing.py
  ${CMAKE_CURRENT_BINARY_DIR}/test_compressed_sensing.py COPY_ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/testTeuchos_ParameterList.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/testTeuchos_ParameterList.py COPY_ONLY)
