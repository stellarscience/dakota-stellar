#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# CMakeLists.txt for util unit tests
#
include_directories(${Pecos_SOURCE_DIR}/util/src)
set(pecos_util_unit_tests
  linear_algebra_test.cpp
  linear_solvers_test.cpp
  math_tools_test.cpp
  )
add_executable(pecos_util_unit_tests 
  teuchos_unit_test_driver.cpp
  ${pecos_util_unit_tests}
  )

add_executable(options_unit_tests  options.cpp)
target_link_libraries(options_unit_tests pecos_util ${Teuchos_LIBRARIES})
add_test(options_unit_tests options_unit_tests)

target_link_libraries(pecos_util_unit_tests pecos_util ${Teuchos_LIBRARIES})
add_test(pecos_util_unit_tests pecos_util_unit_tests)
