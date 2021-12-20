#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# CMakeLists.txt for linear_solvers "unit" tests
message("${DakotaSurrogates_SOURCE_DIR}")

include_directories(${DakotaSurrogates_SOURCE_DIR}/models/src
                    ${Teuchos_INCLUDE_DIRS})

remove_definitions("-DHAVE_CONFIG_H")

option(BOOST_DISABLE_ASSERTS "Toggle BoostAssert, default is AssertDisabled" ON)
if (BOOST_DISABLE_ASSERTS)
  add_definitions("-DBOOST_DISABLE_ASSERTS")
endif(BOOST_DISABLE_ASSERTS)

if (${CMAKE_CXX_COMPILER_ID} MATCHES "XL" OR
    ${CMAKE_CXX_COMPILER_ID} MATCHES "VisualAge")
  add_definitions("-DBOOST_MPL_CFG_ASSERT_USE_RELATION_NAMES")
endif()

if(HAVE_SPARSE_GRID)
  include_directories(${VPISparseGrid_SOURCE_DIR}/src)
endif(HAVE_SPARSE_GRID)

set(polynomial_approx_test
  polynomial_approximation.cpp
)

set(directANN_test
  directANN_approximation.cpp
)

add_executable(polynomial_approx_test teuchos_unit_test_driver.cpp ${polynomial_approx_test} )
target_link_libraries(polynomial_approx_test pecos_models pecos_util ${Pecos_TPL_LIBS})

add_executable(directANN_test teuchos_unit_test_driver.cpp ${directANN_test} )
target_link_libraries(directANN_test pecos_models pecos_util ${Pecos_TPL_LIBS})

add_test(polynomial_approx_test polynomial_approx_test)
add_test(directANN_test directANN_test)
