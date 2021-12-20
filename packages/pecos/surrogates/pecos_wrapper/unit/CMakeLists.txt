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
                    ${DakotaSurrogates_SOURCE_DIR}/pecos_wrapper/src
                    ${Pecos_SOURCE_DIR}/src
                    ${Boost_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIRS}
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

set(pce_approx_test
  test_polynomial_chaos.cpp
)

add_executable(pce_approx_test teuchos_unit_test_driver.cpp ${pce_approx_test} )
target_link_libraries(pce_approx_test
  pecos_wrapper pecos_src pecos_models pecos_util
  ${Pecos_PKG_LIBS} ${Pecos_TPL_LIBS}
  )


add_test(pce_approx_test pce_approx_test)
