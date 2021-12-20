#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# CMakeLists.txt for pecos "unit" tests

include_directories(${Pecos_SOURCE_DIR}/src
                    ${Pecos_SOURCE_DIR}/util/src
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

# BMA: Historical rationale, prior to removing libpecos (empty.cxx). Leaving
# the explicit linkage to lhs to be safe, though removed libpecos. It seems
# just linking ${Pecos_LINK_LIBRARIES} would suffice.
#
# Ideally the test executables need to link to only "pecos" because it
# already links to its dependencies.  However, the lhs library defines
# at least one symbol (SSPEV) that is also provided by lapack.  The
# lhs library must come before lapack even though lhs does not link to
# lapack.  This ensures that the linker prefers symbols from lhs over
# conflicting symbols in lapack.
set(pecos_test_libs pecos_src lhs)

# Complete these for each specific test.
add_executable(pecos_ifft_g pecos_ifft_g.cpp)
target_link_libraries(pecos_ifft_g ${pecos_test_libs})
add_test(pecos_ifft_g pecos_ifft_g)

# Complete these for each specific test.
add_executable(pecos_ifft_sd pecos_ifft_sd.cpp)
target_link_libraries(pecos_ifft_sd ${pecos_test_libs})
add_test(pecos_ifft_sd pecos_ifft_sd)
set_tests_properties(pecos_ifft_sd PROPERTIES DEPENDS pecos_ifft_g)

# Complete these for each specific test.
add_executable(boost_test_dist boost_test_dist.cpp)
target_link_libraries(boost_test_dist ${pecos_test_libs})
add_test(boost_test_dist boost_test_dist)

# Complete these for each specific test.
add_executable(boost_test_rng boost_test_rng.cpp)
target_link_libraries(boost_test_rng ${pecos_test_libs})
add_test(boost_test_rng boost_test_rng)

# Complete these for each specific test.
add_executable(pecos_int_driver pecos_int_driver.cpp)
target_link_libraries(pecos_int_driver ${pecos_test_libs})
add_test(pecos_int_driver pecos_int_driver)

# Complete these for each specific test.
add_executable(pecos_gsg_driver pecos_gsg_driver.cpp)
target_link_libraries(pecos_gsg_driver ${pecos_test_libs})
add_test(pecos_gsg_driver pecos_gsg_driver)

# Complete these for each specific test.
add_executable(pecos_lhs_driver pecos_lhs_driver.cpp)
target_link_libraries(pecos_lhs_driver ${pecos_test_libs})
add_test(pecos_lhs_driver pecos_lhs_driver)

set(pecos_teuchos_unit_tests_src
  pecos_pochhammer.cpp
  pecos_discrete_poly.cpp
  pecos_linear_solvers.cpp
  )
add_executable(pecos_teuchos_unit_tests teuchos_unit_test_driver.cpp ${pecos_teuchos_unit_tests_src} )
target_link_libraries(pecos_teuchos_unit_tests ${pecos_test_libs} ${Teuchos_LIBRARIES})
add_test(pecos_teuchos_unit_tests pecos_teuchos_unit_tests)

install(TARGETS pecos_gsg_driver DESTINATION bin)
