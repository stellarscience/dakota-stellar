#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# CMakeLists.txt for pecos source (WJB: Note the source file GLOB below)

include_directories(${Teuchos_INCLUDE_DIRS})

# TODO: when configuring each package, do we want to add to
# include_directories, or be explicit in each package?
include_directories(${Pecos_SOURCE_DIR}/util/src)

remove_definitions("-DHAVE_CONFIG_H")

option(BOOST_DISABLE_ASSERTS "Toggle BoostAssert, default is AssertDisabled" ON)
if (BOOST_DISABLE_ASSERTS)
  add_definitions("-DBOOST_DISABLE_ASSERTS")
endif(BOOST_DISABLE_ASSERTS)

if (${CMAKE_CXX_COMPILER_ID} MATCHES "XL" OR
    ${CMAKE_CXX_COMPILER_ID} MATCHES "VisualAge")
  add_definitions("-DBOOST_MPL_CFG_ASSERT_USE_RELATION_NAMES")
endif()


# --- Dependency configuration based on previously defined packages ---

if(HAVE_LHS)
  include_directories(${LHS_BINARY_DIR}/include)
endif()

if(HAVE_SPARSE_GRID)
  include_directories(${VPISparseGrid_SOURCE_DIR}/src)
endif()

if(HAVE_FFT)
  if(HAVE_FFTW)
    include_directories("${FFTW_SOURCE_DIR}/api")
    link_directories("${FFTW_BINARY_DIR}/.libs")
  endif()
endif()


## Pecos library sources.
file(GLOB Pecos_src *.cpp)
file(GLOB Pecos_hdr *.hpp)

add_library(pecos_src ${Pecos_src})
target_link_libraries(pecos_src pecos_util ${Pecos_PKG_LIBS} ${Pecos_TPL_LIBS}
		      Boost::boost)
if (HAVE_FFTW)
  add_dependencies(pecos_src fftw)
endif()
install(TARGETS pecos_src EXPORT ${ExportTarget} DESTINATION lib)
install(FILES ${Pecos_hdr} DESTINATION include)
