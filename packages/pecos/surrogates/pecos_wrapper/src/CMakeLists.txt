#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

SET(SOURCES PolynomialChaosExpansionWrap.cpp
            PolynomialChaosExpansion.cpp
            OrthogonalPolynomialBasis.cpp)
SET(HEADERS PolynomialChaosExpansionWrap.hpp
            PolynomialChaosExpansion.hpp
            OrthogonalPolynomialBasis.hpp)
add_library(pecos_wrapper ${SOURCES} ${HEADERS})
include_directories(${DakotaSurrogates_SOURCE_DIR}/pecos_wrapper/src
		    ${DakotaSurrogates_SOURCE_DIR}/models/src
		    ${Pecos_SOURCE_DIR}/src
		    ${Teuchos_INCLUDE_DIRS})
if(HAVE_SPARSE_GRID)
  include_directories(${VPISparseGrid_SOURCE_DIR}/src)
  # BMA: Is this needed?
  list(APPEND Pecos_libs sparsegrid)
endif()


target_link_libraries(pecos_wrapper pecos_src pecos_models pecos_util
  ${Pecos_PKG_LIBS} ${Pecos_TPL_LIBS})
