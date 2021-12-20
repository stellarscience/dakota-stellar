#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

SET(SOURCES Function.cpp CppFunction.cpp
	    Approximation.cpp	
	    PolyApproximation.cpp Monomial.cpp DirectANN.cpp
            Variables.cpp BoundedVariables.cpp
	    VariableTransformation.cpp AffineVariableTransformation.cpp
	    RegressionBuilder.cpp
      DataScaler.cpp)
SET(HEADERS Function.hpp CppFunction.hpp
	    Approximation.hpp
	    PolyApproximation.hpp Monomial.hpp DirectANN.hpp
            Variables.hpp BoundedVariables.hpp
	    VariableTransformation.hpp AffineVariableTransformation.hpp
	    RegressionBuilder.hpp
      DataScaler.hpp)
add_library(pecos_models ${SOURCES} ${HEADERS})
include_directories(${DakotaSurrogates_SOURCE_DIR}/models/src
	            ${Teuchos_INCLUDE_DIRS})

target_link_libraries(pecos_models pecos_util ${Pecos_TPL_LIBS} Boost::boost)
