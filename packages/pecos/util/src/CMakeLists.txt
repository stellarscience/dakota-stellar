#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# Pecos utilities that should depend only on STL and Teuchos

set(UTIL_SOURCES
  least_angle_regression.cpp
  linear_algebra.cpp
  linear_solvers.cpp
  math_tools.cpp
  sparse_matrices.cpp
  CrossValidationIterator.cpp
  LinearSystemCrossValidationIterator.cpp
  LSQCrossValidationIterator.cpp
  CrossValidatedSolver.cpp
)

set(UTIL_HEADERS
  CrossValidatedSolver.hpp
  CrossValidationIterator.hpp
  EqConstrainedLSQSolver.hpp
  LARSolver.hpp
  least_angle_regression.hpp
  linear_algebra.hpp
  linear_solvers.hpp
  LinearSystemCrossValidationIterator.hpp
  LinearSystemSolver.hpp
  LSQCrossValidationIterator.hpp
  LSQSolver.hpp
  math_tools.hpp
  OMPSolver.hpp
  OptionsList.hpp
  orthogonal_matching_pursuit.hpp
  sparse_matrices.hpp
  teuchos_data_types.hpp
  )

add_library(pecos_util ${UTIL_SOURCES} ${UTIL_HEADERS})
#SET(CMAKE_INCLUDE_CURRENT_DIR ON)
target_link_libraries(pecos_util ${Pecos_TPL_LIBS} Boost::boost)
install(TARGETS pecos_util EXPORT ${ExportTarget} DESTINATION lib)
install(FILES ${UTIL_HEADERS} DESTINATION include)
