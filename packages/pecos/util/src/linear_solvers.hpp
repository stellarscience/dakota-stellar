/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_LINEAR_SOLVERS_HPP
#define PECOS_UTIL_LINEAR_SOLVERS_HPP

#include "LinearSystemSolver.hpp"
#include "OMPSolver.hpp"
#include "LARSolver.hpp"
#include "LSQSolver.hpp"
#include "EqConstrainedLSQSolver.hpp"
#include "CrossValidatedSolver.hpp"

namespace Pecos {
namespace util {

/**\
 * \brief Initialize a linear system solver from a list of options and return 
 * a shared pointer to the base class.
 *
 * \param[in] opts a list of options
 *
 * opts (required parameters)
 * -------------------------
 * 
 * "regression_type" : RegressionType
 *     The regression-type. Accepted values are SVD_LEAST_SQ_REGRESSION, 
 *     QR_LEAST_SQ_REGRESSION, LU_LEAST_SQ_REGRESSION, ORTHOG_MATCH_PURSUIT
 *     LEAST_ANGLE_REGRESSION, LASSO_REGRESSION, EQ_CONS_LEAST_SQ_REGRESSION
 *
 * opts (optional parameters)
 * -------------------------
 * 
 * "use-cross-validation" : boolean default=false
 *     If true return a CrossValidatedSolver which wraps a standard linear 
 *     solver of "regression_type"
 */
std::shared_ptr<LinearSystemSolver> regression_solver_factory(OptionsList &opts);

std::shared_ptr<OMPSolver> cast_linear_system_solver_to_ompsolver(std::shared_ptr<LinearSystemSolver> &solver);

std::shared_ptr<LARSolver> cast_linear_system_solver_to_larssolver(std::shared_ptr<LinearSystemSolver> &solver);

std::shared_ptr<LSQSolver> cast_linear_system_solver_to_lsqsolver(std::shared_ptr<LinearSystemSolver> &solver);

std::shared_ptr<EqConstrainedLSQSolver> cast_linear_system_solver_to_eqconstrainedlsqsolver(std::shared_ptr<LinearSystemSolver> &solver);
  

}  // namespace util
}  // namespace Pecos

#endif  // include guard
