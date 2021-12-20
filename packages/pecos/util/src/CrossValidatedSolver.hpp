/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_CROSS_VALIDATED_SOLVER_HPP
#define PECOS_UTIL_CROSS_VALIDATED_SOLVER_HPP

#include "LinearSystemSolver.hpp"
#include "LinearSystemCrossValidationIterator.hpp"
#include <memory>

namespace Pecos {
namespace util {


/**
 * \class CrossValidatedSolver
 * \brief A wrapper for generating cross validated solutions to lineaer 
 * systems. 
 */
class CrossValidatedSolver : public LinearSystemSolver{
protected:
  /// The solutions for each rhs. matrix (num_coeff x num_rhs)
  RealMatrix solutions_;

  /// The l2 norm of the residuals for each rhs. matrix (num_coeff x num_rhs)
  RealVector residualNorms_;
  
  /// The cross validation iterator
  std::shared_ptr<LinearSystemCrossValidationIteratorBase> cvIterator_;

  /**\brief Generate the solutions corresponding to the best residual
   * tolerances.
   * \param[out] matrix  (num_coeff x num_rhs)
   *     The best solutions for each rhs
   */
  void generate_best_solutions(const RealMatrix &A, const RealMatrix &B,
			       RealMatrix &best_solutions,
			       OptionsList & opts);

  /**\copydoc LinearSystemSolver::unnormalize_coefficients()*/
  void unnormalize_coefficients(const RealVector &column_norms);
  
public:
  /// Default constructor
  CrossValidatedSolver();
  
  /// Destructor
  ~CrossValidatedSolver();

  /* \brief Specify the type of regression solver to use.
   *
   * This function creates a derived instance of a 
   * LinearSystemCrossValidationIteratorBase with a call to
   * linear_system_cross_validation_iterator_factory(opts) The returned
   * object will have a means to solve Ax=B.
   */
  void set_linear_system_solver(RegressionType regression_type);

  /**
   * \brief Find solutions to \f$AX \approx B\f$ 
   * The solutions are stored internally and not returned.
   *
   * \param[in] A matrix (num_rows x num_cols)
   *     The matrix A
   * \param[in] B vector (num_cols x num_rhs)
   *     The rhs B
   * \param[in] opts  
   *     List of options
   *
   * LinearSystemSolver::solve preconditions normalizes matrix
   * We do not want to do this here because it will be done by the
   * linear solver this class contains
   */
  void solve(const RealMatrix &A, const RealMatrix &B, OptionsList & opts );


  /**\brief
   *
   * opts (required parameters)
   * --------------------------
   *
   * "cv-opts" : OptionsList
   *     The cross validation parameters. See documentation of cross validation
   *     iterators
   *
   * opts (optional parameters)
   * --------------------------
   *
   * Parameters associated with running the linear system solvers.
   * See documentation of linear system solver, se.g OMPSolver,  
   * or if requesting a form of least squares, e.g. 
   * SVD_LEAST_SQ_REGRESSION, QR_LEAST_SQ_REGRESSION, LU_LEAST_SQ_REGRESSION 
   * then refer to the documentation of the LSQCrossValidationIterator
   */ 
  void multi_rhs_solve(const RealMatrix &A, const RealMatrix &B,
		       OptionsList & opts);

  /**\copydoc LinearSystemSolver::get_solutions_for_all_regularization_params()*/ 
  void get_solutions_for_all_regularization_params(RealMatrix &result_0,
						   int rhs_num) const;

  /**\copydoc LinearSystemSolver::get_residuals_for_all_regularization_params()*/ 
  void get_residuals_for_all_regularization_params(
      RealVector &result_0, int rhs_num) const;

  /**\copydoc LinearSystemSolver::get_final_solutions()*/
  void get_final_solutions(RealMatrix &result_0) const;

    /**\copydoc LinearSystemSolver::get_final_residuals()*/
  void get_final_residuals(RealVector &result_0) const;

  /**\brief return the best cross validation scores for each RHS
   * \param[out] result vector (num_rhs x 1)
   *     The best cross validation scores for each RHS
   */
  void get_best_scores(RealVector &result) const;

  /// Return a shared pointer to the cross validation iterator
  std::shared_ptr<LinearSystemCrossValidationIteratorBase> get_cross_validation_iterator();
};

std::shared_ptr<LinearSystemCrossValidationIteratorBase> linear_system_cross_validation_iterator_factory(OptionsList &opts);
  
/// Cast from LinearSystemSolver to a CrossValidatedSolver
std::shared_ptr<CrossValidatedSolver> cast_to_cross_validated_solver(std::shared_ptr<LinearSystemSolver> &solver);
 
}  // namespace util
}  // namespace Pecos

#endif  // include guard
