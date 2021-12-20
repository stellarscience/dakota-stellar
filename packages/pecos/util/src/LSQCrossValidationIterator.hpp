/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_LSQ_CROSS_VALIDATION_ITERATOR_HPP
#define PECOS_UTIL_LSQ_CROSS_VALIDATION_ITERATOR_HPP

#include "LinearSystemCrossValidationIterator.hpp" 

namespace Pecos {
namespace util {

 /**
  * \class LSQCrossValidationIterator
  * \brief Peform k-folds cross validation for over-determined least squares
  * problems using a fast analytical formula.
  */
class LSQCrossValidationIterator : public LinearSystemCrossValidationIteratorBase{
protected:
  /// Coefficients x obtained by solving Ax=b using least squares regression
  /// on entire data set
  /// matrix (num_coeffs x num_rhs)
  RealMatrix coeffs_;

  /// The l2 norm of the residuals B-Ax obtained by solving Ax=b using
  /// least squares regression on entire data set
  /// matrix (num_coeffs x num_rhs)
  RealVector residualNorms_;

public:

  /**\brief Compute the leave K out residuals on the validation points of each 
   * fold using an analytical formula.
   *
   * The analytical formula requires, for rach fold, 
   * the inversion of a small matrix with the
   * size the size of the number rows of the subset of A formed on the 
   * validation indices. If number of folds is K and number or rows of A is M
   * then size of matrix to invert is (M/K x M/K)
   *
   * \param[in] A matrix (num_rows x num_coeffs)
   * \param[in] B matrix (num_coeffs x num_rhs)
   * \param[in] residuals matrix (num_rows x num_rhs) residuals=B-AX
   *
   * \param[out] result list of num_folds matrices (1 x num_rhs) which contain 
   *             the residual at the left out validation point for each fold
   */
  void lsq_lko_cross_validation(const RealMatrix &A, const RealMatrix &AtA_inv,
				const RealMatrix &residuals,
				std::vector<RealMatrix> &result);

  
  /**\brief Compute the leave one out residuals on the validation point of each 
   * fold using an analytical formula.
   *
   * This is a faster specialization of leave k out cross validation for k=1.
   *
   * \param[in] A matrix (num_rows x num_coeffs)
   * \param[in] B matrix (num_coeffs x num_rhs)
   * \param[in] residuals matrix (num_rows x num_rhs) residuals=B-AX
   *
   * \param[out] result list of num_folds matrices (1 x num_rhs) which contain 
   *             the residual at the left out validation point for each fold
   */
  void lsq_loo_cross_validation(const RealMatrix &A, const RealMatrix &AtA_inv,
				const RealMatrix &residuals, 
				std::vector<RealMatrix> &result);
  
  /**\brief Run K fold cross validation.
   *
   * Obtain cross validation scores for least squares solutions X of AX=B.
   * Solutions X are calculated as a by product of cross validation and are
   * stored internally and can be retrieved using get_coefficients()
   * 
   * \param[in] A matrix (num_rows x num_coeffs)
   * \param[in] B matrix (num_coeffs x num_rhs)
   * \param[in] opts list of options needed for cross validation
   */
  void run(const RealMatrix &A, const RealMatrix &B,
	   OptionsList& opts);

  
  /**\brief Get the coefficients x obtained by solving Ax=b 
   * using least squares regression. 
   *
   * The coefficients are obtained as a by product of producing the cross
   * validation estimates.
   *
   * \param[out] result coefficients x matrix (num_coeffs x num_rhs)
   */
  void get_coefficients(RealMatrix &result) const;

  /**\copydoc LinearSystemCrossValidationIteratorBase::get_best_scores()*/
  void get_best_scores(RealVector &result) const;

  /**\copydoc LinearSystemCrossValidationIteratorBase::generate_best_solutions()*/
  void generate_best_solutions(const RealMatrix &A, const RealMatrix &B,
			       RealMatrix &best_solutions,
			       RealVector &residuals,
			       OptionsList & opts);
};

std::shared_ptr<LSQCrossValidationIterator> cast_to_least_squares_cross_validation_iterator(std::shared_ptr<LinearSystemCrossValidationIteratorBase> &solver);


}  // namespace util
}  // namespace Pecos

#endif  // include guard
