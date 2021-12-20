/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_LINEAR_SYSTEM_CROSS_VALIDATION_ITERATOR_HPP
#define PECOS_UTIL_LINEAR_SYSTEM_CROSS_VALIDATION_ITERATOR_HPP

#include "CrossValidationIterator.hpp"
#include "LinearSystemSolver.hpp"
#include <memory>

namespace Pecos {
namespace util {

/**
 * \class LinearSystemCrossValidationIteratorBase
 * \brief Base class for classes peforming k-folds cross validation of
 * solutions obtained using linear system solvers.
 */
class LinearSystemCrossValidationIteratorBase : public CrossValidationIterator {
protected:

  /// The cross validation scores at each of the unique tolerances and RHS
  /// List of num RHS vectors of size (num_unique_tolsx1)
  std::vector<RealVector> scores_;

  /// The difference between the validation value and the approximation 
  /// at each validation point and for each solution generated on the 
  /// solution path for last RHS column accessed
  std::vector< RealMatrix > foldDiffs_;

  /// Set the cross validation parameters specified in opts 
  void parse_options(const RealMatrix &A, const RealMatrix &B,
		     OptionsList &opts);

public:

  /// Default Constructor
  LinearSystemCrossValidationIteratorBase();

  /// Destructor
  virtual ~LinearSystemCrossValidationIteratorBase();

  /**
   * \brief Set the number of equations per point
   * 
   * If we are cross validating a linear system used to build a polynomial 
   * surrogate, the number equations per point will be one if just using 
   * function values. However if we are also using derivatives then the number
   * of equations per point will be num_vars+1

   * \param num_eq - The number of equations per point
   */
  void set_num_equations_per_point(int num_eq);

  /**
   * \brief Extract a subset of rows from a matrix.
   * 
   * \param A matrix of size (num_rows x num_cols) - 
   *     A set of equations for each point being used in cross validation.
   *	 This function assumes that there is a set of primary equations stored 
   *      0...numPts_-1, then the non-primary equations are stored 
   *	  1...numEquationsPerPoint_-1 for first point then all non-primary 
   *	  equations for second point and so on.
   * \param indices vector of length (l) - 
   *     The chosen indices that define the chosen rows of A 
   * 
   * \return result matrix of size (l*numEquationsPerPoint_ x num_cols) - 
   *     Subset of rows of the matrix A
   */
  void extract_matrix(const RealMatrix &A,
		      const IntVector &indices, 
		      RealMatrix &result) const;

  /**
   * \brief Extract the subset of rows of a matrix A and the corresponding rows
   * from the vector b.
   * 
   * \param A matrix of size (numPts_*numEquationsPerPoint_ x num_cols) - 
   *     A set of equations for each point being used in cross validation.
   *     This function assumes that there is a set of primary equations stored 
   *     0...numPts_-1, then the non-primary equations are stored 
   *	 1...numEquationsPerPoint_-1 for first point then all non-primary 
   *	 equations for second point and so on.
   * \param vector of size (numPts_*numEquationsPerPoint_) - 
   *     A set of data for each point being used in cross validation.
   *     This function assumes that there is a set of primary equations stored 
   *     0...numPts_-1, then the non-primary equations are stored 
   *	 1...numEquationsPerPoint_-1 for first point then all non-primary 
   *	 equations for second point and so on.
   * \param indices vector of length (l) - 
   *     The chosen indices that define the chosen rows of A 
   * 
   * \return result_0 matrix of size (l*numEquationsPerPoint_ x num_cols) - 
   *     Subset of rows of the matrix A
   * \return result_1 vector of size (l*numEquationsPerPoint_) - 
   *     Subset of rows of the vector b
   */
  void extract_linear_system(const RealMatrix &A, const RealVector &b,
			     const IntVector &indices, 
			     RealMatrix &result_0,
			     RealVector &result_1) const;
  /**
   * \brief Run the cross validation.
   */
  virtual void run(const RealMatrix &A, const RealMatrix &B,
		   OptionsList& opts)=0;

  /**\brief Get the cross validation scores for each tolerance considered
   * in cross validation
   */
  void get_scores(std::vector<RealVector> &result) const;

  /**\brief
   */
  void get_fold_validation_residuals(std::vector< RealMatrix > &result) const;

  /**\brief Get the best cross validation scores for each RHS
   * \param[out] result vector (num_rhs x 1)
   *    The best cross validation scores for each RHS 
   */
  void get_best_scores(RealVector &result) const;

  /**\brief Get the indices in scores_ of the best cross validation scores 
   * for each RHS
   * \param[out] scores vector (num_rhs x 1)
   *    The index of the best cross validation scores for each RHS 
   */
  void get_best_score_indices(IntVector & result) const;

  
  /**\brief Generate the solutions using the best regularization parameters
   * as chosen by cross validation.
   * \param[in] A matrix (num_rows x num_cols)
   *     The matrix A
   * \param[in] B vector (num_cols x num_rhs)
   *     The RHS B
   * \param[inout] opts  
   *     List of options. On exit opts will contain the best residual tolerances
   * \param[out] best_solutions matrix (num_cols x num_rhs)
   *     The best solutions for each RHS
   * \param[out] best_solutions vector (num_rhs x 1)
   *     The l2 norm of the resdidual for the best solution for each RHS
   */
  virtual void generate_best_solutions(const RealMatrix &A, const RealMatrix &B,
				       RealMatrix &best_solutions,
				       RealVector & residuals,
				       OptionsList & opts)=0;
};

  
/**
 * \class LinearSystemCrossValidationIterator
 * \brief Peform k-folds cross validation of solutions obtained 
 * using linear system solvers.
 */
class LinearSystemCrossValidationIterator : public LinearSystemCrossValidationIteratorBase {
protected:
  
  /// The maximum number of unique tolerances for which cross validation scores
  /// will be computed
  int maxNumUniqueTols_;

  /// The unique tolerances for each RHS
  /// List of num RHS vectors of size (num_unique_tolsx1)
  std::vector<RealVector> uniqueTols_;

  
  /// The solver tolerances corresponding to each solution generated on the 
  /// solution path for last RHS column accessed. List of num RHS vectors of
  /// size (num_unique_toleraces x 1)
  std::vector< RealVector > foldTols_;

  /// The cross validation scores ( a function of foldDiffs_ ) for each fold
  /// and solution generated on solution path for last RHS column accessed
  std::vector< RealVector > foldScores_;

  
  /// Pointer to the solver for solving the linear system
  std::shared_ptr<LinearSystemSolver> linearSystemSolver_;

  
  /**
   * \brief Once cross validation has been run, compute the cross validation
   * scores.
   */
  void compute_scores(RealVector &scores, RealVector &unique_tols);

  /**
   * \brief Return a subset of the tolerances (used to terminate method for 
   * solving linear system) computed during cross validation.

   * There are often thousands of unique parameter values so take 
   * a thinned subset of these. The size of the subset is controlled by
   * max_num_unique_tols
   */
  void define_unique_tolerances(RealVector &unique_tols);

  
  /**
   * \brief Run the cross validation for a single RHS
   */
  void run_single_rhs(const RealMatrix &A, const RealVector &b,
		      OptionsList &opts);
  
public:

  /// Default Constructor
  LinearSystemCrossValidationIterator();

  /// Destructor
  ~LinearSystemCrossValidationIterator();
  
  /**
   * \brief Set the maximum number of unique tolerances for which cross 
   * validation scores will be computed. See define_unique_tolerances.
   */
  void set_max_num_unique_tolerances(int max_num_tols);
  
  
  /**\brief
   */
  void get_fold_tolerances(std::vector< RealVector > &result) const;


  /**\brief
   */
  void get_fold_scores(std::vector< RealVector > &result) const;

  /**\brief Get the unique tolerances for which the cross validation scores in
     get_scores() correspond.
   */
  void get_unique_tolerances(std::vector<RealVector> &result) const;

  /**
   * \brief Set the linear system solver which will be used to solve Ax~b.
   *
   * The solver must have a member function with the API
   * solve_single_rhs(A,b,x,metrics,opts);
   */
  void set_linear_system_solver(const std::shared_ptr<LinearSystemSolver> &solver);
  
  /**
   * \brief Run the cross validation.
   */
  void run(const RealMatrix &A, const RealMatrix &B,
	   OptionsList& opts);

  /**\brief Get the l2 norm of the residuals associated with the solutions
     with the best cross validation scores for each RHS
   *
   * Needs to return a matrix (not a vector) because of use with 
   * OptionsList where options assumes a matrix not a a vector
   *
   * \param[out] result matrix (num_rhs x 1)
   *    The best cross validation scores for each RHS 
   */
  void get_best_residual_tolerances(RealMatrix &result) const;

  /**\brief Get a shared pointer to the linear system solver being used 
   * by the cross validation
   */
  std::shared_ptr<LinearSystemSolver> get_linear_system_solver();

  /**\copydoc LinearSystemCrossValidationIteratorBase::generate_best_solutions()*/
  void generate_best_solutions(const RealMatrix &A, const RealMatrix &B,
			       RealMatrix &result_0,
			       RealVector &result_1,
			       OptionsList & opts);

  /** \brief Get the adjusted best residual tolerances which will
   * be used to find the solution on the entire data set.
   * We adjust residual tolerances to account for fact they were found
   * using a subset of the data. We assume the squared l2 residual_norm
   * on the training data grows linearly with the size of the training data.
   *
   * Needs to return a matrix (not a vector) because of use with 
   * OptionsList where options assumes a matrix not a a vector
   *
   * \param[out] result matrix (num_rhs x 1)
   *    The best cross validation scores for each RHS adjusted by multiplying by
   *    num_folds/(num-folds-1)
   */
  void get_adjusted_best_residual_tolerances(RealMatrix &result) const;
 
};

std::shared_ptr<LinearSystemCrossValidationIterator> cast_to_linear_system_cross_validation_iterator(std::shared_ptr<LinearSystemCrossValidationIteratorBase> &solver);

}  // namespace util
}  // namespace Pecos

#endif  // include guard
