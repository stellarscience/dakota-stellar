/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_EQ_CONSTRAINED_LSQ_SOLVER_HPP
#define PECOS_UTIL_EQ_CONSTRAINED_LSQ_SOLVER_HPP

#include "LinearSystemSolver.hpp"

namespace Pecos {
namespace util {

/**
 *\class EqConstrainedLSQSolver
 *Solve the linear equality constrained least squares problem
 * \f$ \min_x \lVert Ax - b \rVert\f$ subject to \f$ Bx = d \f$
 * where A is an m-by-n matrix, B is a p-by-n matrix, b is an m vector, 
 * and d is a p vector, with $p \leq n \leq p+m$.   
 */
class EqConstrainedLSQSolver : public LinearSystemSolver{
protected:
  std::vector<RealVector> solutions_;
  RealVector residuals_;

  /**\copydoc LinearSystemSolver::unnormalize_coefficients()*/
  void unnormalize_coefficients(const RealVector &column_norms){
    for (size_t i=0; i<solutions_.size(); ++i)
      adjust_coefficients(column_norms, solutions_[i]);
  };
  
public:
  
  /// Default constructor
  EqConstrainedLSQSolver(): LinearSystemSolver() {};
  /// Destructor
  ~EqConstrainedLSQSolver(){};

  /**\copydoc LinearSystemSolver::get_solutions_for_all_regularization_params()*/
  void get_solutions_for_all_regularization_params(
      RealMatrix &result_0, int rhs_num) const{
    shape_uninitialized(result_0, solutions_[rhs_num].numRows(), 1);
    for (int i=0; i<solutions_[rhs_num].numRows(); ++i)
      result_0(i,0)=solutions_[rhs_num][i];
  };

  /**\copydoc LinearSystemSolver::get_residuals_for_all_regularization_params()*/ 
  void get_residuals_for_all_regularization_params(
      RealVector &result_0, int rhs_num) const{
    size_uninitialized(result_0, 1);
    result_0[0]=residuals_[rhs_num];
  };

  /**\copydoc LinearSystemSolver::get_final_solutions()*/
  void get_final_solutions(RealMatrix &result_0) const{
    shape_uninitialized(result_0, solutions_[0].length(), (int)solutions_.size());
    for (size_t j=0; j<solutions_.size(); ++j){
      for (int i=0; i<solutions_[j].length(); ++i)
	result_0(i,j) = solutions_[j][i];
    }
  };

  /**\copydoc LinearSystemSolver::get_final_residuals()*/  
  void get_final_residuals(RealVector &result_0) const{
    result_0 = residuals_;
  };

  /**\copydoc LinearSystemSolver::multi_rhs_solve()*/  
  void multi_rhs_solve(const RealMatrix &A, const RealMatrix &B,
		       OptionsList & opts){
    RealVector b;
    size_uninitialized(residuals_, B.numCols());
    if( solutions_.empty() )
      solutions_.resize(B.numCols());
    for (int i=0; i<B.numCols(); ++i){
      b = Teuchos::getCol(Teuchos::View, const_cast<RealMatrix &>(B), i);
      residuals_[i]=single_rhs_solve(A, b, opts, solutions_[i]);
    }
  };

  /**\brief  Solve the linear equality constrained least squares problem
   * \f$ \min_x \Vert Ax - b \Vert\f$ subject to \f$Bx = d \f$
   * where A is an m-by-n matrix, B is a p-by-n matrix, b is an m vector, 
   * and d is a p vector, with $p \leq n \leq p+m$. 
   *
   * \param[in] A matrix (m+p x n)
   *    The A and B matrices. A is stored in first rows then B is stored in 
   *    remaining rows.
   * \param[in] b vector (m+p x 1)
   *    The b and d vectors. b is stored in first rows then d is stored in 
   *    remaining rows.
   * \param[in] opts  
   *     List of options
   *
   * \param[out] solution vector (n x 1)
   *    The equality constrined solution
   *
   * opts (required parameters)
   * --------------------------
   * 
   * "num-primary-equations" : integer 
   *     The number of rows (p) in the vector d
   * 
   */
  Real single_rhs_solve(const RealMatrix &A, const RealVector &b,
			OptionsList& opts, RealVector &solution){
        
    int num_primary_eqs = opts.get<int>("num-primary-equations");
    if ( num_primary_eqs <= 0 ){
      std::string msg;
      msg = "EqConstrainedLSQSolver::solve() set num-primary-equations";
      throw(std::runtime_error(msg));
    }

    if ( (num_primary_eqs>A.numCols()) ){
      std::string msg = "EqConstrainedLSQSolver::solve() ";
      msg += "num-primary-equations is not smaller than the number of columns ";
      msg += "in A";
      throw(std::runtime_error(msg));
    }

    if ( num_primary_eqs > A.numRows() ){
      std::string msg = "EqConstrainedLSQSolver::solve() ";
      msg += "num-primary-equations is larger than the number of rows in A";
      throw(std::runtime_error(msg));
    }

    if (A.numRows() < A.numCols()){
      std::string msg = "EqConstrainedLSQSolver::solve() A is ";
	msg += "underdetermined";
      throw(std::runtime_error(msg));
    }

    RealMatrix solutions, metrics;
	
    RealMatrix C_eq( Teuchos::View, A, num_primary_eqs, A.numCols(), 0, 0 );
    RealMatrix A_eq( Teuchos::View, A, A.numRows() - num_primary_eqs,
		     A.numCols(), num_primary_eqs, 0 );
    RealVector d_eq( Teuchos::View, b.values(), num_primary_eqs );
    RealVector b_eq( Teuchos::View, b.values() + num_primary_eqs,
		     b.numRows() - num_primary_eqs );
    equality_constrained_least_squares_solve( A_eq, b_eq, C_eq, d_eq,
					      solution );

    RealVector residual( b_eq );
    residual.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS,
    		       -1.0, A_eq, solutions, 1.0 );
    return residual.normFrobenius();

  };
};


}  // namespace util
}  // namespace Pecos

#endif  // include guard
