/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_LSQ_SOLVER_HPP
#define PECOS_UTIL_LSQ_SOLVER_HPP

#include "LinearSystemSolver.hpp"

namespace Pecos {
namespace util {

//----------------------------------------------------------------
/**
 *\class LSQSolver
 *\brief Solve \f$\min\lVert Ax-b\rVert_2\f$ using either SVD, QR, or LU 
 * factorization.
 */
class LSQSolver : public LinearSystemSolver {
protected:

  /// The solutions X to AX=B. X is a matrix (A.numCols() x B.numCols())
  RealMatrix solutions_;
  /// The l2 norm of the residuals AX-B. residuals_ is a vector (B.numCols() x 1)
  RealVector residuals_;

  /**\copydoc LinearSystemSolver::unnormalize_coefficients()*/
  void unnormalize_coefficients(const RealVector &column_norms){
    adjust_coefficients(column_norms, solutions_);
  };
  
public:

  /// Default constructor
  LSQSolver() : LinearSystemSolver() {};

  /// Destructor
  ~LSQSolver(){};

  /**\copydoc LinearSystemSolver::get_solutions_for_all_regularization_params()*/
  void get_solutions_for_all_regularization_params(
      RealMatrix &result_0, int rhs_num) const{
    shape_uninitialized(result_0, solutions_.numRows(), 1);
    for (int i=0; i<solutions_.numRows(); ++i)
      result_0(i,0)=solutions_(i, rhs_num);
  };

  /**\copydoc LinearSystemSolver::get_residuals_for_all_regularization_params()*/ 
  void get_residuals_for_all_regularization_params(
      RealVector &result_0, int rhs_num) const{
    size_uninitialized(result_0, 1);
    result_0[0]=residuals_[rhs_num];
  };

  /**\copydoc LinearSystemSolver::get_final_solutions()*/  
  void get_final_solutions(RealMatrix &result_0) const{
    result_0 = solutions_;
  };

  /**\copydoc LinearSystemSolver::get_final_residuals()*/  
  void get_final_residuals(RealVector &result_0) const{
    result_0 = residuals_;
  };

  /**
   * \brief Solve \f$\min\lVert Ax-b\rVert_2\f$ using either SVD, QR, or LU 
   * factorization.
   *
   * \param[in] A matrix (m x n)
   *    The A matrix. 
   * \param[in] b vector (m x num-rhs)
   *    The rhs. 
   * \param[in] opts  
   *     List of options
   *
   * \param[out] solution vector (n x num-rhs)
   *    The least-sqaures solutions for each rhs
   * 
   * opts (optional parameters)
   * -------------------------
   * 
   * "regression_type" : RegressionType default=SVD_LEAST_SQ_REGRESSION
   *     The factorization method used to solve the least squares problem
   *     Accepted values are SVD_LEAST_SQ_REGRESSION, QR_LEAST_SQ_REGRESSION,
   *     LU_LEAST_SQ_REGRESSION
   */
  void multi_rhs_solve(const RealMatrix &A, const RealMatrix &B,
                       OptionsList& opts ){
    std::string name = "regression_type";
    RegressionType regression_type =
      get_enum_apply_default(opts, name, SVD_LEAST_SQ_REGRESSION);

    if ( A.numRows() < A.numCols() )
      std::cout << "LSQSolver::solve() Warning A is under-determined. " <<
	"M = " << A.numRows() << " N = " << A.numCols() <<
	". Returning minimum norm solution\n";

    switch (regression_type){
    case SVD_LEAST_SQ_REGRESSION:{
      int rank(0);
      RealVector singular_values;
      Real solver_tol = opts.get("rcond_tol", 1.e-6);
      svd_solve( A, B, solutions_, singular_values, rank, solver_tol );
      break;
    }
    case QR_LEAST_SQ_REGRESSION:{
      Teuchos::ETransp trans = opts.get("transpose", Teuchos::NO_TRANS);
      qr_solve( A, B, solutions_, trans );
      break;
    }
    case LU_LEAST_SQ_REGRESSION:{
      Teuchos::ETransp trans = opts.get("transpose", Teuchos::NO_TRANS);
      lu_solve(const_cast<RealMatrix&>(A), B, solutions_, true, trans );
      break;
    }
    default: {
      throw(std::runtime_error("Incorrect regression type"));
    }
    }
    RealMatrix residual( B );
    residual.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS,
		       -1.0, A, solutions_, 1.0 );
    size_uninitialized(residuals_, B.numCols());
    for (int j=0; j<B.numCols(); ++j)
      residuals_[j] = residual.normFrobenius();
      
  };
};


}  // namespace util
}  // namespace Pecos

#endif  // include guard
