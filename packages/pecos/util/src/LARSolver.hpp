/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_LAR_SOLVER_HPP
#define PECOS_UTIL_LAR_SOLVER_HPP

#include "LinearSystemSolver.hpp"
#include "least_angle_regression.hpp"

namespace Pecos {
namespace util {

/**
 *\class LARSolver
 *Find the solution to 
 *   \f$\min \lVert x\rVert_0 \mathrm{s.t} \lVert Ax - b\rVert_2 < \epsilon\f$
 * using least angle regression 
 *   \f$\min \lVert x\rVert_1 \mathrm{s.t} \lVert Ax - b\rVert_2 < \epsilon\f$
 * using least angle regression with the LASSO modification. This optimization
 * problem as written is not the lasso problem which is
 *   \f$\min \lVert Ax - b\rVert_2  \mathrm{s.t} \lVert x\rVert_1 < delta\f$
 * however these problems are equivalent if \f$\delta\f$ and \f$\epsilon\f$
 * are chosen correctly. We take advantage of this to provide a stopping 
 * criterion on $\espilon$ which is more intuitive than using $\delta$.
 */
class LARSolver : public SparseSolver{
private:

  /// Flag specifying whether to use LAR with or without LASSO modification
  int solver_;
  
  Real delta_;
  
public:
  
  /// Default constructor
  LARSolver() : SparseSolver(), solver_( LEAST_ANGLE_REGRESSION ),
		 delta_( 0.0 ){};

  /// Destructor
  ~LARSolver(){clear();};

  /// Reset the solver
  void clear(){
    solver_ = LEAST_ANGLE_REGRESSION; delta_ = 0.0;
  };

  /// Specify whether to solve least angle regression with or without the
  /// lasso modification. TODO move to an opts passed to solve()
  // This will allow simplification of sparsesolver multi_rhs_solve
  // so that weights can be treated correctly, e.g set weights or scale A with
  // weights
  void set_sub_solver( int solver_id ){
    if ( ( solver_id != LASSO_REGRESSION ) &&
	 ( solver_id != LEAST_ANGLE_REGRESSION ) ){
      std::stringstream msg;
      msg << "set_sub_solver() solver id must be either: " << LASSO_REGRESSION
	  << " or " << LEAST_ANGLE_REGRESSION << "\n";
      throw( std::runtime_error( msg.str() ) );
    }
    solver_ = solver_id;
  };

  /// Deprecated do not use
  void set_delta( Real delta ){
    delta_ = delta;
  }
  
  /**
   * \brief Find the solution to 
   *   \f$\min \lVert x\rVert_0 \mathrm{s.t} \lVert Ax - b\rVert_2 < \epsilon\f$
   * using least angle regression 
   *   \f$\min \lVert x\rVert_1 \mathrm{s.t} \lVert Ax - b\rVert_2 < \epsilon\f$
   * using least angle regression with the LASSO modification. This optimization
   * problem as written is not the lasso problem which is
   *   \f$\min \lVert Ax - b\rVert_2  \mathrm{s.t} \lVert x\rVert_1 < delta\f$
   * however these problems are equivalent if \f$\delta\f$ and \f$\epsilon\f$
   * are chosen correctly. We take advantage of this to provide a stopping 
   * criterion on $\espilon$ which is more intuitive than using $\delta$.
   *
   * \param[in] A matrix (m x n)
   *    The A matrix. 
   * \param[in] b vector (m x 1)
   *    The rhs. 
   * \param[in] opts  
   *     List of options
   *
   * \param[out] solution vectors (n x 1) or (n x num-iters)
   *    The LAR/LASSO solutions. If store-history is true then
   *    solutions  will be (n x num-iters), where num-iters is the number
   *    of iterations taken by the algorithm, otherwise solutions will 
   *    be (n x 1) and contain only the solution from the final step of the 
   * algorithm
   *
   * opts (optional parameters)
   * -------------------------
   * "verbosity" : integer in [0,infinity) default=0
   *     controls amount of print statements
   *
   * "normalize-choice" : boolean  default=true
   *     If true base choice of next column to add by finding
   *     maximum correlation of correlation between column and residual 
   *     normalized by l2 norm of each column
   *
   * "max-iters" : integer default=A.numCols()
   *     The maximum number of iterations of the LAR algorithm. If 
   *     solver-type=LARS_REGRESSION then final solution will have 
   *     a maximum number of non-zeros=max-iters. However if 
   *     solver-type=LASSO_REGRESSION the maximum number of non-zeros in the 
   *     final solution will be bounded by max-iters.
   *
   * "residual-tolerance" : integer default=A.numCols()
   *    The minimum allowable value of \f$\epsilon\f$. When the tolerance
   *    is exceeded the algorithm will terminate. The final solution obtained
   *    will be the fist solution which produced a resdiual less than or equal
   *    to the tolerance.  This value will be overwritten by the value in the
   *    vector residual-tols-single-rhs if it is specified.
   *
   * "store-history" : boolean default=true
   *     If true store the solution at each iteration of the OMP algorithm.
   *     If false only store the final solution
   *
   * "non-negative" : boolean default=true
   *     If true enforce the non-zero elements of the solution are all positive.
   *     If false solve the standard problem.
   * 
   * "weights" : Teuchos::SerialDenseVector (A.numCols x 1) default=empty
   *     Non-negative weights used to solve weighted minimization problem
   *     \f$\min \lVert x\rVert_0 such that \lVert AWx - b\rVert_2 < eps\f$
   *     where W is a diagonal matrix whose non-zero entries are given by
   *     the vector. The default is empty which corresponds to all weights=1.
   *
   * "residual-tols-single-rhs" : Teuchos::SerialDenseVector (1 x 1) default=empty
   *    The minimum allowable value of \f$\epsilon\f$. When the tolerance
   *    is exceeded the algorithm will terminate. The final solution obtained
   *    will be the fist solution which produced a resdiual less than or equal
   *    to the tolerance. If specified this value will overwrite the one given 
   *    by residual-tolerance.
   *
   * "memory-chunk-size" : integer default = min(500,"max-iters")
   *    Algorithm internally allocates memory for solutions in 
   *    chunks to avoid allocating unnecessary memory for sparsity << A.numCols()
   *    This is the size of those chunks
   */
  void single_rhs_solve( const RealMatrix &A, const RealVector &b,
			 OptionsList& opts,
			 RealMatrix &result_0, RealVector &result_1){
    int verbosity         = opts.get("verbosity", 0);
    bool normalize_choice = opts.get("normalize-choice", false);
    int max_iters         = opts.get("max-iters", 10*A.numRows());
    Real residual_tol     = opts.get("residual-tolerance", 0.);
    bool store_history     = opts.get("store-history",true);
    bool non_negative     = opts.get("non-negative",false);
    int maxNNZ            =
      opts.get("max-num-non-zeros",std::min(A.numRows(),A.numCols()));
    int memory_chunk_size = opts.get("memory-chunk-size",
				     std::min(500, max_iters));

       // todo consider creating stepbasedsparsesolver and share the common code
    // here and in larssolver
    if (opts.isType<RealVector>("residual-tols-single-rhs")){
      RealVector residual_tols=
        opts.get<RealVector>("residual-tols-single-rhs");
      if (residual_tols.length()!=1)
        throw(std::runtime_error("residual_tols vector must only have one entry"));
      residual_tol = (residual_tols)[0];
    }
    
    bool use_weights=false;
    RealVector weights;
    if (opts.isType<RealVector>("weights")){
      weights  = opts.get<RealVector>("weights");
      use_weights = (weights.length()==A.numCols());
    }
    RealMatrix A_(Teuchos::View, A.values(), A.stride(), A.numRows(),A.numCols());
    if (use_weights){
      A_.assign(A);
      apply_weights_to_matrix(weights, A_);
    }

    RealMatrix metrics;
    least_angle_regression( A, b, result_0, metrics,
			    residual_tol, solver_, delta_,
			    max_iters, maxNNZ,
			    verbosity, normalize_choice,
			    non_negative, store_history, memory_chunk_size);

    if (use_weights) adjust_coefficients(weights, result_0); 

    size_uninitialized(result_1,metrics.numCols());
    for (int i=0; i<result_1.length(); ++i)
      result_1[i] = metrics(0,i);
  };
};

}  // namespace util
}  // namespace Pecos

#endif  // include guard
