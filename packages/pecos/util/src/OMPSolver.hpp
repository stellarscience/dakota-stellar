/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_OMP_SOLVER_HPP
#define PECOS_UTIL_OMP_SOLVER_HPP

#include "LinearSystemSolver.hpp"
#include "orthogonal_matching_pursuit.hpp"

namespace Pecos {
namespace util {
 
/**
 *\class OMPSolver
 *\brief Greedily solve \f$\text{argmin} \lVert x\rVert_0\quad \text{s.t} \lVert Ax-b\rVert_2<\epsilon\f$ using orthogonal matching pursuit.
 */
class OMPSolver : public SparseSolver{
protected:
  IntVector ordering_; // enforce a set of columns to be chosen first
  
public:
  
  /// Default constructor
  OMPSolver() : SparseSolver() {};

  /// Destructor
  ~OMPSolver(){};

  /**
   * \brief Find a greedy solution to 
   *     \f$\min \lVert x\rVert_0 such that \lVert Ax - b\rVert_2 < \epsilon\f$
   * using orthogonal matching pursuit
   *
   * \param[in] A matrix (m x n)
   *    The A matrix. 
   * \param[in] b vector (m x 1)
   *    The rhs. 
   * \param[in] opts  
   *     List of options
   *
   * \param[out] solution vectors (n x 1) or (n x num-iters)
   *    The OMP solutions. If store-history is true then
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
   *     The maximum number of iterations of the OMP algorithm. For OMP
   *     this also coincides with the maxium number of non-zeros.
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
  void single_rhs_solve(const RealMatrix &A, const RealVector &b,
			OptionsList& opts,
			RealMatrix &result_0, RealVector &result_1){
    int verbosity         = opts.get("verbosity", 0);
    bool normalize_choice = opts.get("normalize-choice", true);
    int max_iters         = opts.get("max-iters", A.numCols());
    Real residual_tol     = opts.get("residual-tolerance", 0.);
    bool store_history    = opts.get("store-history",true);
    int memory_chunk_size = opts.get("memory-chunk-size",
				     std::min(A.numRows(),A.numCols()));

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
    orthogonal_matching_pursuit<Real>( A_, b, result_0, metrics,
				       residual_tol, max_iters, verbosity,
				       ordering_, normalize_choice,
				       store_history, memory_chunk_size );
    
    if (use_weights) adjust_coefficients(weights, result_0); 

    size_uninitialized(result_1,metrics.numCols());
    for (int i=0; i<result_1.length(); ++i)
      result_1[i] = metrics(0,i);
  };
    
  /**
   *\brief Set the ordering of the first set of columns of A that must be
   * selected before any other column can be selected.
   * \param[in] ordering vector (num_enforced_basis_terms)
   *    The columns of A which must be added before any other columns can be 
   *    added. Ensure num_enforced_basis_terms <= A.numCols().
   */
  void set_ordering( IntVector &ordering )
  {
    ordering_.resize( ordering.length() );
    ordering_.assign( ordering );
  };
};


}  // namespace util
}  // namespace Pecos

#endif  // include guard
