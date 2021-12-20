/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_LINEAR_SYSTEM_SOLVER_HPP
#define PECOS_UTIL_LINEAR_SYSTEM_SOLVER_HPP

#include "teuchos_data_types.hpp"
#include "OptionsList.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "linear_algebra.hpp"

#include <boost/numeric/conversion/cast.hpp>

namespace Pecos {
namespace util {

  enum RegressionType{
    SVD_LEAST_SQ_REGRESSION, EQ_CONS_LEAST_SQ_REGRESSION,
    ORTHOG_MATCH_PURSUIT, LASSO_REGRESSION, LEAST_ANGLE_REGRESSION,
    BASIS_PURSUIT, BASIS_PURSUIT_DENOISING, QR_LEAST_SQ_REGRESSION,
    LU_LEAST_SQ_REGRESSION
  };
  

  /**
   *\class LinearSystemSolver
   *\brief The base class of all wrapper of methods for solving regression type problems, e.g least squares, compressive sensing, etc.
   */
  class LinearSystemSolver{
  protected:
  
    /**\brief adjust the stored solutions coefficients to account for 
     * normalization of the columns of A when solving Ax=B
     * \param[in] column_norms vector (A.numRows() x 1)
     *    The column norms of A
     */
    virtual void unnormalize_coefficients(const RealVector &column_norms)=0;
  
  public:

    /// Default constructor
    LinearSystemSolver() { };

    /// Destructor
    virtual ~LinearSystemSolver() { };

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
     */
    virtual void solve(const RealMatrix &A, const RealMatrix &B,
                       OptionsList & opts){
      bool normalize_inputs = opts.get("Normalize Inputs", false);
      bool precondition = opts.get("precondition", false);

      if ((A.numCols()==0)||(A.numRows()==0)||(B.numCols()==0)||(B.numRows()==0))
        throw(std::runtime_error("A and/or B had zero columns and/or rows"));

      RealMatrix A_copy(Teuchos::Copy, A, A.numRows(), A.numCols());
      RealMatrix B_copy(Teuchos::Copy, B, B.numRows(), B.numCols());
      RealVector column_norms;

      if (precondition)
        precondition_linear_system(A_copy, B_copy);
    
      if ( normalize_inputs )
        normalise_columns(A_copy, column_norms);

      multi_rhs_solve(A_copy, B_copy, opts);

      if (normalize_inputs)
        unnormalize_coefficients(column_norms);
    };
  
    /**
     * \brief Find solutions to \f$AX \approx B\f$ 
     * The solutions are stored internally and not returned. 
     * Most users should not call this function and use solve() instead.
     * This function implements efficient multiple rhs solves and/or
     * extraction and application of options in opts specific to derived class.
     *
     * \param[in] A matrix (num_rows x num_cols)
     *     The matrix A
     * \param[in] B vector (num_cols x num_rhs)
     *     The rhs B
     * \param[in] opts  
     *     List of options
     */
    virtual void multi_rhs_solve( const RealMatrix &A, const RealMatrix &B,
                                  OptionsList & opts )=0;


    /**\brief Normalize the columns of A to have unit l2 norm
     * \param[in] A matrix (num_rows x num_cols)
     * \param[out] result vector (num_cols x 1)
     *    The l2 norms of the columns of A
     */
    void normalise_columns( RealMatrix &A, RealVector &result ){
      int M = A.numRows(), N = A.numCols();
      result.sizeUninitialized( N );
      for (int i=0; i<N; i++){
        RealVector col( Teuchos::View, A[i], M );
        result[i] = col.normFrobenius();
        col *= 1./result[i];
      }
    };

    /**\brief If the columns of A have been normalized
     * before solving Ax=B adjust coefficients to account for this normalization
     * \param[in] result vector (num_cols x 1)
     *    The factor used to normalize A such that new A = A/factors
     * \param[out] result matrix (num_cols x num_solutions)
     *    The adjusted solutions
     */
    void adjust_coefficients( const RealVector &normalisation_factors, 
                              RealMatrix &result ) const{
      int num_coeff = result.numRows(), num_qoi = result.numCols();
      for (int i=0; i<num_qoi; i++){
        for ( int j = 0; j < num_coeff; j++ )
          result(j,i) /= normalisation_factors[j];
      } 
    };

    /**\brief Apply preconditioning to the system Ax=B
     * \param[inout] A matrix (num_rows x num_cols)
     *     On entry the matrix A on exit the preconditioned matrix W*A
     * \param[inout] B matrix (num_cols x num_rhs)
     *     On entry the rhs B on exit the preconditioned matrix W*B
     */
    void precondition_linear_system(RealMatrix &A, RealMatrix &B) const{
      // call precondition Object
      //preconditioner_.apply(A,B);
    };

    //void set_preconditioner(LinearSystemPreconditioner & preconditioner){
    // preconditioner_ = preconditioner
    //};
  
    /**\brief Get the solutions x to Ax=B corresponding to all values of
       the regularization parameter(s) considered
       * \param[in] rhs_num
       *     The column of B that we want the solutions for
       * \param[out] result_0 (num_coeff x num_reg_params) matrix
       *     The solutions for each of the regularization parameters
       */
    virtual void get_solutions_for_all_regularization_params(
                                                             RealMatrix &result_0, int rhs_num) const=0;

    /**\brief Get the l2 norm of the residuals Ax-B corresponding to all values of
       the regularization parameter(s) considered
       * \param[in] rhs_num
       *     The column of B that we want the solutions for
       * \param[out] result_0 (num_reg_paramsx1) vector
       *     The l2 norm of the residual for each of the regularization parameters
       */
    virtual void get_residuals_for_all_regularization_params(
                                                             RealVector &result_0, int rhs_num) const=0;

    /**\brief Get the solutions to Ax-B corresponding to the final values of
       the regularization parameter(s) for each rhs
       * \param[out] result (num_coeff x num_rhs) vector
       *     The final solutions for each RHS
       */
    virtual void get_final_solutions(RealMatrix &result_0) const=0;

    /**\brief Get the l2 norm of the residuals to Ax-B 
     * corresponding to the final values of the regularization parameter(s) 
     * for each rhs
     * \param[out] result (num_coeff x num_rhs) vector
     *     The final l2 norm of the residuals for each RHS
     */
    virtual void get_final_residuals(RealVector &result_0) const=0;
  };

  /**
   *\class SparseSolver
   *\brief A wrapper for all compressive sensing solvers.
   */
  class SparseSolver : public LinearSystemSolver{
  protected:
    /// Solutions of Ax=b num_rhs x num_coeff x num_regularization parameters
    std::vector<RealMatrix> solutions_;
  
    /// Residuals (Ax-b) num_rhs x num_regularization parameters
    std::vector<RealVector> residuals_;

    /// The weights used to solve weighted l1 minimization problem.
    /// matrix of size (num_coeff x num_rhs)
    RealMatrix sparsityWeights_;

    /// The residual tolerances t used as a stopping criteria when solving
    /// min \f$\lVert x\rVert_1\f$ subject to \f$\lVert Ax-b\rVert_2<t\f$
    /// matrix of size (num_rhs x num_tolerances)
    RealMatrix residualTolerances_;

    /**\copydoc LinearSystemSolver::unnormalize_coefficients()*/
    void unnormalize_coefficients(const RealVector &column_norms){
      for (size_t i=0; i<solutions_.size(); ++i)
        adjust_coefficients(column_norms, solutions_[i]);
    };
  
  public:

    /// Default constructor
    SparseSolver() : LinearSystemSolver() {};

    /// Destructor
    virtual ~SparseSolver() {};

    /**\brief Set the weights used to solve weighted l1 minimization problem.
     * \param[in] weights matrix (num_coeff x num_rhs)
     *    The weights used to solve weighted l1 minimization problem.
     */
    void set_sparsity_weights(const RealMatrix &weights){
      sparsityWeights_ = weights;
    };

    /*\brief Set the residual tolerances t used as a stopping criteria when 
     * solving min \f$\lVert x\rVert_1\f$ subject to \f$\lVert Ax-b\rVert_2<t\f$
     * \param[in] tolerances matrix of size (num_rhs x num_tolerances)
     *     The tolerances on the l2 residual norm
     */
    void set_residual_tolerances(const RealMatrix &tolerances){
      residualTolerances_ = tolerances;
    }

    /**\copydoc LinearSystemSolver::get_solutions_for_all_regularization_params()*/  
    void get_solutions_for_all_regularization_params(
                                                     RealMatrix &result_0, int rhs_num) const{
      result_0=solutions_[rhs_num];
    };

    /**\copydoc LinearSystemSolver::get_residuals_for_all_regularization_params()*/  
    void get_residuals_for_all_regularization_params(
                                                     RealVector &result_0, int rhs_num) const{
      result_0=residuals_[rhs_num];
    };

    /**
     * \brief Find solutions to \f$AX \approx b\f$ where \f$b\f$ 
     * for a set of regularization parameters and only one RHS. 
     * The regularization parameters are either contained in opts or determined 
     * internally by the algorithm implemented.
     *
     * \param[in] A matrix (num_rows x num_cols)
     *     The matrix A
     * \param[in] b vector (num_cols x 1)
     *     The rhs b
     * \param[in] opts  
     *     List of options
     *
     * \param[out] result_0 matrix (num_cols x num_reg_params)
     *     The solutions orresponding to all values ofthe regularization 
     *     parameter(s) considered
     * \param[out] result_1 maxtrix (2 x num_reg_params)  
     *     First row is l2 norm of residual for each regularization parameter
     *     Second row is number of non-zeros
     */
    virtual void single_rhs_solve( const RealMatrix &A, const RealVector &b,
                                   OptionsList & opts,
                                   RealMatrix &result_0, RealVector &result_1)=0;
  
    /**
     * \brief Find (possibly weighted) sparse solutions to \f$AX \approx B\f$ 
     * for a set of regularization parameters and possibly many RHS. 
     * The regularization parameters are either contained in opts or 
     * determined internally by the algorithm implemented. The solutions are
     * stored internally and not returned.
     * Most users should not use this function and use solve() instead.
     *
     * \param[in] A matrix (num_rows x num_cols)
     *     The matrix A. Not const because may need to modify it
     * \param[in] B vector (num_cols x num_rhs)
     *     The rhs B
     * \param[in] opts  
     *     List of options
     */
    virtual void multi_rhs_solve( const RealMatrix &A, const RealMatrix &B,
                                  OptionsList & opts){
      // Set weights.
      // Allow for situation where weights were set already with
      // set_sparsity_weights. if they exist in opts weights_ will be overwritten
      bool use_weights = false;
      if (opts.isType<RealMatrix>("sparsity-weights")){
        RealMatrix weights = opts.get<RealMatrix>("sparsity-weights");
        set_sparsity_weights(weights);
        if ((sparsityWeights_.numRows()>0)&&(sparsityWeights_.numCols()>0)){
          if ((sparsityWeights_.numRows()!=A.numCols())||
              ((sparsityWeights_.numCols()!=B.numCols()))){
            std::string msg = "shape of weights is inconsistent with A and B";
            throw(std::runtime_error(msg));
          }
          use_weights=true;
        }
      }

      bool use_residual_tolerances = false;
      if (opts.isType<RealMatrix>("residual-tolerances")){
        RealMatrix residual_tolerances =
          opts.get<RealMatrix>("residual-tolerances");
        set_residual_tolerances(residual_tolerances);
        if ((residualTolerances_.numRows()>0)&&(residualTolerances_.numCols()>0)){
          if (residualTolerances_.numRows()!=B.numCols()){
            std::string msg = "shape of weights is inconsistent with B.numCols()";
            throw(std::runtime_error(msg));
          }
          use_residual_tolerances = true;
        }
      }
    
      RealVector b, // single rhs
        weights_i,  // weights for solving weighted sparse minimization
        tols_i;     // tolerances to stop minimization
      if( solutions_.size()!=(size_t)B.numCols() )
        solutions_.resize(B.numCols());
      if( residuals_.size()!=(size_t)B.numCols() )
        residuals_.resize(B.numCols());
      for (int i=0; i<B.numCols(); ++i){
        if (use_weights){
          weights_i = Teuchos::getCol(Teuchos::View, sparsityWeights_, i);
          opts.set("weights", weights_i);
        }
        if (use_residual_tolerances){
          size_uninitialized(tols_i, residualTolerances_.numCols());
          for (int j=0; j<tols_i.length(); ++j)
            tols_i[j] = residualTolerances_(i,j);
          opts.set("residual-tols-single-rhs", tols_i);
        }
        b = Teuchos::getCol(Teuchos::View, const_cast<RealMatrix &>(B), i);
        single_rhs_solve(A, b, opts, solutions_[i], residuals_[i]);
      }
    };

    void apply_weights_to_matrix(const RealVector &weights, RealMatrix &A) const{
      int M = A.numRows(), N = A.numCols();
      for (int i=0; i<N; i++){
        RealVector col( Teuchos::View, A[i], M );
        col *= weights[i];
      }
    }
  
    /**\copydoc  LinearSystemSolver::get_final_solutions()*/
    virtual void get_final_solutions(RealMatrix &result_0) const{
      shape_uninitialized(result_0,solutions_[0].numRows(),boost::numeric_cast<int>(solutions_.size()));
      for (size_t j=0; j<solutions_.size(); ++j){
        int num_solutions = solutions_[j].numCols();
        for (size_t i=0; i<(size_t)solutions_[j].numRows(); ++i)
          result_0(i,j) = solutions_[j](i,num_solutions-1);
      }
    };

    /**\copydoc  LinearSystemSolver::get_final_residuals()*/
    virtual void get_final_residuals(RealVector &result_0) const{
      size_t num_rhs = residuals_.size();
      size_uninitialized(result_0, boost::numeric_cast<int>(num_rhs));
      for (size_t i=0; i<num_rhs; ++i){
        int num_solutions = residuals_[i].length();
        result_0[i] = residuals_[i][num_solutions-1];
      }
    };

    void get_residual_tolerances(RealMatrix &result_0) const{
      result_0 = residualTolerances_;
    };
  };
  


  /*
  **
  *\class WeightedSparseSolver
  *\brief A wrapper for all iterative reweighting for sparse solvers.
  *
  class WeightedSparseSolver : public SparseSolver{
  public:
  
  WeightedSparseSolver() : SparseSolver() {};

  virtual ~WeightedSparseSolver() {};
  
  
  void cross_validated_unweighted_solve(const RealMatrix &A, const RealVector &b, RealMatrix &result_0, RealMatrix &result_1, OptionsList & opts, OptionsList & out){
  if (opts.isParameter("cv_opts")){ 	
  cv_iterator = MultipleSolutionLinearModelCrossValidationIterator();
  cv_iterator.set_solve_function(unweighted_solve);
  if (opts.isParameter("fault_data")){
  FaultData fault_data = opts.get<FaultData>("fault_data");
  cv_iterator.set_fault_data(fault_data);
  }
  OptionsList cv_opts;
  if (opts.isParameter("cv_opts"))
  cv_opts=opts.get<OptionsList>("cv_opts");
  cv_iterator.run(cv_opts);
  RealVector residual_tols(1,false);
  residual_tols[0] = cv_iterator.get_best_residual_tolerance();
  opts.set("residual tols",residual_tols);
  }
  unweighted_solve(A, b, result_0, result_1, opts, out);
  };

  **\brief Solve problem using iterative reweighting with repeated calls to
  unweighted_solve(). 1 reweighting iters is the special case of single unweighetd solve.
  *
  void iterative_reweighting_solve(const RealMatrix &A, const RealVector &b, RealMatrix &result_0, RealMatrix &result_1, OptionsList & opts){
    
  int num_rows = Amatrix.numRows(), num_cols = Amatrix.numCols();
  int max_reweighting_iters = opts.get("max_reweighting_iters", 1)
    
  int it = 0;
  RealVector weights(num_cols,1.);
  RealMatrix weighted_A(A);
  Real prev_cv_error = std::numeric_limits<float>::max();
  RealMatrix coeff, best_coeff;
  while (it<max_reweighting_iters){
  for (int j=0; j<num_cols; ++j){
  for (int i=0; i<num_rows; ++i)
  weighted_A(i,j) /= weights[j];
  }
  OptionsList out;
  unweighted_solve_using_single_rhs(A, b, coeff, metrics, opts, out);
  for (int j=0; j<num_cols; ++j){
  coeff[j] /= weights[j];
  weights[j] = 1./(reweighting_delta+std::abs(coeff[j]));
  }
  if (out.isType<Real>("CV error")){
  cv_error = get<Real>("CV error");
  if (cv_error >= prev_cv_error){
  coeff = best_coeff;
  if (verbosity_>1){
  std::cout << 
  "Reweighting degraded the approximation stopping at " <<
  "iter " << it << " of " << max_reweighting_iters << std::endl;
  }
  break;
  }else{
  best_coeff = coeff;
  prev_cv_error = cv_error;
  }
  }
  }
  };
  };*/

}  // namespace util
}  // namespace Pecos

#endif  // include guard
