/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "CrossValidatedSolver.hpp"
#include "linear_solvers.hpp"
#include "LSQCrossValidationIterator.hpp"

namespace Pecos {
namespace util {

void CrossValidatedSolver::unnormalize_coefficients(const RealVector &column_norms){
 throw(std::runtime_error("We never need to do this because this is handled internally by call to the linear solver we are wrapping."));
}

CrossValidatedSolver::CrossValidatedSolver() : LinearSystemSolver() {}
  
CrossValidatedSolver::~CrossValidatedSolver(){}

void CrossValidatedSolver::
set_linear_system_solver(RegressionType regression_type){
  OptionsList opts;
  opts.set("regression_type", regression_type);
  cvIterator_ = linear_system_cross_validation_iterator_factory(opts);
};

void CrossValidatedSolver::
solve(const RealMatrix &A, const RealMatrix &B, OptionsList & opts ){
  RealMatrix A_copy(Teuchos::Copy, A, A.numRows(), A.numCols());
  RealMatrix B_copy(Teuchos::Copy, B, B.numRows(), B.numCols());
  multi_rhs_solve(A_copy, B_copy, opts);
}

void CrossValidatedSolver::
multi_rhs_solve(const RealMatrix &A, const RealMatrix &B,
		OptionsList & opts){
  if (!cvIterator_)
    throw(std::runtime_error("Must call set_linear_system_solver"));
  
  OptionsList cv_opts;
  if (opts.isType<OptionsList>("cv-opts"))
    cv_opts = opts.get<OptionsList>("cv-opts");
  else
    throw(std::runtime_error("Must specify cv-opts"));

  // We want to allow CrossValidatedSolver to be called like other solvers
  // where opts just contains all regression_opts. We also want to
  // respect that CVIterators take opts consiting of cv_opts and sublist
  // consisting of regression_opts. So do internal manipulation of opts
  // to achieve this goal here.
  OptionsList regression_opts(opts);
  regression_opts.erase("cv-opts");
  OptionsList cv_opts_copy(cv_opts);
  cv_opts_copy.set("regression-opts", regression_opts);

  cvIterator_->run(A, B, cv_opts_copy);
  cvIterator_->generate_best_solutions(A, B, solutions_, residualNorms_,
                                       regression_opts);
}

void CrossValidatedSolver::get_best_scores(RealVector &result) const{
  cvIterator_->get_best_scores(result);
}

void CrossValidatedSolver::
get_solutions_for_all_regularization_params(RealMatrix &result_0,
					    int rhs_num) const{
  IntVector column_indices(1,false); column_indices[0]=0;
  extract_submatrix_from_column_indices(solutions_, column_indices, result_0);
}

void CrossValidatedSolver::get_residuals_for_all_regularization_params(
      RealVector &result_0, int rhs_num) const{
  throw(std::runtime_error("This function does not make sense for all solvers. consider refining class hierarchy"));
}

void CrossValidatedSolver::get_final_solutions(RealMatrix &result_0) const{
  result_0 = solutions_;
}

void CrossValidatedSolver::get_final_residuals(RealVector &result_0) const{
  result_0 = residualNorms_;
}

std::shared_ptr<LinearSystemCrossValidationIteratorBase> CrossValidatedSolver::
get_cross_validation_iterator(){
  return cvIterator_;
}
  
std::shared_ptr<LinearSystemCrossValidationIteratorBase> linear_system_cross_validation_iterator_factory(OptionsList &opts){
  std::string name = "regression_type";
  RegressionType regression_type =
    get_enum_enforce_existance<RegressionType>(opts, name);
  
  switch (regression_type){
  case ORTHOG_MATCH_PURSUIT : case LEAST_ANGLE_REGRESSION :
  case LASSO_REGRESSION : case EQ_CONS_LEAST_SQ_REGRESSION : {
    std::shared_ptr<LinearSystemSolver> solver =
      regression_solver_factory(opts);
    std::shared_ptr<LinearSystemCrossValidationIterator>
      cv_iterator(new LinearSystemCrossValidationIterator);
    cv_iterator->set_linear_system_solver(solver);
    return cv_iterator;
    break;
  }
  /// Use fast leasts squares cross validation using cholesky factorization
  /// of normal equations, regardless of type of factorization specified.
  case SVD_LEAST_SQ_REGRESSION : case LU_LEAST_SQ_REGRESSION :
  case QR_LEAST_SQ_REGRESSION :{
    std::shared_ptr<LSQCrossValidationIterator>
      cv_iterator(new LSQCrossValidationIterator);
    return cv_iterator;
    break;
  }
  default: {
    throw(std::runtime_error("Incorrect \"regression-type\""));
  }
  }
}

std::shared_ptr<CrossValidatedSolver> cast_to_cross_validated_solver(std::shared_ptr<LinearSystemSolver> &solver){
  std::shared_ptr<CrossValidatedSolver> solver_cast =
    std::dynamic_pointer_cast<CrossValidatedSolver>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to CrossValidatedSolver shared_ptr"));
  return solver_cast;
}


}  // namespace util
}  // namespace Pecos
