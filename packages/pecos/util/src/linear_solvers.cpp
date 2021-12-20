/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "linear_solvers.hpp"

namespace Pecos {
namespace util {

std::shared_ptr<LinearSystemSolver> regression_solver_factory(OptionsList &opts){
  std::string name = "regression_type";
  RegressionType regression_type =
    get_enum_enforce_existance<RegressionType>(opts, name);
  
  bool use_cross_validation = opts.get("use-cross-validation", false);
  if (use_cross_validation){
    std::shared_ptr<CrossValidatedSolver> cv_solver(new CrossValidatedSolver);
    cv_solver->set_linear_system_solver(regression_type);
    return cv_solver;
  }
  
  switch (regression_type){
  case ORTHOG_MATCH_PURSUIT : {
    std::shared_ptr<OMPSolver> omp_solver(new OMPSolver);
      return omp_solver;
  }
  case LEAST_ANGLE_REGRESSION : {
    std::shared_ptr<LARSolver> lars_solver(new LARSolver);
    lars_solver->set_sub_solver(LEAST_ANGLE_REGRESSION);
    return lars_solver;
  }
  case LASSO_REGRESSION : {
    std::shared_ptr<LARSolver> lars_solver(new LARSolver);
    lars_solver->set_sub_solver(LASSO_REGRESSION);
    return lars_solver;
  }
  case EQ_CONS_LEAST_SQ_REGRESSION : {
    std::shared_ptr<EqConstrainedLSQSolver>
      eqlsq_solver(new EqConstrainedLSQSolver);
    return eqlsq_solver;
  }
  case SVD_LEAST_SQ_REGRESSION: case LU_LEAST_SQ_REGRESSION:
  case QR_LEAST_SQ_REGRESSION:
  {
    //\todo add set_lsq_solver to lsqsolver class so we can switch
    // between svd, qr and lu factorization methods
    std::shared_ptr<LSQSolver> lsq_solver(new LSQSolver);
    return lsq_solver;
  }
  default: {
    throw(std::runtime_error("Incorrect \"regression-type\""));
  }
  }
}

std::shared_ptr<OMPSolver> cast_linear_system_solver_to_ompsolver(std::shared_ptr<LinearSystemSolver> &solver){
  std::shared_ptr<OMPSolver> solver_cast =
    std::dynamic_pointer_cast<OMPSolver>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to OMPSolver shared_ptr"));
  return solver_cast;
}

std::shared_ptr<LARSolver> cast_linear_system_solver_to_larssolver(std::shared_ptr<LinearSystemSolver> &solver){
  std::shared_ptr<LARSolver> solver_cast =
    std::dynamic_pointer_cast<LARSolver>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to LARSolver shared_ptr"));
  return solver_cast;
}

std::shared_ptr<LSQSolver> cast_linear_system_solver_to_lsqsolver(std::shared_ptr<LinearSystemSolver> &solver){
  std::shared_ptr<LSQSolver> solver_cast =
    std::dynamic_pointer_cast<LSQSolver>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to LSQSolver shared_ptr"));
  return solver_cast;
}

std::shared_ptr<EqConstrainedLSQSolver> cast_linear_system_solver_to_eqconstrainedlsqsolver(std::shared_ptr<LinearSystemSolver> &solver){
  std::shared_ptr<EqConstrainedLSQSolver> solver_cast =
    std::dynamic_pointer_cast<EqConstrainedLSQSolver>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to EqConstrainedLSQSolver shared_ptr"));
  return solver_cast;
}

}  // namespace util
}  // namespace Pecos
