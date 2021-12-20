/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "LinearSystemCrossValidationIterator.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

namespace Pecos {
namespace util {

LinearSystemCrossValidationIteratorBase::
LinearSystemCrossValidationIteratorBase(){};

LinearSystemCrossValidationIteratorBase::
~LinearSystemCrossValidationIteratorBase(){};

void LinearSystemCrossValidationIteratorBase::
extract_matrix( const RealMatrix &A, const IntVector &indices, 
		RealMatrix &result ) const {
  if ( A.numRows() != numPts_ * numEquationsPerPoint_ )
    throw( std::runtime_error("extract_matrix: num pts and num equations per point are inconsistent with A") );

  int num_indices=indices.length();
  // We assume that indices will not contain any primary equation faults
  // as these indices are training or validation indices that have already
  // been parsed for primary equation failure
  int num_indices_with_secondary_eq = 0;
  for (int i=0; i<num_indices; i++){
    if (!failedRespData_[indices[i]])
      num_indices_with_secondary_eq++;
  }
  int num_submatrix_rows = indices.length() +
    (numEquationsPerPoint_-1)*num_indices_with_secondary_eq;
  shape_uninitialized(result, num_submatrix_rows, A.numCols());
  
  for (int j=0; j<A.numCols(); j++){
    int k=0;
    for (int i=0; i<num_indices; i++){
      // assumes there is a set of primary equations stored 
      // 0...b.length()-1
      // in b. Then the non-primary equations are stored 
      // 1...numEquationsPerPoint_ for first point then all non-primary 
      // equations for second point and so on. 
      result(i,j)=A(indices[i],j);
      // Determine if there are any faults associated with indices
      if (!failedRespData_[indices[i]]){
	int shift=numPts_ + indices[i] * (numEquationsPerPoint_-1);
	for (int m=0; m<numEquationsPerPoint_-1; m++){
	  result(num_indices+k,j)=A(shift+m,j);
	  k++;
	}
      }
    }
  }
}

void LinearSystemCrossValidationIteratorBase::
extract_linear_system( const RealMatrix &A, const RealVector &b,
		       const IntVector &indices, 
		       RealMatrix &result_0, RealVector &result_1 ) const {
  extract_matrix( A, indices, result_0 );
  extract_values( b, indices, result_1 );
}

void LinearSystemCrossValidationIteratorBase::
set_num_equations_per_point(int num_eq){
  numEquationsPerPoint_=num_eq;
}

void LinearSystemCrossValidationIteratorBase::
parse_options(const RealMatrix &A, const RealMatrix &B,
	      OptionsList &opts){
  set_num_points(opts.get<int>("num-points"));
  set_num_folds(opts.get("num-folds",std::min(10,numPts_)));
  set_seed(opts.get("seed",-1));
  int num_eq_per_pt = A.numRows()/numPts_;
  if (A.numRows()%numPts_!=0){
    std::string msg = "run: num rows of A must be a integer multiple of num pts";
    throw(std::runtime_error(msg));
  }
  set_num_equations_per_point(num_eq_per_pt);
  
  if ( numPts_ <= 0 )
    throw( std::runtime_error("run: num pts not set") );
  if ( A.numRows() != numPts_ * numEquationsPerPoint_ ){
    std::string msg;
    msg ="run: num pts and num equations per point are inconsistent with A";
    throw( std::runtime_error(msg) );
  }
  if ( B.numRows() != A.numRows() )
    throw( std::runtime_error("run: A and B are inconsistent") );
}

  
void LinearSystemCrossValidationIteratorBase::
get_scores(std::vector<RealVector> &result) const{
  result = scores_;
}

void LinearSystemCrossValidationIteratorBase::
get_fold_validation_residuals(std::vector< RealMatrix > &result) const{
  result = foldDiffs_;
};

void LinearSystemCrossValidationIteratorBase::
get_best_score_indices(IntVector &best_indices) const{
  int num_rhs=scores_.size(); 
  size_uninitialized(best_indices, num_rhs);
  for (int i=0; i<num_rhs; i++){
    best_indices[i] = 0;
    Real best_score = scores_[i][0];
    for (int j=1; j<scores_[i].length(); j++){
      if (scores_[i][j]<best_score){
	best_score = scores_[i][j];
	best_indices[i] = j;
      }
    }
  }
}
  
void LinearSystemCrossValidationIteratorBase::
get_best_scores(RealVector & scores) const{
  IntVector best_indices;
  get_best_score_indices(best_indices);
  size_uninitialized(scores, best_indices.length());
  for (int i=0; i<best_indices.length(); i++){
    scores[i] = scores_[i][best_indices[i]];
  }
}


LinearSystemCrossValidationIterator::LinearSystemCrossValidationIterator() :
  maxNumUniqueTols_(100) {};

LinearSystemCrossValidationIterator::~LinearSystemCrossValidationIterator(){};

  
void LinearSystemCrossValidationIterator::
set_linear_system_solver(const std::shared_ptr<LinearSystemSolver> &solver){
  linearSystemSolver_ = solver;
}

void LinearSystemCrossValidationIterator::
get_fold_tolerances(std::vector< RealVector > &result) const{
  result = foldTols_;
};

void LinearSystemCrossValidationIterator::
get_fold_scores(std::vector< RealVector > &result) const{
  result = foldScores_;
};

void LinearSystemCrossValidationIterator::
get_unique_tolerances(std::vector<RealVector> &result) const{
  result = uniqueTols_;
}
  
void LinearSystemCrossValidationIterator::
run(const RealMatrix &A, const RealMatrix &B, OptionsList &opts){

  parse_options(A, B, opts);
  create_partitions();

  if (!opts.isType<OptionsList>("regression-opts")){
    std::string msg = "Parameter List \"regression-opts\" is required.";
    throw(std::runtime_error(msg));
  }
  
  OptionsList regression_opts =  opts.get<OptionsList>("regression-opts");
 
  int num_rhs = B.numCols();
  scores_.resize(num_rhs);
  uniqueTols_.resize(num_rhs);
  for (int i=0; i<num_rhs; ++i){
    RealVector b = Teuchos::getCol(Teuchos::View, const_cast<RealMatrix &>(B), i);
    run_single_rhs(A, b, regression_opts);
    compute_scores(scores_[i], uniqueTols_[i]);
  }
}

void LinearSystemCrossValidationIterator::
run_single_rhs(const RealMatrix &A, const RealVector &b,
	       OptionsList &solver_opts){
  
  foldDiffs_.resize( num_folds() );
  foldTols_.resize( num_folds() );
  foldScores_.resize( num_folds() );
  RealMatrix A_train, A_valid;
  RealVector b_train, b_valid;
  IntVector training_indices, validation_indices;
  for (int iter=0; iter<num_folds(); iter++){
    // Must operate on B one column at a time because we may want to support each
    // column having different faults.
    // Right now we assume faults are the same for each of the RHS. If this
    // assumption is made permanent then we can increase efficiency by only
    // extracting A matrix once
    get_fold_indices( iter, training_indices, validation_indices );
    extract_matrix( A, training_indices, A_train );
    extract_matrix( A, validation_indices, A_valid );
    extract_values( b, training_indices, b_train );
    extract_values( b, validation_indices, b_valid );

    linearSystemSolver_->solve( A_train, b_train, solver_opts );

    RealMatrix coeff;
    RealVector residuals;
    linearSystemSolver_->get_solutions_for_all_regularization_params(coeff,0);
    linearSystemSolver_->get_residuals_for_all_regularization_params(residuals,0);
    int num_path_steps=coeff.numCols();

    // FIXME (BMA/JDJ): The following num_validation_primary_eqs
    // is incorrect in the case of mixed function and gradient
    // data with failures.  Want to compute cross validation
    // differences w.r.t. function values only, but that set may
    // be empty.  The number of valid function vs. gradient rows
    // needs to be tracked from remove_faulty_data.

    // TODO: Add test for above case

    // only keep values associated with primary equations.
    // assumes if faulty data exists then all data associated with 
    // the primary equation is removed. E.g. If the primary data
    // is a function value then it and all the gradients are removed
    // even if some gradients are fine.
    int num_validation_primary_eqs=
      b_valid.numRows() / numEquationsPerPoint_;

    foldDiffs_[iter].shapeUninitialized( num_validation_primary_eqs, 
					 num_path_steps );
    for (int i=0; i<num_validation_primary_eqs; i++){
      for (int j=0; j<num_path_steps; j++)
	foldDiffs_[iter](i,j)=b_valid[i];
    }

    // Speed up by only multiplying with columns of A that correspond
    // to non zero coeff
    for (int k=0; k<num_path_steps; k++){
      for (int j=0; j<A_valid.numCols(); j++){
	Real coeff_jk=coeff(j,k);
	Real *A_valid_j=A_valid[j], 
	  *fold_diffs_k=foldDiffs_[iter][k];
	if (std::abs( coeff_jk ) > std::numeric_limits<double>::epsilon()){
	  for (int i=0; i<num_validation_primary_eqs; i++)
	    fold_diffs_k[i] -= A_valid_j[i] * coeff_jk;
	}
      }
    }

    foldTols_[iter].shapeUninitialized( coeff.numCols(), 1 );
    for (int i=0; i<num_path_steps; i++)
      foldTols_[iter][i]=residuals[i]; 
    compute_fold_score( foldDiffs_[iter], foldScores_[iter] );
  }
}

void LinearSystemCrossValidationIterator::
set_max_num_unique_tolerances( int max_num_tols ){
  maxNumUniqueTols_ = max_num_tols;
}

  
void LinearSystemCrossValidationIterator::define_unique_tolerances(RealVector &unique_tols){
  int max_num_unique_tols = 0;
  std::set<Real> unique_tols_set;
  for ( int iter = 0; iter < num_folds(); iter++ ){
    int num_path_steps = foldDiffs_[iter].numCols();
    unique_tols_set.insert(foldTols_[iter].values(), 
			   foldTols_[iter].values() + 
			   foldTols_[iter].length());
    max_num_unique_tols = std::max( max_num_unique_tols, 
				    num_path_steps );
  }

  max_num_unique_tols = std::min( max_num_unique_tols, maxNumUniqueTols_ );
  // There are often thousands of unique parameter values so take 
  // a thinned subset of these. The size of the subset is controlled by
  // max_num_unique_tols

  int num_unique_tols = unique_tols_set.size();
  int stride = num_unique_tols / max_num_unique_tols;
  num_unique_tols = unique_tols_set.size() / stride;
  if (  unique_tols_set.size() % stride != 0 )
    num_unique_tols += 1;

  unique_tols.sizeUninitialized( num_unique_tols );
  std::set<Real>::iterator it = unique_tols_set.begin();
  int i = 0, j = 0;
  for (it=unique_tols_set.begin(); it!=unique_tols_set.end(); ++it){
    if ( j % stride == 0 ){
      unique_tols[i] = *it; 
      i++;
    }
    j++;
  }
}

void LinearSystemCrossValidationIterator::compute_scores(RealVector &scores,
							 RealVector &unique_tols){
  // Each path will have a different set of parameters. Consequently
  // we need to use interpolation to compute the parameters at 
  // a set of common points. These common points are contained in 
  // unique_tols

  define_unique_tolerances(unique_tols);  
  int num_unique_tols = unique_tols.length();
  RealMatrix unique_fold_scores( num_unique_tols, num_folds(), false );
  for (int iter=0; iter<num_folds(); iter++){
    // LinearInterpolant1D requires increasing x values
    // unique_tols_ is in ascending order
    // but foldTols_[iter] and foldScores_[iter] are in descending order
    // so we must reverse their order
    reverse( foldTols_[iter] );
    reverse( foldScores_[iter] );
    LinearInterpolant1D interp( foldTols_[iter], foldScores_[iter] );
    RealVector unique_fold_scores_col( Teuchos::View, 
				       unique_fold_scores[iter],
				       num_unique_tols );
    interp.interpolate( unique_tols, unique_fold_scores_col );
  }
  accumulate_fold_scores(unique_fold_scores, scores);
}

void LinearSystemCrossValidationIterator::
get_best_residual_tolerances(RealMatrix & residual_tols) const{
  IntVector best_indices;
  get_best_score_indices(best_indices);
  int num_rhs = best_indices.length();
  shape_uninitialized(residual_tols, num_rhs,1);
  for (int i=0; i<num_rhs; i++){
    residual_tols(i,0) = uniqueTols_[i][best_indices[i]];
  }
}

std::shared_ptr<LinearSystemSolver> LinearSystemCrossValidationIterator::
get_linear_system_solver(){
  return linearSystemSolver_;
}

void LinearSystemCrossValidationIterator::
generate_best_solutions(const RealMatrix &A, const RealMatrix &B,
                        RealMatrix &best_solutions, RealVector &residuals,
                        OptionsList & opts){
  if (!opts.isType<RealMatrix>("residual-tolerances")){
    RealMatrix best_residual_tolerances;
    get_adjusted_best_residual_tolerances(best_residual_tolerances);
    opts.set("residual-tolerances", best_residual_tolerances);
  }
  linearSystemSolver_->multi_rhs_solve(A, B, opts);
  linearSystemSolver_->get_final_solutions(best_solutions);
  linearSystemSolver_->get_final_residuals(residuals);
}



void LinearSystemCrossValidationIterator::
get_adjusted_best_residual_tolerances(RealMatrix &result) const{
  get_best_residual_tolerances(result);
  result *= (Real)num_folds()/((Real)num_folds()-1.0);
}

std::shared_ptr<LinearSystemCrossValidationIterator> cast_to_linear_system_cross_validation_iterator(std::shared_ptr<LinearSystemCrossValidationIteratorBase> &solver){
  std::shared_ptr<LinearSystemCrossValidationIterator> solver_cast =
    std::dynamic_pointer_cast<LinearSystemCrossValidationIterator>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to LinearSystemCrossValidationIterator shared_ptr"));
  return solver_cast;
}

}  // namespace util
}  // namespace Pecos
