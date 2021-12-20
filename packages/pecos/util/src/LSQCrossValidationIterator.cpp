/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "LSQCrossValidationIterator.hpp"

namespace Pecos {
namespace util {

void LSQCrossValidationIterator::
lsq_lko_cross_validation(const RealMatrix &A, const RealMatrix &AtA_inv,
			 const RealMatrix &residuals,
			 std::vector<RealMatrix> &fold_diffs) {
  
  int num_rows = residuals.numRows(), num_rhs = residuals.numCols();
  fold_diffs.resize(num_folds());
  
  IntVector validation_indices;
  RealMatrix cv_diffs, A_valid, H_valid, tmp, validation_residuals;
  for (int iter=0; iter<num_folds(); iter++){
    get_fold_validation_indices(iter, validation_indices);
    int num_validation_indices = validation_indices.length();
    extract_matrix(A, validation_indices, A_valid);
    extract_values(residuals, validation_indices, validation_residuals);
    
    // Compute A_valid * (A'A)_inv
    shape_uninitialized(tmp, A_valid.numRows(), AtA_inv.numCols());
    tmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A_valid, AtA_inv,0.0);
    
    // Compute - A_valid * (A'A)_inv * A_valid.T
    shape_uninitialized(H_valid, A_valid.numRows(),A_valid.numRows());
    H_valid.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, -1.0, tmp, A_valid, 0.0);
    // Compute H = I - A_valid * (A'A)_inv * A_valid.T
    for (int i=0; i<A_valid.numRows(); i++)
      H_valid(i,i) += 1.0;

    RealMatrix chol_factor, H_valid_inv;
    cholesky(H_valid, chol_factor, Teuchos::LOWER_TRI, true);
    cholesky_inverse(chol_factor, H_valid_inv, Teuchos::LOWER_TRI);

    // WARNING when activating faults num_validation_indices will not
    // equal A_valid.numRows(). Want to compute diffs for all rows of A
    // then down select to primary equations
    fold_diffs[iter].shapeUninitialized(num_validation_indices, num_rhs);
    fold_diffs[iter].multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0,
			      H_valid_inv, validation_residuals, 0.0);
    }
}

void LSQCrossValidationIterator::
lsq_loo_cross_validation(const RealMatrix &A, const RealMatrix &AtA_inv,
			 const RealMatrix &residuals, 
			 std::vector<RealMatrix> &fold_diffs){
  
  int num_rows = residuals.numRows(), num_rhs = residuals.numCols();
  fold_diffs.resize( num_rows );
  
  // Compute A * (A'A)_inv
  RealMatrix tmp( num_rows, AtA_inv.numCols(), false );
  tmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, AtA_inv, 0.0);
  // Compute A * (A'A)_inv * A'
  RealMatrix H(num_rows, num_rows, false);  
  H.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, tmp, A, 0.0);
  // Compute cross validation differences
  IntVector validation_indices;
  for (int iter = 0; iter<num_rows; iter++){
    get_fold_validation_indices(iter, validation_indices);
    fold_diffs[iter].shapeUninitialized(1, num_rhs);
    for (int n=0; n<num_rhs; n++){
      int index = validation_indices[0];
      fold_diffs[iter](0,n) = residuals(index,n) / (1. - H(index,index));
    }
  }
}

void LSQCrossValidationIterator::
run(const RealMatrix &A, const RealMatrix &B, OptionsList& opts){

  parse_options(A, B, opts);
  create_partitions();
  
  int num_rows = A.numRows(), num_cols = A.numCols();

  // TODO: Make work with faults
  // Extract full A and B without faults
  // Compute X and residuals using faultless A and B
  // pass faultless AtA_inv to loo and lko functions
  // compute cross validation scores only on primary equations
  
  // Compute A't*B 
  RealMatrix AtB(A.numCols(), B.numCols(), false);
  AtB.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, B, 0.0);
  
  // Compute A'A
  RealMatrix AtA(A.numCols(), A.numCols(), false);
  AtA.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, A, 0.0);

  // Compute (A'A)_inv
  RealMatrix chol_factor, AtA_inv;
  cholesky(AtA, chol_factor, Teuchos::LOWER_TRI, true);
  cholesky_inverse(chol_factor, AtA_inv, Teuchos::LOWER_TRI);

  // Solve the linear system Ax=b
  coeffs_.shapeUninitialized(num_cols,B.numCols());
  coeffs_.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, 
	     AtA_inv, AtB, 0.0);
  
  // Compute residual at each build point
  RealMatrix residuals(num_rows, B.numCols(), false);
  residuals.assign(B);
  residuals.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0, A,coeffs_,1.0);

  // Compute cross validation differences
  if (num_folds()==num_points())
    lsq_loo_cross_validation(A, AtA_inv, residuals, foldDiffs_);
  else
    lsq_lko_cross_validation(A, AtA_inv, residuals, foldDiffs_);

  int num_rhs = B.numCols();
  scores_.resize(num_rhs);
  RealMatrix fold_scores(1, num_folds(), false);
  for (int n=0; n<num_rhs; n++){
    fold_scores = 0.0;
    for (int iter=0; iter<num_folds(); iter++){
      for (int i=0; i<foldDiffs_[iter].numRows(); i++){
	fold_scores(0,iter) += foldDiffs_[iter](i,n) * foldDiffs_[iter](i,n);
      }
    }
    accumulate_fold_scores(fold_scores, scores_[n]);
  }

  size_uninitialized(residualNorms_, num_rhs);
  for (int n=0; n<num_rhs; n++){
    RealVector residuals_col = Teuchos::getCol(Teuchos::View, residuals, n);
    residualNorms_[n] = residuals_col.normFrobenius();
  }
}

void LSQCrossValidationIterator::
get_coefficients(RealMatrix &coeffs) const{
  coeffs = coeffs_;
}

void LSQCrossValidationIterator::
get_best_scores(RealVector & scores) const{
  int num_rhs = (int)scores_.size();
  size_uninitialized(scores, num_rhs); 
  for (int i=0; i<num_rhs; ++i)
    scores[i] = scores_[i][0];
}

void LSQCrossValidationIterator::
generate_best_solutions(const RealMatrix &A, const RealMatrix &B,
			RealMatrix &best_solutions,
			RealVector &residuals,
			OptionsList & opts){
  get_coefficients(best_solutions);
  residuals = residualNorms_;
}

std::shared_ptr<LSQCrossValidationIterator> cast_to_least_squares_cross_validation_iterator(std::shared_ptr<LinearSystemCrossValidationIteratorBase> &solver){
  std::shared_ptr<LSQCrossValidationIterator> solver_cast =
  std::dynamic_pointer_cast<LSQCrossValidationIterator>(solver);
  if (!solver_cast)
    throw(std::runtime_error("Could not cast to LSQCrossValidationIterator shared_ptr"));
  return solver_cast;
}
  

}  // namespace util
}  // namespace Pecos
