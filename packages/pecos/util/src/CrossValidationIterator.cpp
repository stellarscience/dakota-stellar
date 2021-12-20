/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "CrossValidationIterator.hpp"

namespace Pecos {
namespace util {

CrossValidationIterator::CrossValidationIterator() : 
  numFolds_(0), numPts_(0), seed_(0), numEquationsPerPoint_(1) {}

CrossValidationIterator::~CrossValidationIterator(){ clear(); }

void CrossValidationIterator::get_point_indices(IntVector &result) const{
  result = indices_;
}

void CrossValidationIterator::set_num_folds(int num_folds){
  numFolds_ = num_folds;
  if (numPts_ == 0){
    std::string msg = "set_num_points() Please set numPts_";
    throw(std::runtime_error(msg));
  }
    
  if (numFolds_ > numPts_){
    std::string msg = "set_num_points() Ensure numFolds_ <= numPts_";
    throw(std::runtime_error(msg));
  }
}

void CrossValidationIterator::create_partitions(){
  
  if (failedRespData_.length()==0){
    // Set default of no failures by initializing to zero
    failedRespData_.size(numPts_);
  }
  
  int num_active_pts=0;
  for (int i=0; i<numPts_; ++i){
    if (!(failedRespData_[i] &1))
      num_active_pts++;
  }

  if (num_active_pts<numFolds_){
    std::string msg;
    msg = "Number of folds exceeds number of points after removing faults";
    throw(std::runtime_error(msg));
  }
  
  IntVector permutation_indices;
  if (seed_<0)
    range(permutation_indices, 0, num_active_pts, 1);
  else if (seed_ == 0)
    random_permutation(num_active_pts, 1, (unsigned int)std::time(0),
		       permutation_indices);
  else
    random_permutation(num_active_pts, 1, (unsigned int)seed_,
		       permutation_indices);
  int j=0;
  indices_.sizeUninitialized(num_active_pts);
  for (int i=0; i<numPts_; ++i){
    if (!(failedRespData_[i] &1)){
      indices_[permutation_indices[j]]=i;
      j++;
    }
  }
  
  foldStartingIndices_.sizeUninitialized(numFolds_);
  int max_fold_size = num_active_pts / numFolds_;
  if (num_active_pts % numFolds_ != 0) 
    max_fold_size++;
  foldStartingIndices_[0] = 0;
  for (int i=0; i<numFolds_-1; i++){
    int fold_size = max_fold_size-1;
    if ((i+1)*(max_fold_size) <= num_active_pts-((numFolds_-i-1)*(max_fold_size-1)))
      fold_size = max_fold_size;
    foldStartingIndices_[i+1] = foldStartingIndices_[i] + fold_size;
  }
  
}

void CrossValidationIterator::set_num_points(int num_points){
  numPts_ = num_points;
}

void CrossValidationIterator::set_seed(int seed){
  seed_ = seed;
}

void CrossValidationIterator::clear(){
  numFolds_ = 0; numPts_ = 0; indices_.sizeUninitialized(0);
  seed_ = 0;  numEquationsPerPoint_ = 0;
  failedRespData_.sizeUninitialized(0);
}

void CrossValidationIterator::copy(const CrossValidationIterator &source){
  set_num_folds(source.numFolds_); 
  seed_ = source.seed_;
  numPts_ = source.numPts_;
  indices_.sizeUninitialized(source.indices_.length());
  indices_.assign(source.indices_);
  numEquationsPerPoint_ = source.numEquationsPerPoint_;
}

void CrossValidationIterator::
get_fold_validation_indices(int iter, IntVector &validation_indices) const{
  int num_training_indices, num_validation_indices;
  get_fold_size(iter, num_training_indices, num_validation_indices);
  size_uninitialized(validation_indices, num_validation_indices);
  for (int j=0; j<num_validation_indices; j++)
    validation_indices[j] = indices_[foldStartingIndices_[iter]+j];
}

void CrossValidationIterator::
get_fold_training_indices(int iter, IntVector &training_indices) const{
  int num_training_indices, num_validation_indices;
  get_fold_size(iter, num_training_indices, num_validation_indices);
  int validation_end = foldStartingIndices_[iter] + num_validation_indices;
  size_uninitialized(training_indices, num_training_indices);
  int j=0;
  for (j=0; j<foldStartingIndices_[iter]; j++)
    training_indices[j] = indices_[j];
  int num_active_pts = indices_.length();
  for (int k=0; k<num_active_pts-validation_end; k++)
    training_indices[j+k] = indices_[validation_end+k];
}

void CrossValidationIterator::
get_fold_indices(int iter, 
		 IntVector &training_indices, 
		 IntVector &validation_indices) const{
  //int num_training_indices, num_validation_indices;
  get_fold_validation_indices(iter, validation_indices);
  get_fold_training_indices(iter, training_indices);
}

int CrossValidationIterator::num_folds() const{
  return numFolds_;
}

int CrossValidationIterator::num_points() const{
  return numPts_;
}

void CrossValidationIterator::
get_fold_size(int iter, int &training_size, int& validation_size) const{
  int num_active_pts = indices_.length();
  if (iter<numFolds_ - 1)
    validation_size = foldStartingIndices_[iter+1] - foldStartingIndices_[iter];
  else
    validation_size = num_active_pts - foldStartingIndices_[iter];
  training_size = num_active_pts - validation_size;
}
  
void CrossValidationIterator::
extract_values(const RealMatrix &B,
	       const IntVector &indices, RealMatrix &result_0) const{

  if (B.numRows() != numPts_ * numEquationsPerPoint_)
    throw(std::runtime_error("extract_values: num pts and num equations per point are inconsistent with b"));

  int num_indices = indices.length(), num_rhs = B.numCols();
  int num_indices_with_secondary_eq = 0;
  for (int i=0; i<num_indices; i++){
    if (!failedRespData_[indices[i]])
      num_indices_with_secondary_eq++;
  }
  int num_submatrix_rows = indices.length() +
    (numEquationsPerPoint_-1)*num_indices_with_secondary_eq;

  shape_uninitialized(result_0, num_submatrix_rows, num_rhs);
  for (int j=0; j<num_rhs; j++){
    int k = 0;
    for (int i=0; i<num_indices; i++){
      // assumes there is a set of primary equations stored 0...b.numRows()-1
      // in b. Then the non-primary equations are stored 
      // 1...numEquationsPerPoint_ for first point then all non-primary 
      // equations for second point and so on.
      result_0(i,j) = B(indices[i],j);
      if (!failedRespData_[indices[i]]){
	int shift = numPts_ + indices[i] * (numEquationsPerPoint_-1);
	for (int m=0; m<numEquationsPerPoint_-1; m++){
	  result_0(num_indices+k,j) = B(shift+m,j);
	  k++;
	}
      }
    }
  }
}
  
void CrossValidationIterator::extract_points(const RealMatrix &points,
					     const IntVector &training_indices, 
					     RealMatrix &result_0) const{
  shape_uninitialized(result_0, points.numRows(), training_indices.length());
  for (int j=0; j<training_indices.length(); j++)
    for (int i=0; i<points.numRows(); i++)
      result_0(i,j) = points(i,training_indices[j]);
}

  
void CrossValidationIterator::
compute_fold_score(const RealMatrix &fold_diffs, RealVector &result) const{
  size_uninitialized(result, fold_diffs.numCols());
  result=0.;
  for ( int n = 0; n < fold_diffs.numCols(); n++ ){
    for ( int m = 0; m < fold_diffs.numRows(); m++ )
      result[n] += fold_diffs(m,n) * fold_diffs(m,n);
  }
}

   
void CrossValidationIterator::
accumulate_fold_scores(const RealMatrix &fold_scores, RealVector &scores) const{
  int num_scores_per_fold = fold_scores.numRows();
  size_uninitialized(scores, num_scores_per_fold);
  scores=0.;
  for (int i=0; i<num_scores_per_fold; i++){
    for (int iter=0; iter<num_folds(); iter++)
      scores[i] += fold_scores(i,iter);
    scores[i] /= num_points();
  }
}

void CrossValidationIterator::
set_fault_data(const IntVector &failed_response_data){
  if (failed_response_data.length()!=numPts_){
    std::string msg = "failed response data is not consistent with numPts_";
    throw(std::runtime_error(msg));
  }
  failedRespData_ = failed_response_data;
}
  
  /*
void run(const RealMatrix &points, const RealVector &values,
	 SurrogateBuilder &builder, OptionsList &build_opts){
  RealMatrix pts_train, pts_valid;
  RealVector value;
  for (int iter=0; iter<num_folds(); iter++){
    extract_points(points, training_indices, pts_train);
    extract_points(points, validation_indices, pts_valid);
    extract_values(b, validation_indices, b_valid);
    OptionsList result;
    std::shared_ptr<Approximation> approx =
      builder->build(pts_train, b_train, build_opts, result);
    approx.value(pts_valid,approx_values);
    for (int i=0; i<values_valid.length(); i++){
      foldDiffs_[iter](i,0) = values_valid[i] - approx_values[i];
    }
  }
  }*/

}  // namespace util
}  // namespace Pecos
