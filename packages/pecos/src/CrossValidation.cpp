#include "CrossValidation.hpp"
#include "MathTools.hpp"

namespace Pecos {

CrossValidationIterator::CrossValidationIterator() : 
  numFolds_( 0 ), numPts_( 0 ), seed_( 0 ), dataType_( 0 ),
  numEquationsPerPoint_( 0 ), faultInfoActive_(false) {};

CrossValidationIterator::~CrossValidationIterator()
{ clear(); }

void CrossValidationIterator::get_point_indices( IntVector &result )
{
  result = indices_;
}

void CrossValidationIterator::set_num_folds( int num_folds )
{
  numFolds_ = num_folds;
};

void CrossValidationIterator::set_num_points( int num_points )
{
  numPts_ = num_points;

  if ( numFolds_ > numPts_ )
    {
      std::string msg = "set_num_points() Ensure numFolds_ <= numPts_";
      throw( std::runtime_error( msg ) );
    }
  if ( numFolds_ == 0 )
    {
      std::string msg = "set_num_points() Please set numFolds_";
      throw( std::runtime_error( msg ) );
    }

  foldStartingIndices_.sizeUninitialized( numFolds_ );
  int max_fold_size = numPts_ / numFolds_;
  if ( numPts_ % numFolds_ != 0 ) 
    max_fold_size++;
  foldStartingIndices_[0] = 0;
  for ( int i = 0; i < numFolds_-1; i++ )
    {
      int fold_size = max_fold_size-1;
      if ( (i+1) *  ( max_fold_size ) <= 
	   numPts_ - ( (numFolds_-i-1)*(max_fold_size-1) ) )
	fold_size = max_fold_size;
      foldStartingIndices_[i+1] = foldStartingIndices_[i] + fold_size;
    }

  if ( seed_ < 0 )
    range( indices_, 0, numPts_, 1 );
  else if ( seed_ == 0 )
    RNG_.permutation( indices_, numPts_, 1, (unsigned int)std::time(0) );
  else
    //get_permutations( indices_, numPts_, 1, (unsigned int)seed_ );
    RNG_.permutation( indices_, numPts_, 1, (unsigned int)seed_ );
};

void CrossValidationIterator::set_num_equations_per_point( int num_eq )
{
  numEquationsPerPoint_ = num_eq;
}

void CrossValidationIterator::set_seed( int seed )
{
  seed_ = seed;
  // if permutations have already computed then recompute
  if ( numPts_ > 0 )
    set_num_points( numPts_ );
}

void CrossValidationIterator::clear()
{
  numFolds_ = 0; numPts_ = 0; indices_.sizeUninitialized( 0 );
  seed_ = 0;  numEquationsPerPoint_ = 0;
};

void CrossValidationIterator::copy( const CrossValidationIterator &source )
{
  set_num_folds( source.numFolds_ ); 
  // shallow copy indices
  seed_ = source.seed_;
  numPts_ = source.numPts_;
  indices_ = source.indices_;
  numEquationsPerPoint_ = source.numEquationsPerPoint_;
};

void CrossValidationIterator::get_fold_indices( int iter, 
					        IntVector &training_indices, 
						IntVector &validation_indices )
{
  int num_training_indices, num_validation_indices;
  get_fold_size( iter, num_training_indices, num_validation_indices );
  validation_indices.sizeUninitialized( num_validation_indices );
  for ( int j = 0; j < num_validation_indices; j++ )
    validation_indices[j] = indices_[foldStartingIndices_[iter]+j];
  
  int validation_end = foldStartingIndices_[iter] + num_validation_indices;
  training_indices.sizeUninitialized( num_training_indices );
  int j = 0;
  for ( j = 0; j < foldStartingIndices_[iter]; j++ )
    training_indices[j] = indices_[j];
  for ( int k = 0; k < numPts_-validation_end; k++ )
    training_indices[j+k] = indices_[validation_end+k];
};

int CrossValidationIterator::num_folds()
{
  return numFolds_;
};

int CrossValidationIterator::num_pts()
{
  return numPts_;
}

void CrossValidationIterator::set_data_type( int data_type )
{
  dataType_ = data_type;
}

void CrossValidationIterator::get_fold_size( int iter, int &training_size, 
					     int& validation_size )
{
  if ( iter < numFolds_ - 1 )
    validation_size = foldStartingIndices_[iter+1] - foldStartingIndices_[iter];
  else
    validation_size = numPts_ - foldStartingIndices_[iter];
  training_size = numPts_ - validation_size;
}
  
void CrossValidationIterator::extract_values( RealVector &b,
					      IntVector &indices,
					      RealVector &result_0 )
{

  if ( b.numRows() != numPts_ * numEquationsPerPoint_ )
    throw( std::runtime_error("extract_values: num pts and num equations per point are inconsistent with b") );

  int k = 0;
  int num_indices = indices.length();
  result_0.sizeUninitialized( num_indices * numEquationsPerPoint_  );
  for ( int i = 0; i < num_indices; i++ )
    {
      // assumes there is a set of primary equations stored 0...b.length()-1
      // in b. Then the non-primary equations are stored 
      // 1...numEquationsPerPoint_ for first point then all non-primary 
      // equations for second point and so on.
      result_0[i] = b[indices[i]];
      int shift = numPts_ + indices[i] * (numEquationsPerPoint_-1);
      for ( int m = 0; m < numEquationsPerPoint_-1; m++ )
	{
	  result_0[num_indices+k] = b[shift+m];
	  k++;
	}
    }
}

void CrossValidationIterator::extract_matrix( RealMatrix &A,
					      IntVector &indices, 
					      RealMatrix &result )
{
  if ( A.numRows() != numPts_ * numEquationsPerPoint_ )
    throw( std::runtime_error("extract_matrix: num pts and num equations per point are inconsistent with A") );

  int num_indices = indices.length();
  result.shapeUninitialized( indices.length() * numEquationsPerPoint_, 
			     A.numCols() );
  for ( int j = 0; j < A.numCols(); j++ )
    {
      int k = 0;
      for ( int i = 0; i < num_indices; i++ )
	{
	  // assumes there is a set of primary equations stored 
	  // 0...b.length()-1
	  // in b. Then the non-primary equations are stored 
	  // 1...numEquationsPerPoint_ for first point then all non-primary 
	  // equations for second point and so on.
	  result(i,j) = A(indices[i],j);
	  int shift = numPts_ + indices[i] * (numEquationsPerPoint_-1);
	  for ( int m = 0; m < numEquationsPerPoint_-1; m++ )
	    {
	      result(num_indices+k,j) = A(shift+m,j);
	      k++;
	    }
	}
    }
}

void CrossValidationIterator::extract_linear_system( RealMatrix &A, 
						     RealVector &b,
						     IntVector &indices, 
						     RealMatrix &result_0,
						     RealVector &result_1 )
{
  extract_matrix( A, indices, result_0 );
  result_1.sizeUninitialized( indices.length() );
  extract_values( b, indices, result_1 );
}

void CrossValidationIterator::extract_points( RealMatrix &points,
					      IntVector &training_indices, 
					      RealMatrix &result_0 )
{
  result_0.shapeUninitialized( points.numRows(),training_indices.length());
  for ( int j = 0; j < training_indices.length(); j++ )
    for ( int i = 0; i < points.numRows(); i++ )
      result_0(i,j) = points(i,training_indices[j]);
}

void LinearModelCrossValidationIterator::get_scores( RealVector &result )
{
  result = scores_;
}

void LinearModelCrossValidationIterator::get_std_errors( RealVector &result )
{
  result = stdErrors_;
}

void LinearModelCrossValidationIterator::get_unique_tolerances( RealVector &result )
{
  result = uniqueTols_;
}

void LinearModelCrossValidationIterator::get_fold_errors( std::vector< RealMatrix > &result )
{
  result = foldDiffs_;
};

void LinearModelCrossValidationIterator::get_fold_tolerances( std::vector< RealVector > &result )
{
  result = foldTols_;
};

Real LinearModelCrossValidationIterator::get_best_residual_tolerance()
{
  return bestResidualTol_;
}

void LinearModelCrossValidationIterator::set_solver( LinearSolver_ptr solver )
{
  solver_ = solver;
};

LinearSolver_ptr LinearModelCrossValidationIterator::get_solver()
{
  return solver_;
};

void LinearModelCrossValidationIterator::copy_solver( LinearSolver_ptr solver )
{
  if ( !solver )
    {
      std::string msg = "copy_solver() source is an empty pointer";
      throw( std::runtime_error( msg ) );
    }
  solver_->copy( *solver );
};

void LinearModelCrossValidationIterator::compute_fold_score( 
							    RealMatrix &fold_diffs, RealVector &result )
{
  result.size( fold_diffs.numCols() );
  for ( int n = 0; n < fold_diffs.numCols(); n++ )
    {
      for ( int m = 0; m < fold_diffs.numRows(); m++ )
	result[n] += fold_diffs(m,n) * fold_diffs(m,n);
      //result[n] = std::sqrt( result[n] / (double)fold_diffs.numRows() );
    }
};

void MultipleSolutionLinearModelCrossValidationIterator::set_max_num_unique_tolerances( int max_num_tols )
{
  maxNumUniqueTols_ = max_num_tols;
}

void MultipleSolutionLinearModelCrossValidationIterator::collect_fold_data()
{
  //std::cout << "Processor " << processor_id() << " collecting data\n";
  if ( !is_master() )
    { 
#ifdef HEAT_ENABLE_MPI
      send( foldTols_, master_processor_id(), MPICommunicator_ ); 
      send( foldErrors_, master_processor_id(), MPICommunicator_ );
      send( foldCoefficientStats_, master_processor_id(), MPICommunicator_ );
#endif
    }
  else 
    { 
      for ( int p = 0; p < num_processors(); p++ )
	{
	  if ( p != master_processor_id() )
	    {
	      std::vector<RealVector> fold_tols, fold_errors;
	      std::vector<RealMatrix> fold_coefficient_stats;
	      
	      //std::cout << "Master process receiving data from processor "
	      //	  << p << std::endl;

#ifdef HEAT_ENABLE_MPI //hack 
	      receive( fold_tols, p, MPICommunicator_ ); 
	      receive( fold_errors, p, MPICommunicator_ );
	      receive( fold_coefficient_stats, p, MPICommunicator_ );
#endif
	      for ( int iter = 0; iter < numFolds_; iter++ )
		{
		  if ( ( (iter+1) % num_processors() ) == p )
		    {
		      foldTols_[iter] = fold_tols[iter];
		      foldErrors_[iter] = fold_errors[iter];
		      foldCoefficientStats_[iter] = fold_coefficient_stats[iter];
		    }
		}
	    }
	}
    }
  //std::cout << "Processor " << processor_id() << " finished collecting data\n";
};

void MultipleSolutionLinearModelCrossValidationIterator::define_unique_tolerances()
{
  int max_num_unique_tols = 0;
  std::set<Real> unique_tols_set;
  if ( is_master() )
    {
      for ( int iter = 0; iter < num_folds(); iter++ )
	{
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

      uniqueTols_.sizeUninitialized( num_unique_tols );
      std::set<Real>::iterator it = unique_tols_set.begin();
      int i = 0, j = 0;
      for ( it = unique_tols_set.begin(); it != unique_tols_set.end(); ++it )
	{
	  if ( j % stride == 0 )
	    {
	      uniqueTols_[i] = *it; 
	      i++;
	    }
	  j++;
	}
    }
}

void MultipleSolutionLinearModelCrossValidationIterator::compute_scores()
{
  if ( is_master() )
    {
      // Each path will have a different set of parameters. Consequently
      // we need to use interpolation to compute the parameters at 
      // a set of common points. These commone points are contained in 
      // uniqueTols_

      define_unique_tolerances();
    
      int num_unique_tols = uniqueTols_.length();
      RealMatrix unique_cv_errors( num_unique_tols, num_folds(), false );
      for ( int iter = 0; iter < num_folds(); iter++ )
	{
	  // LinearInterpolant1D requires increasing x values
	  // uniqueTols_ is in ascending order
	  // but foldTols_[iter] and foldErrors_[iter] are in descending order
	  // so we must reverse their order
	  reverse( foldTols_[iter] );
	  reverse( foldErrors_[iter] );
	  LinearInterpolant1D interp( foldTols_[iter], foldErrors_[iter] );
	  RealVector unique_cv_errors_col( Teuchos::View, 
					   unique_cv_errors[iter],
					   num_unique_tols );
	  interp.interpolate( uniqueTols_, unique_cv_errors_col );
	}

      scores_.size( num_unique_tols );
      stdErrors_.sizeUninitialized( num_unique_tols );
      for ( int i = 0; i < num_unique_tols; i++ )
	{
	  Real s2 = 0., s1 = 0.;
	  for ( int iter = 0; iter < num_folds(); iter++ )
	    {
	      int training_size, validation_size;
	      get_fold_size( iter, training_size, validation_size );
	      scores_[i] += unique_cv_errors(i,iter);
	      Real incr = ( unique_cv_errors(i,iter) / (Real)validation_size );
	      s2 += incr * incr;
	      s1 += incr;
	    }
	  scores_[i] /= num_pts();
	  // sd = sqrt((sum x^2 - barx^2*n)/(n-1))
	  stdErrors_[i] = std::sqrt( ( s2 - s1*s1/(Real)num_folds()) /
				     (Real)(num_folds()-1) );
	  // se = sd / n
	  stdErrors_[i] = stdErrors_[i] / std::sqrt( (Real)num_folds() );
	}
    }
}

Real MultipleSolutionLinearModelCrossValidationIterator::run_cross_validation( RealMatrix &A, RealVector &b )
{
  if ( !solver_ )
    {
      std::string msg = "run_cross_validation: Ensure set_solver() ";
      msg += "has been envoked()";
      throw( std::runtime_error( msg ) );
    }

  // must clear solver because it may contain residualTol_ information
  // from a previous run
  // solver_->set_residual_tolerance( 0.0 );
  // Now the above is obsolete must best tol is stored in iterator

  //set_num_points( b.numRows() );
  if ( numPts_ <= 0 )
    throw( std::runtime_error("run_cross_validation: num pts not set") );
  if ( ( dataType_ == 0 ) && ( A.numRows() != numPts_ * numEquationsPerPoint_ ) )
    throw( std::runtime_error("run_cross_validation: num pts and num equations per point are inconsistent with A") );
    
  foldDiffs_.resize( num_folds() );
  foldTols_.resize( num_folds() );
  foldErrors_.resize( num_folds() );
  foldCoefficientStats_.resize(num_folds());
  for ( int iter = 0; iter < num_folds(); iter++ )
    {
      if ( ( (iter+1) % num_processors() ) == processor_id() ) 
	{
	  RealMatrix A_train, A_valid;
	  RealVector b_train, b_valid;
	  IntVector training_indices, validation_indices;
	  get_fold_indices( iter, training_indices, validation_indices );
	  extract_values( b, training_indices, b_train );
	  extract_values( b, validation_indices, b_valid );
	  //int num_validation_indices = validation_indices.length();
	  RealMatrix coeff, metrics;
	  if ( dataType_ == 0 )
	    {
	      // A is linear system
	      extract_matrix( A, training_indices, A_train );
	      extract_matrix( A, validation_indices, A_valid );

	      RealMatrix points_dummy;
	      if (faultInfoActive_) {
		remove_faulty_data( A_train, b_train, points_dummy, 
				    training_indices,
				    faultInfo_, failedRespData_ );
		remove_faulty_data( A_valid, b_valid, points_dummy, 
				    validation_indices,
				    faultInfo_, failedRespData_ );
	      }
	      solver_->solve( A_train, b_train, coeff, metrics );
	    }
	  else
	    {
	      // A is coordinates of build points
	      RealMatrix pts_train;
	      extract_points( A, training_indices, pts_train );
	      solver_->solve_using_points( pts_train, b_train, 
					   coeff, metrics );
	      // construct A_valid
	      RealMatrix A_all;
	      // dont forget A is points when dataType_ == 0
	      solver_->build_matrix( A, A_all );
	      extract_matrix( A_all, validation_indices, A_valid );
	    }

	  int num_path_steps = coeff.numCols();

	  foldCoefficientStats_[iter].shapeUninitialized(num_path_steps,
							   1);
	  for ( int j = 0; j < num_path_steps; j++ )
	    foldCoefficientStats_[iter](j,0) = coeff(0,j);

	  // FIXME (BMA/JDJ): The following num_validation_primary_eqs
	  // is incorrect in the case of mixed function and gradient
	  // data with failures.  Want to compute cross validation
	  // differences w.r.t. function values only, but that set may
	  // be empty.  The number of valid function vs. gradient rows
	  // needs to be tracked from remove_faulty_data.

	  // only keep values associated with primary equations.
	  // assumes if faulty data exists then all data associated with 
	  // the primary equation is removed. E.g. If the primary data
	  // is a function value then it and all the gradients are removed
	  // even if some gradients are fine.
	  int num_validation_primary_eqs = 
	    b_valid.numRows() / numEquationsPerPoint_;

	  foldDiffs_[iter].shapeUninitialized( num_validation_primary_eqs, 
					       num_path_steps );
	  for ( int i = 0; i < num_validation_primary_eqs; i++ ){
	    for ( int j = 0; j < num_path_steps; j++ )
	      foldDiffs_[iter](i,j) = b_valid[i];
	  }

	  // Speed up by only multiplying with columns of A that correspond
	  // to non zero coeff
	  for ( int k = 0; k < num_path_steps; k++ )
	    {
	      for ( int j = 0; j < A_valid.numCols(); j++ )
		{
		  Real coeff_jk = coeff(j,k);
		  Real *A_valid_j = A_valid[j], 
		    *fold_diffs_k = foldDiffs_[iter][k];
		  if ( std::abs( coeff_jk ) > 
		       std::numeric_limits<double>::epsilon() )
		    {
		      for ( int i = 0; i < num_validation_primary_eqs; i++ )
			fold_diffs_k[i] -= A_valid_j[i] * coeff_jk;
		    }
		}
	    }

	  foldTols_[iter].shapeUninitialized( coeff.numCols(), 1 );
	  for ( int i = 0; i < num_path_steps; i++ )
	    foldTols_[iter][i] = metrics(0,i);

	  compute_fold_score( foldDiffs_[iter], foldErrors_[iter] );
	}	
    }

  collect_fold_data();
    
  compute_scores();
    
  if ( is_master() )
    {
      int argmin_index = 0;
      argmin_index = argmin( scores_.length(), scores_.values()  );
      bestResidualTol_ = uniqueTols_[argmin_index];
      int training_size, validation_size;
      get_fold_size( 0, training_size, validation_size );
      bestResidualTol_ *= std::sqrt( (Real)num_pts() / (Real)training_size );

      //solver_->set_residual_tolerance( best_tol );

      // for each fold find coefficient stat that corresponds to the
      // residual which is closest to the chosen residual. The chosen
      // residual is interpolated so the choice may not be the same for
      // each fold.
      coefficientStats_.shapeUninitialized(num_folds(),1);
      for ( int iter = 0; iter < num_folds(); iter++ ){
	int chosen_idx = -1;
	Real resid_distance = std::numeric_limits<Real>::max();
	for ( int i = 0; i < foldCoefficientStats_[iter].numRows(); i++ ){
	  Real dist=std::abs(
		foldCoefficientStats_[iter](i,0)-bestResidualTol_);
	  if (dist <= resid_distance){
	    resid_distance = dist;
	    chosen_idx = i;
	  }else{
	    coefficientStats_(iter,0) =
	      foldCoefficientStats_[iter](chosen_idx,0);
	    break;
	  }
	}
      }

      
      return scores_[argmin_index];
    }
  else return -1;
}

boost::shared_ptr<LinearModelCrossValidationIterator> MultipleSolutionLinearModelCrossValidationIterator::copy()
{
  boost::shared_ptr<LinearModelCrossValidationIterator> cv_iterator( new MultipleSolutionLinearModelCrossValidationIterator() );
  cv_iterator->CrossValidationIterator::copy((CrossValidationIterator)(*this));
  cv_iterator->copy_solver( solver_ );
    
  return cv_iterator;
};

void MultipleSolutionLinearModelCrossValidationIterator::get_coefficient_stats(RealMatrix &coeff_stats){
  coeff_stats = coefficientStats_;
}

} // namespace Pecos
