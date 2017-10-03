#ifndef CROSS_VALIDATION_HPP
#define CROSS_VALIDATION_HPP

#include "MathTools.hpp"
#include "RandomNumberGenerator.hpp"
#include "LinearSolver.hpp"
#include "RuntimeEnvironment.hpp"
#include "FaultTolerance.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {

class CrossValidationIterator : public ParallelObject
{
protected:
  int numFolds_;
  
  int numPts_;

  IntVector foldStartingIndices_;

  IntVector indices_;

  int seed_;

  int dataType_;

  int numEquationsPerPoint_;

  RandomNumberGenerator RNG_;
  
  // Dakota specific member variables
  bool faultInfoActive_;
  FaultInfo faultInfo_;
  SizetShortMap failedRespData_;

public:
  
  CrossValidationIterator();

  ~CrossValidationIterator();

  void set_num_folds( int num_folds );

  void set_num_points( int num_points );

  void set_num_equations_per_point( int num_eq );

  void set_seed( int seed );

  int num_folds();

  int num_pts();

  void clear();

  void get_point_indices( IntVector &result );

  void get_fold_indices( int iter, IntVector &result_0,
			 IntVector &result_1 );

  void copy( const CrossValidationIterator &source );

  void set_data_type( int data_type );

  void get_fold_size( int iter, int &training_size, int& validation_size );
  
  void extract_values( RealVector &b,
		       IntVector &indices,
		       RealVector &result_0 );

  void extract_matrix( RealMatrix &A,
		       IntVector &indices, 
		       RealMatrix &result );

  void extract_linear_system( RealMatrix &A, RealVector &b,
			      IntVector &indices, 
			      RealMatrix &result_0,
			      RealVector &result_1 );

  void extract_points( RealMatrix &points,
		       IntVector &training_indices, 
		       RealMatrix &result_0 );

  // dakota specific functions
  void set_fault_data( FaultInfo &fault_info,
		       const SizetShortMap& failed_resp_data )
  {
    faultInfoActive_ = true;
    faultInfo_ = fault_info;
    failedRespData_ = failed_resp_data;
  };
};

class LinearModelCrossValidationIterator : public CrossValidationIterator
{
protected:

  LinearSolver_ptr solver_;

  RealVector scores_;

  RealVector stdErrors_;

  RealVector uniqueTols_;

  /// The difference between the validation value and the approximation 
  /// at each validation point and for each solution generated on the 
  /// solution path
  std::vector< RealMatrix > foldDiffs_;

  /// The solver tolerances corresponding to each solution generated on the 
  /// solution path
  std::vector< RealVector > foldTols_;

  /// The errors ( a function of foldDiffs_ ) for each solution generated on the 
  /// solution path
  std::vector< RealVector > foldErrors_;

  /// The residual tolerance that minimizes the cross validation score
  Real bestResidualTol_;

  /// statistics of the coefficients generated on the solution path. This is
  /// useful if using cross validation to estimate variance in a moment, e.g
  /// the variance of the mean of a PCE 
  std::vector< RealMatrix > foldCoefficientStats_;

  /// statistics of the coefficients generated corresponding to the best
  /// residual chosen by cross validation.
  RealMatrix coefficientStats_;
  
public:

  LinearModelCrossValidationIterator() : bestResidualTol_( 0.0 ) {};

  ~LinearModelCrossValidationIterator()
  {
    foldDiffs_.clear(); foldTols_.clear(); bestResidualTol_ = 0.;
  }

  void get_scores( RealVector &result );

  void get_std_errors( RealVector &result );

  void get_unique_tolerances( RealVector &result );

  void get_fold_errors( std::vector< RealMatrix > &result );

  void get_fold_tolerances( std::vector< RealVector > &result );

  Real get_best_residual_tolerance();
  
  virtual Real run_cross_validation( RealMatrix &A, RealVector &b ) = 0;


  virtual void set_solver( LinearSolver_ptr solver );

  LinearSolver_ptr get_solver();


  void copy_solver( LinearSolver_ptr solver );

  /// Copy only the solver of the cross validation iterator;
  virtual boost::shared_ptr<LinearModelCrossValidationIterator> copy() = 0;

  virtual void compute_fold_score( RealMatrix &fold_diffs, RealVector &result );

};

typedef boost::shared_ptr<LinearModelCrossValidationIterator>  LinearModelCVIterator_ptr;

class MultipleSolutionLinearModelCrossValidationIterator : public LinearModelCrossValidationIterator
{
protected:

  int maxNumUniqueTols_;

public:
  
  MultipleSolutionLinearModelCrossValidationIterator() : 
    maxNumUniqueTols_( std::numeric_limits<int>::max() ) {};

  ~MultipleSolutionLinearModelCrossValidationIterator()
  {
    maxNumUniqueTols_=std::numeric_limits<int>::max();
  };

  void set_max_num_unique_tolerances( int max_num_tols );

  void collect_fold_data();

  void define_unique_tolerances();

  void compute_scores();

  Real run_cross_validation( RealMatrix &A, RealVector &b );

  void get_coefficient_stats(RealMatrix &coeff_stats);

   /// Copy only the solver of the cross validation iterator;
  boost::shared_ptr<LinearModelCrossValidationIterator> copy();
};

}  // namespace Pecos

#endif //CROSS_VALIDATION_HPP
