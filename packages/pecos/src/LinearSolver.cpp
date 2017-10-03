#include "LinearSolver.hpp"
#include "MathTools.hpp"

namespace Pecos {

void normalise_columns( RealMatrix &A, RealVector &result )
{
  int M = A.numRows(), N = A.numCols();
  result.sizeUninitialized( N );
  for ( int i = 0; i < N; i++ )
    {
      RealVector col( Teuchos::View, A[i], M );
      result[i] = col.normFrobenius();
      col.scale(1./result[i]);
    }
};

 void CompressedSensingTool::solve( RealMatrix &A, RealMatrix &B, 
				    RealMatrixArray &solutions, 
				    CompressedSensingOptions &opts,
				    CompressedSensingOptionsList &opts_list )
{
  int M( A.numRows() ), N( A.numCols() ), num_rhs( B.numCols() );

  set_linear_solver( opts );

  //  multiple solution matrices are needed when creating pce of gradients
  // this will cause B to have multiple columns
  solutions.resize( num_rhs );
  opts_list.resize( num_rhs );
  
  if (opts.solver == SVD_LEAST_SQ_REGRESSION )
    {
      RealMatrix single_solution;
      RealVector singular_values;
      int rank(0);
      svd_solve( A, B, single_solution, singular_values, rank, 
		 -1 );
      for ( int k = 0; k < num_rhs; k++ )
	{
	  reshape( solutions[k], N, 1 );
	  RealMatrix solution( Teuchos::View, single_solution, N, 1, 0, k );
	  solutions[k].assign( solution );

	  // Store solution metrics in CompressedSensingOptions format.
	  opts_list[k].resize( 1 );
	  // First copy the original options
	  opts_list[k][0] = opts;
	  // There is no need to update the needed metrics as they will not
	  // effect the solution produced if used again.
	}      
    }
  else
    {
      for ( int k = 0; k < num_rhs; k++ )
	{
	  // ( M x 1 ) rhs vector of Ax = b
	  RealVector b( Teuchos::View, B[k], M );
	  RealMatrix solution_metrics;
	  linearSolver_->solve( A, b, solutions[k], 
				solution_metrics );
	  // Store solution metrics in CompressedSensingOptions format
	  int num_solutions( solutions[k].numCols() );

	  opts_list[k].resize( num_solutions );
	  for ( int l = 0; l < num_solutions; l++ )
	    {
	      // First copy the original options
	      opts_list[k][l] = opts;
	      // Then update the needed metrics
	      opts_list[k][l].epsilon = solution_metrics(0,l);
	      opts_list[k][l].maxNumIterations = 
		solution_metrics(1,l);
	  
	    }
	  if ( !opts.storeHistory )
	    {
	      // Only return final solution 
	      int num_solutions( solutions[k].numCols() );
	      RealMatrix final_solution( Teuchos::Copy, solutions[k], 
					 N, 1, 0, num_solutions - 1 );
	      solutions[k].reshape( N, 1 );
	      solutions[k] = final_solution;
	      // Only return the options that produced the final solution
	      opts_list[k][0] = opts_list[k][num_solutions-1];
	      opts_list[k].resize( 1 );
	    }
    }
  }
};

void CompressedSensingTool::set_linear_solver( CompressedSensingOptions &opts )
{
  short solver( opts.solver );
  Real solver_tolerance( opts.solverTolerance );
  // override opts.standardizeInputs
  int standardize = true;
  //int standardize = false;
  if ( solver == SVD_LEAST_SQ_REGRESSION )
    standardize = false;

  switch( solver )
    {
      case SVD_LEAST_SQ_REGRESSION:
	{
	  LinearSolver_ptr lsq_solver( new LSQSolver() );
	  //linearSolver_->set_normalise_inputs( false );
	  //linearSolver_->set_solver_tolerance( solver_tolerance );
	  linearSolver_ = lsq_solver;
	  break;
	}
      case BASIS_PURSUIT:
	{
	  LinearSolver_ptr bp_solver( new BPSolver() );
	  linearSolver_ = bp_solver;
	  linearSolver_->set_normalise_inputs( standardize );
	  linearSolver_->set_solver_tolerance( solver_tolerance );
	  linearSolver_->set_conjugate_gradient_tolerance( 
				     opts.conjugateGradientsTolerance );
	  linearSolver_->set_verbosity( opts.verbosity+1 );
	  break;
	}
    case BASIS_PURSUIT_DENOISING:
      {
	LinearSolver_ptr bpdn_solver( new BPDNSolver() );
	linearSolver_ = bpdn_solver;
	RealVector residual_tols( 1 ); residual_tols[0] = opts.epsilon;
	linearSolver_->set_residual_tolerances( residual_tols );
	linearSolver_->set_normalise_inputs( standardize );
	linearSolver_->set_solver_tolerance( solver_tolerance );
	linearSolver_->set_conjugate_gradient_tolerance( 
				opts.conjugateGradientsTolerance );
	linearSolver_->set_verbosity( opts.verbosity+1 );
	break;
      }
    case ORTHOG_MATCH_PURSUIT:
      {
	LinearSolver_ptr omp_solver( new OMPSolver() );
	linearSolver_ = omp_solver;
	linearSolver_->set_normalise_inputs( standardize );
	linearSolver_->set_residual_tolerance( opts.epsilon );
	linearSolver_->set_solver_tolerance( solver_tolerance );
	linearSolver_->set_verbosity( opts.verbosity+1 );
	linearSolver_->set_max_iters( opts.maxNumIterations );	
	break;
      }
    case LASSO_REGRESSION:
      {
	LinearSolver_ptr lasso_solver( new LARSSolver() );
	linearSolver_ = lasso_solver;
	linearSolver_->set_normalise_inputs( standardize );
	boost::static_pointer_cast<LARSSolver>(linearSolver_)->set_sub_solver( LASSO_REGRESSION );
	linearSolver_->set_residual_tolerance( opts.epsilon );
	linearSolver_->set_solver_tolerance( solver_tolerance );
	linearSolver_->set_verbosity( opts.verbosity+1 );
	linearSolver_->set_max_iters( opts.maxNumIterations );
	boost::static_pointer_cast<LARSSolver>(linearSolver_)->set_delta( opts.delta );
	break;
      }
    case LEAST_ANGLE_REGRESSION:
      {
	LinearSolver_ptr lars_solver( new LARSSolver() );
	linearSolver_ = lars_solver;
	linearSolver_->set_normalise_inputs( standardize );
	boost::static_pointer_cast<LARSSolver>(linearSolver_)->set_sub_solver( LEAST_ANGLE_REGRESSION );
	linearSolver_->set_residual_tolerance( opts.epsilon );
	linearSolver_->set_solver_tolerance( solver_tolerance );
	linearSolver_->set_verbosity( opts.verbosity+1 );
	linearSolver_->set_max_iters( opts.maxNumIterations );
	boost::static_pointer_cast<LARSSolver>(linearSolver_)->set_delta( opts.delta );
	break;
      }
    case EQ_CON_LEAST_SQ_REGRESSION:
      {
	LinearSolver_ptr eqlsq_solver( new EqualityConstrainedLSQSolver() );
	linearSolver_ = eqlsq_solver;
	linearSolver_->set_num_primary_equations( opts.numFunctionSamples );
	linearSolver_->set_normalise_inputs( standardize );
	linearSolver_->set_verbosity( opts.verbosity+1 );
	break;
      }
    default:
      {
	std::string msg = "CompressedSensingTool::solve() ";
	msg += " incorrect solver specified";
	PCout << "solver = " << solver << std::endl;
	throw( std::runtime_error( msg ) );
	break;
      }
    }
}

LinearSolver_ptr CompressedSensingTool::get_linear_solver()
{
  return linearSolver_;
};


}  // namespace Pecos
