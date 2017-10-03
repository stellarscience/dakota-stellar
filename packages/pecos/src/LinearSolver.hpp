#ifndef LINEAR_SOLVERS_HPP
#define LINEAR_SOLVERS_HPP

#include "compressed_sensing.hpp"
#include <boost/shared_ptr.hpp>

namespace Pecos {

void normalise_columns( RealMatrix &A, RealVector &result );

class LinearSolver
{
protected:

  RealVector residualTols_;

  int maxIters_;

  int verbosity_;

  bool normaliseInputs_;

  Real solverTol_;
  
  Real conjugateGradientTol_;

  int numPrimaryEqs_;

public:

  LinearSolver() : 
    maxIters_(std::numeric_limits<int>::max()), verbosity_(0),
    normaliseInputs_(false), solverTol_(1.e-6), conjugateGradientTol_(-1.), 
    numPrimaryEqs_( -1 )
  { residualTols_.size( 1 ); };

  ~LinearSolver()
  {
    clear();
  };

  void clear()
  {
    maxIters_ = std::numeric_limits<int>::max();
    verbosity_ = 0; solverTol_ = 1.e-6; conjugateGradientTol_ = -1.;
    numPrimaryEqs_ = -1; residualTols_.size( 1 );
  }

  /**
   * \brief Find a regularized solution to AX = B
   */
  virtual void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
		      RealMatrix &result_1 )
  {
    std::string msg = "solve() Has not been implemented for ";
    msg += "this class.";
    throw( std::runtime_error( msg ) );
  };
  
  virtual void solve_using_points( RealMatrix &build_points, RealMatrix &B, 
				   RealMatrix &result_0, 
				   RealMatrix &result_1 )
  {
    std::string msg = "solve_using_points() Has not been implemented for ";
    msg += "this class.";
    throw( std::runtime_error( msg ) );
  };

  void set_verbosity( int verbosity )
  {
    verbosity_ = verbosity;
  };
  
  void set_residual_tolerance( Real tolerance )
  {
    RealVector residual_tols( 1, false );
    residual_tols[0] = tolerance;
    set_residual_tolerances( residual_tols );
  };

  void set_residual_tolerances( const RealVector &residual_tols )
  {
    residualTols_.sizeUninitialized( residual_tols.length() );
    residualTols_.assign( residual_tols );
  };


  void set_solver_tolerance( Real tolerance )
  {
    solverTol_  = tolerance;
  };

  void set_conjugate_gradient_tolerance( Real tolerance )
  {
    conjugateGradientTol_  = tolerance;
  };

  void set_num_primary_equations( int num_primary_eqs )
  {
    numPrimaryEqs_  = num_primary_eqs;
  };

  Real get_residual_tolerance() const
  {
    return residualTols_[0];
  };

  void get_residual_tolerances( RealVector &result_0 )
  {
    result_0.sizeUninitialized( residualTols_.length() );
    result_0.assign( residualTols_ );
  };

  void set_max_iters( int max_iters )
  {
    maxIters_ = max_iters;
  };

  void set_normalise_inputs( bool normalise )
  {
    normaliseInputs_ = normalise;
  };

  void copy( const LinearSolver &source )
  {
    set_residual_tolerances( source.residualTols_ );
    set_max_iters( source.maxIters_ );
    set_verbosity( source.verbosity_ );
  };

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

  void adjust_coefficients( const RealVector &normalisation_factors, 
			    RealMatrix &coefficients )
  {
    int num_coeff = coefficients.numRows(), num_qoi = coefficients.numCols();
    for ( int i = 0; i < num_qoi; i++ )
      {
	for ( int j = 0; j < num_coeff; j++ )
	  coefficients(j,i) /= normalisation_factors[j];
      } 
  };

  /**
   * For adaptive methods such as orthogonal least interpolation
   * it is sometimes necessary to reconstruct the Matrix A of Ax=b
   * for instance when performing cross validation
   */
  virtual void build_matrix( const RealMatrix &build_points,
			     RealMatrix &result_0 )
  {
    std::string msg = "linear_solver::build_matrix() Not implemented.";
    throw( std::runtime_error( msg ) );
  }
};

class BPSolver : public LinearSolver
{
protected:

public:

  BPSolver(){};

  ~BPSolver(){};

  /**
   * \brief Find the solution min ||x||_0 such that |AX = B||_2 == 0
   */
  void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
	      RealMatrix &result_1 )
  {
    if ( B.numCols() != 1 )
      throw( std::runtime_error(" BPSolver::solve() B must be a vector") );

    RealVector b( Teuchos::View, B[0], B.numRows() );
    RealMatrix A_copy( A );

    RealVector column_norms;
    if ( normaliseInputs_ )
      normalise_columns( A_copy, column_norms );

    BP_primal_dual_interior_point_method( A_copy, b, 
					  result_0, 
					  solverTol_, 
					  conjugateGradientTol_,
					  verbosity_ );

    if ( normaliseInputs_ )
      adjust_coefficients( column_norms, result_0 );

    result_1.shapeUninitialized( 2, 1 );
    result_1(0,0) = 0.;
    int num_non_zeros = 0;
    for ( int i = 0; i < result_0.numRows(); i++ )
      if ( std::abs( result_0(i,0) ) > std::numeric_limits<double>::epsilon() )
	num_non_zeros++;
    result_1(1,0) = num_non_zeros;
  };
};

class BPDNSolver : public LinearSolver
{ 
public:
  
  BPDNSolver(){};
  
  ~BPDNSolver(){};

  /**
   * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
   */
  void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
	      RealMatrix &result_1 )
  {
    if ( residualTols_.length() <= 0 )
      throw( std::runtime_error(" BPDNSolver::solve() set residual tols") );

    if ( B.numCols() != 1 )
      throw( std::runtime_error(" BPDNSolver::solve() B must be a vector") );

    RealVector b( Teuchos::View, B[0], B.numRows() );
    RealMatrix A_copy( A );

    RealVector column_norms;
    if ( normaliseInputs_ )
      normalise_columns( A_copy, column_norms );

    result_0.shapeUninitialized( A.numCols(), residualTols_.length() );
    result_1.shapeUninitialized( 2, residualTols_.length() );
    for ( int j = 0; j < residualTols_.length(); j++ )
      {
	RealMatrix x;
	BPDN_log_barrier_interior_point_method( A_copy, b, 
						x, 
						residualTols_[j],
						solverTol_, 
						conjugateGradientTol_,
						verbosity_ );
	
	if ( normaliseInputs_ )
	  adjust_coefficients( column_norms, x );
	
	result_1(0,j) = residualTols_[j];
	int num_non_zeros = 0;
	for ( int i = 0; i < result_0.numRows(); i++ )
	  {
	    if ( std::abs( x(i,0) ) > std::numeric_limits<double>::epsilon() )
	      num_non_zeros++;
	    result_0(i,j) = x(i,0);
	  }
	result_1(1,j) = num_non_zeros;
      }
  };
};

class OMPSolver : public LinearSolver
{
protected:
  IntVector ordering_; // enforce a set of columns to be chosen first

public:

  OMPSolver(){};

  ~OMPSolver(){};

  /**
   * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
   */
  void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
	      RealMatrix &result_1 )
  {
    if ( B.numCols() != 1 )
      throw( std::runtime_error(" OMPSolver::solve() B must be a vector") );

    RealVector b( Teuchos::View, B[0], B.numRows() );
    RealMatrix A_copy( A );

    RealVector column_norms;
    if ( normaliseInputs_ )
      normalise_columns( A_copy, column_norms );

    orthogonal_matching_pursuit( A_copy, b, result_0, result_1, 
				 residualTols_[0], maxIters_, verbosity_,
				 ordering_ );

    if ( normaliseInputs_ )
      adjust_coefficients( column_norms, result_0 );
  };

  void set_ordering( IntVector &ordering )
  {
    ordering_.resize( ordering.length() );
    ordering_.assign( ordering );
  };
};

class LARSSolver : public LinearSolver
{
private:

  int solver_;
  
  Real delta_;

public:
  LARSSolver() : solver_( LEAST_ANGLE_REGRESSION ), delta_( 0.0 ){};

  ~LARSSolver(){clear();};

  void clear()
  {
    LinearSolver::clear();
    solver_ = LEAST_ANGLE_REGRESSION; delta_ = 0.0;
  };

  void set_sub_solver( int solver_id )
  {
    if ( ( solver_id != LASSO_REGRESSION ) && 
	 ( solver_id != LEAST_ANGLE_REGRESSION ) )
      {
	std::stringstream msg;
	msg << "set_sub_solver() solver id must be either: " << LASSO_REGRESSION
	    << " or " << LEAST_ANGLE_REGRESSION << "\n";
	throw( std::runtime_error( msg.str() ) );
      }
    solver_ = solver_id;
  };

  void set_delta( Real delta )
  {
    delta_ = delta;
  }

  /**
   * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
   */
  void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
	      RealMatrix &result_1 )
  {
    if ( B.numCols() != 1 )
      throw( std::runtime_error(" LARSSolver::solve() B must be a vector") );

    RealVector b( Teuchos::View, B[0], B.numRows() );
    RealMatrix A_copy( A );
    
    RealVector column_norms;
    if ( normaliseInputs_ )
      normalise_columns( A_copy, column_norms );

    least_angle_regression( A_copy, b, result_0, result_1, 
			    residualTols_[0], solver_, delta_, 
			    maxIters_, verbosity_ );

    if ( normaliseInputs_ )
      adjust_coefficients( column_norms, result_0 );
  };
};

class COSAMPSolver : public LinearSolver
{
private:
  int sparsity_;

public:
  COSAMPSolver() : sparsity_( 0 ) {};

  ~COSAMPSolver(){clear();};

  void clear()
  {
    LinearSolver::clear();
    sparsity_ = 0;
  };

  void set_sparsity( int sparsity )
  {
    sparsity_ = sparsity;
  }

  /**
   * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
   */
  void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
	      RealMatrix &result_1 )
  {
    if ( B.numCols() != 1 )
      throw( std::runtime_error(" COSAMPSolver::solve() B must be a vector") );

    RealVector b( Teuchos::View, B[0], B.numRows() );
    RealMatrix A_copy( A );
    
    RealVector column_norms;
    if ( normaliseInputs_ )
      normalise_columns( A_copy, column_norms );

    cosamp( A_copy, b, result_0, result_1, sparsity_,
	    maxIters_, verbosity_ );

    if ( normaliseInputs_ )
      adjust_coefficients( column_norms, result_0 );
  };
};

class LSQSolver : public LinearSolver
{
public:
  LSQSolver(){};

  ~LSQSolver(){};

  /**
   * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
   */
  void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
	      RealMatrix &result_1 )
  {
    if ( B.numCols() != 1 )
      throw( std::runtime_error("LSQSolver::solve() B must be a vector") );

   if ( A.numRows() < A.numCols() )
     std::cout << "LSQSolver::solve() Warning A is under-determined. " <<
       "M = " << A.numRows() << " N = " << A.numCols() <<
       ". Returning minimum norm solution\n";
   

    RealVector b( Teuchos::View, B[0], B.numRows() );
    RealMatrix A_copy( A );
    
    RealVector column_norms;
    if ( normaliseInputs_ )
      normalise_columns( A_copy, column_norms );

    RealVector singular_values;
    int rank(0);
    svd_solve( A_copy, b, result_0, singular_values, rank, 
	       solverTol_ );
    
    result_1.shapeUninitialized( 2, 1 );
    RealVector residual( b );
    residual.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		       -1.0, A_copy, result_0, 1.0 );
    result_1(0,0) = residual.normFrobenius();
    int num_non_zeros = 0;
    for ( int i = 0; i < result_0.numRows(); i++ )
      {
	if ( std::abs( result_0(i,0) ) > std::numeric_limits<double>::epsilon() )
	  num_non_zeros++;
      }
    result_1(1,0) = num_non_zeros;	  

    if ( normaliseInputs_ )
      adjust_coefficients( column_norms, result_0 );
  };
};

class EqualityConstrainedLSQSolver : public LinearSolver
{

public:
  EqualityConstrainedLSQSolver(){};

  ~EqualityConstrainedLSQSolver(){};

  /**
   * \brief Find the solution min ||x||_0 such that |AX = B||_2 < eps
   */
  void solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, 
	      RealMatrix &result_1 )
  {
    if ( B.numCols() != 1 )
      throw( std::runtime_error(" EqualityConstrainedLSQSolver::solve() B must be a vector") );

    if ( numPrimaryEqs_ <= 0 )
      throw( std::runtime_error(" EqualityConstrainedLSQSolver::solve() set num primary equations") );

    if ( numPrimaryEqs_ > A.numCols() )
      throw( std::runtime_error(" EqualityConstrainedLSQSolver::solve() num primary equations is larger than the number of columns in A") );

    if ( A.numRows() < A.numCols() )
      throw( std::runtime_error(" EqualityConstrainedLSQSolver::solve() A is underdetermined") );

    RealMatrix A_copy( A );
    
    RealVector column_norms;
    if ( normaliseInputs_ )
      normalise_columns( A_copy, column_norms );

    RealMatrix C_eq( Teuchos::View, A_copy, numPrimaryEqs_,
		     A_copy.numCols(), 0, 0 );
    RealMatrix A_eq( Teuchos::View, A_copy, 
		     A_copy.numRows() - numPrimaryEqs_,
		     A_copy.numCols(), numPrimaryEqs_, 0 );
    RealVector d_eq( Teuchos::View, B.values(), numPrimaryEqs_ );
    RealVector b_eq( Teuchos::View, 
		     B.values() + numPrimaryEqs_, 
		     B.numRows() - numPrimaryEqs_ );
    equality_constrained_least_squares_solve( A_eq, b_eq, C_eq, d_eq,
					      result_0 );
    
    
    result_1.shapeUninitialized( 2, 1 );
    RealVector residual( b_eq );
    residual.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		       -1.0, A_eq, result_0, 1.0 );
    result_1(0,0) = residual.normFrobenius();
    int num_non_zeros = 0;
    for ( int i = 0; i < result_0.numRows(); i++ )
      {
	if ( std::abs( result_0(i,0) ) > std::numeric_limits<double>::epsilon() )
	  num_non_zeros++;
      }
    result_1(1,0) = num_non_zeros;	  

    if ( normaliseInputs_ )
      adjust_coefficients( column_norms, result_0 );
  };
};

typedef boost::shared_ptr<LinearSolver> LinearSolver_ptr;

/**
 * \brief Specify a set of options for using the CompressedSensingTool
 */
class CompressedSensingOptions
{
public:
  //solverType solver; //!< Specify which regression solver to use. See solverType
  short solver; //!< Specify which regression solver to use: see pecos_global_defs
  Real solverTolerance; //!< Specify the internal tolerance of the solver
  Real epsilon;         //!< Specify the residual tolerance of the solver
  Real delta;           //!< Specify the regularization parameter value
  int maxNumIterations; //!< Specify the maximum number of solver iterations
  bool standardizeInputs;  //!< Specify if the inputs need to be standardized
  bool storeHistory;       //<! Specify if the solution history should be stored 
  Real conjugateGradientsTolerance; //<! Specify wether to use conjugate gradients internally to solve newton step in BP and BPDN.  If < 0 cholesky factorization will be used.
  int verbosity;           //!< The verbosity level. 0: off, 1: warnings on,  2: all print statements on.
  int numFunctionSamples; //!< The number of function samples used to construct A and B. Used when A contains gradient information. If zero then numFunctionSamples = A.numRows()

public:

  CompressedSensingOptions() : 
    solver( DEFAULT_LEAST_SQ_REGRESSION ), solverTolerance( -1. ),
    epsilon( 0.0 ), delta( 0.0 ),
    maxNumIterations( std::numeric_limits<int>::max() ), 
    standardizeInputs( false ), storeHistory( false ), 
    conjugateGradientsTolerance( -1 ), verbosity( 0 ), numFunctionSamples( 0 )
  {};

  ~CompressedSensingOptions(){};

  void print()
  {
    std::cout << "Solver: " << solver << "\n";
    std::cout << "Solver Tolerance: " << solverTolerance << "\n";
    std::cout << "Epsilon: " << epsilon << "\n";
    std::cout << "Delta: " << delta << "\n";
    std::cout << "MaxNumIterations: " << maxNumIterations << "\n";
    std::cout << "StandardizeInputs: " << standardizeInputs << "\n";
    std::cout << "StoreHistory: " << storeHistory << "\n";
    std::cout << "Verbosity: " << verbosity << "\n";
  };

  CompressedSensingOptions& operator=(const CompressedSensingOptions& source )
  {
    if(this == &source)
      return (*this);
    solver = source.solver;
    solverTolerance = source.solverTolerance; 
    epsilon = source.epsilon;         
    delta = source.delta;
    maxNumIterations = source.maxNumIterations; 
    standardizeInputs = source.standardizeInputs;
    storeHistory = source.storeHistory;     
    conjugateGradientsTolerance = source.conjugateGradientsTolerance; 
    verbosity = source.verbosity;
    numFunctionSamples = source.numFunctionSamples;
    return (*this);
  }

  CompressedSensingOptions( const CompressedSensingOptions& source )
  {
    solver = source.solver;
    solverTolerance = source.solverTolerance; 
    epsilon = source.epsilon;         
    delta = source.delta;
    maxNumIterations = source.maxNumIterations; 
    standardizeInputs = source.standardizeInputs;
    storeHistory = source.storeHistory;     
    conjugateGradientsTolerance = source.conjugateGradientsTolerance; 
    verbosity = source.verbosity;
    numFunctionSamples = source.numFunctionSamples;
  }
};

typedef std::vector< std::vector<CompressedSensingOptions> > CompressedSensingOptionsList;

/**
 * \class CompressedSensingTool
 * \brief Tool that implements a number of popular compressed sensing algorithms 
 */
class CompressedSensingTool
{
protected:
  LinearSolver_ptr linearSolver_;

public:

  /// Default constructor
  CompressedSensingTool()
  {   
    std::cout.precision( std::numeric_limits<Real>::digits10 );
    std::cout.setf( std::ios::scientific );
  };

  // Deconstructor
  ~CompressedSensingTool(){};

  //! @name Interface tools.
  //@{ 

  /**
   * \brief Wrapper to call any of the compressed sensing methods
   *
   * \param A ( M x N ) matrix of the linear system AX=B
   *
   * \param B ( M x num_rhs ) matrix of the linear system AX=B
   *
   * \param solutions (output) vector containing multiple solutions
   * to AX=B. Each entry in solutions is a ( N x num_rhs ) matrix
   * opts.solver=LS will return only one solution whilst methods such
   * as OMP, LARS, and LASSO will return a history of solutions
   *
   * \param opts specifies the method options
   *
   * \param opts_list specifies the method options that can be used to
   * reproduce all the solutions found.
   */
  void solve( RealMatrix &A, 
	      RealMatrix &B, 
	      RealMatrixArray &solutions,
	      CompressedSensingOptions &opts,
	      CompressedSensingOptionsList &opts_list );

  void set_linear_solver( CompressedSensingOptions &opts );

  LinearSolver_ptr get_linear_solver();
};

} // namespace Pecos

#endif //LINEAR_SOLVERS_HPP
