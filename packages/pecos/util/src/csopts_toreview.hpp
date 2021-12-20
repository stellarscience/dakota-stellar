/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#error "csopts_toreview.hpp is an inactive header"

#ifndef COMPRESSED_SENSING_HPP
#define COMPRESSED_SENSING_HPP

#include "linear_algebra.hpp"
#include "linear_solvers.hpp"
//#include "pecos_global_defs.hpp"

namespace Pecos {

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

typedef std::shared_ptr<LinearSolver> LinearSolver_ptr;


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

#endif //COMPRESSED_SENSING_HPP
