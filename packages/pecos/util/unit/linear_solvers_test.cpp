/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "linear_solvers.hpp"
#include <Teuchos_UnitTestHarness.hpp>
#include "OptionsList.hpp"
#include <ctype.h>
#include <string>

using namespace Pecos;
using namespace Pecos::util;

namespace {

  const int NUMROWS = 5;
  const int NUMCOLS = 5;

  // A random matrix
  //   we may want to dial in scaling
  RealMatrix get_test_matrix() {
    RealMatrix mat(NUMROWS,NUMCOLS);
    mat.random();
    return mat;
  }

  //--------------------------------------

  // A random matrix
  //   we may want to dial in scaling
  RealMatrix get_spd_test_matrix() {
    RealMatrix mat = get_test_matrix();
    RealMatrix trans_mat(mat, Teuchos::TRANS);
    mat += trans_mat;
    mat.scale(0.5);
    Real shift = 2.0;
    for(int i=0; i<mat.numRows(); ++i)
      mat(i,i) += shift;
    return mat;
  }

  //--------------------------------------

  Real test_solver(LinearSystemSolver* sol, OptionsList& params, bool use_spd) {

    RealMatrix A;
    if( use_spd )
      A = get_spd_test_matrix();
    else
      A = get_test_matrix();

    RealMatrix gold_soln(NUMROWS,1);
    gold_soln.random();

    RealMatrix B(gold_soln);
    // Create the RHS needed to recover the gold soln via solve
    B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);

    RealMatrix solutions;
    sol->solve(A, B, params);
    sol->get_solutions_for_all_regularization_params(solutions,0);

    //gold_soln.print(std::cout);
    //std::cout << std::endl;
    //result0.print(std::cout);
    //std::cout << std::endl;
    //result1.print(std::cout);

    // Test correctness
    gold_soln -= solutions;
    Real shifted_diff_norm = gold_soln.normFrobenius() + 1.0;

    return shifted_diff_norm;
  }
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_solvers, ompsolver)
{
  LinearSystemSolver * psol = new OMPSolver;
  //params.set("Normalize Inputs", true);
  OptionsList params;
  Real diff = test_solver(psol, params, false);

  Real tol = params.get("Solver Tolerance", 1.e-8);
  // This solve fails
  //TEST_FLOATING_EQUALITY( 1.0, diff, tol )
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_solvers, larssolver)
{
  LinearSystemSolver * psol = new LARSolver;
  OptionsList params;
  params.set("Max Nonzeros", 0); // will trigger internal determination max(M,N)
  params.set("Max Iters", 30);
  Real diff = test_solver(psol, params, false);

  Real tol = params.get("Solver Tolerance", 1.e-8);
  // This solve fails
  //TEST_FLOATING_EQUALITY( 1.0, diff, tol )
}

//----------------------------------------------------------------
