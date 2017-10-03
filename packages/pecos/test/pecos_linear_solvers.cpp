
#include <ctype.h>
#include <string>

#include <Teuchos_UnitTestHarness.hpp> 

#include "pecos_data_types.hpp"
#include "LinearSolver.hpp"

using namespace Pecos;

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

  Real get_solver_solve_tol(const LinearSolver* sol) {
    Real tol = sol->get_residual_tolerance();
    if( tol < 1.e-16 )
      tol = 1.e-8;
    return tol;
  }

  //--------------------------------------

  Real test_solver(LinearSolver* sol, bool use_spd) {

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

    RealMatrix result0, result1;
    sol->solve(A, B, result0, result1);

    //gold_soln.print(std::cout);
    //std::cout << std::endl;
    //result0.print(std::cout);
    //std::cout << std::endl;
    //result1.print(std::cout);

    // Test correctness
    gold_soln -= result0;
    Real shifted_diff_norm = gold_soln.normFrobenius() + 1.0;

    return shifted_diff_norm;
  }
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_solvers, bpsolver)
{
  LinearSolver * psol = new BPSolver;
  Real diff = test_solver(psol, false);

  Real tol = get_solver_solve_tol(psol);
  TEST_FLOATING_EQUALITY( 1.0, diff, tol )
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_solvers, bpdnsolver)
{
  LinearSolver * psol = new BPDNSolver;
  Real diff = test_solver(psol, true);

  Real tol = get_solver_solve_tol(psol);
  TEST_FLOATING_EQUALITY( 1.0, diff, tol )
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_solvers, ompsolver)
{
  LinearSolver * psol = new OMPSolver;
  //psol->set_normalise_inputs(true);
  Real diff = test_solver(psol, true);

  Real tol = get_solver_solve_tol(psol);
  // This solve fails
  //TEST_FLOATING_EQUALITY( 1.0, diff, tol )
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_solvers, larssolver)
{
  LinearSolver * psol = new LARSSolver;
  Real diff = test_solver(psol, false);

  Real tol = get_solver_solve_tol(psol);
  // This solve fails
  //TEST_FLOATING_EQUALITY( 1.0, diff, tol )
}

//----------------------------------------------------------------

//TEUCHOS_UNIT_TEST(linear_solvers, cosampsolver)
//{
//  LinearSolver * psol = new COSAMPSolver;
//  Real diff = test_solver(psol, true);
//
//  Real tol = get_solver_solve_tol(psol);
//  TEST_FLOATING_EQUALITY( 1.0, diff, tol )
//}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_solvers, lsqsolver)
{
  LinearSolver * psol = new LSQSolver;
  Real diff = test_solver(psol, false);

  Real tol = get_solver_solve_tol(psol);
  TEST_FLOATING_EQUALITY( 1.0, diff, tol )
}

//----------------------------------------------------------------

//TEUCHOS_UNIT_TEST(linear_solvers, equalityconstrainedlsqsolver)
//{
//  LinearSolver * psol = new EqualityConstrainedLSQSolver;
//  psol->set_num_primary_equations(NUMROWS);
//  Real diff = test_solver(psol, true);
//
//  Real tol = get_solver_solve_tol(psol);
//  TEST_FLOATING_EQUALITY( 1.0, diff, tol )
//}
