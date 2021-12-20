/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <ctype.h>
#include <string>

#include <Teuchos_UnitTestHarness.hpp> 
#include <Teuchos_SerialDenseHelpers.hpp> 

#include "linear_algebra.hpp"

using namespace Pecos;
using namespace Pecos::util;

namespace {

  const int NUMROWS = 5;
  const int NUMCOLS = 5;
  //--------------------------------------
  // Generate a random matrix 
  //   NUMROWS x NUMCOLS
  //   ... using a specified random seed
  RealMatrix get_test_matrix(unsigned int seed) {
    RealMatrix mat(NUMROWS,NUMCOLS);
    Teuchos::ScalarTraits<Real>::seedrandom(seed);
    mat.random();
    return mat;
  }
  //--------------------------------------
  // Generate a random matrix 
  //   NUMROWS x NUMCOLS
  //   we may want to dial in scaling 
  RealMatrix get_test_matrix() {
    RealMatrix mat(NUMROWS,NUMCOLS);
    mat.random();
    return mat;
  }
  //--------------------------------------
  // A random matrix -- positive definite
  //   NUMROWS x NUMCOLS
  //   (a real matrix A is positive definite *iff* its 
  //   symmetric part (A + A^T)/2 is positive definite. 
  //   we may want to dial in scaling 
  
  RealMatrix get_spd_test_matrix() {
    RealMatrix mat = get_test_matrix();
    RealMatrix trans_mat(mat, Teuchos::TRANS);
    mat += trans_mat;
    mat.scale(0.5);
    Real shift = 2.0;
    
    // *** *** *** *** *** *** *** *** ***
    // making larger shifts as matrix size grows
    // mat += NUMROWS*eye(NUMROWS) 
    for(int i=0; i<mat.numRows(); ++i)
      //mat(i,i) += shift;
      mat(i,i) += mat.numRows();
    return mat;
  }
  //--------------------------------------
  // A random matrix -- not positive definite
  // Same as above, but with negative diagonal shifting
  // or we can enforce non positive definite 
  // by composing A with negative eigenvalues 
  // D = diag(-abs(rand(NUMROWS)))
  // V = orth(rand(NUMROWS))
  // A = V*D*V'
  
  RealMatrix get_nspd_test_matrix() {
    RealMatrix mat = get_test_matrix();
    RealMatrix trans_mat(mat, Teuchos::TRANS);
    mat += trans_mat; 
    mat.scale(0.5);
    //for(int i=0; i<mat.numRows(); ++i)
    //	mat(i,i) -= mat.numRows();
    mat(1,1) -= mat.numRows();
    return mat; 
    }

    //-------------------------------------
    // A random matrix -- singular

    RealMatrix get_singular_test_matrix() {
      RealMatrix mat = get_test_matrix();
      for(int i=0; i<mat.numCols(); ++i)
        mat(2,i) = 2.0*mat(1,i);
      return mat;
    }  

    //-------------------------------------
    // A random matrix -- Overdetermined

    RealMatrix get_overd_test_matrix() {
      RealMatrix mat(NUMROWS,NUMCOLS+3);
      mat.random();
      return mat;
    }
 
    //-------------------------------------
    // A random matrix -- Underdetermined

    RealMatrix get_underd_test_matrix() {
      RealMatrix mat(NUMROWS+3,NUMCOLS);
      mat.random();
      return mat;
    }  
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, cholesky)
{

  // **********************************************************
  // ***** Test Cholesky with a positive definite matrix  *****
  // Construct a rnd system matrix
  RealMatrix A = get_spd_test_matrix();
  // Random solution vector
  RealMatrix gold_soln(NUMROWS,1);
  gold_soln.random();
  
  // Create the RHS needed to recover the gold soln via solve
  RealMatrix B(gold_soln);
  B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);

  RealMatrix L, result;
  
  // Colesky to get L then solve using the L
  int info = Pecos::util::cholesky( A, L, Teuchos::LOWER_TRI, false );
  TEST_EQUALITY( 0, info )
  info =  Pecos::util::solve_using_cholesky_factor( L, B, result, Teuchos::LOWER_TRI);
  TEST_EQUALITY( 0, info )

  // Test correctness
  gold_soln -= result;
  Real shifted_diff_norm = gold_soln.normFrobenius() + 1.0;
  Real tol = 1.0e-8; // get_solver_solve_tol(psol);
  TEST_FLOATING_EQUALITY( 1.0, shifted_diff_norm, tol )
  
  // *********************************************************
  // *** Test Cholesky with a not positive definite matrix ***
  A = get_nspd_test_matrix();
  B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);
  info = Pecos::util::cholesky( A, L, Teuchos::LOWER_TRI, false );
  TEST_INEQUALITY( 0, info ) //-- info=2

  // *********************************************************
  // *** Test Cholesky with a singular matrix ***
  A = get_singular_test_matrix();
  B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);
  info = Pecos::util::cholesky( A, L, Teuchos::LOWER_TRI, false );
  TEST_INEQUALITY( 0, info )  //-- info=1
  
  // DEBUG
  //A = get_rect_test_matrix();
  //A.print(std::cout);
  //std::cout << "\n\n";

}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, teuchos_cholesky)
{
  // ===> computes L internally
  
  // ***** Test Cholesky with a positive definite matrix  *****
  RealMatrix A = get_spd_test_matrix();
  RealMatrix gold_soln(NUMROWS,1);
  gold_soln.random();

  // Create the RHS needed to recover the gold soln via solve
  RealMatrix B(gold_soln);
  B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);

  RealMatrix result;
  Real rcond = 1.0;
  int info = Pecos::util::teuchos_cholesky_solve( A, B, result, rcond );
  TEST_EQUALITY( 0, info )

  // Test correctness
  gold_soln -= result;
  Real shifted_diff_norm = gold_soln.normFrobenius() + 1.0;

  Real tol = 1.0e-8; // get_solver_solve_tol(psol);
  TEST_FLOATING_EQUALITY( 1.0, shifted_diff_norm, tol )

}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, qr_solve)
{
  // Might need to try an over (under) determined system for better testing
  RealMatrix A = get_test_matrix();
  RealMatrix gold_soln(NUMROWS,1);
  gold_soln.random();

  // Create the RHS needed to recover the gold soln via solve
  RealMatrix B(gold_soln);
  B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);

  RealMatrix result;
  Pecos::util::qr_solve(A, B, result);

  // Test correctness
  gold_soln -= result;
  Real shifted_diff_norm = gold_soln.normFrobenius() + 1.0;

  Real tol = 1.0e-8; // get_solver_solve_tol(psol);
  TEST_FLOATING_EQUALITY( 1.0, shifted_diff_norm, tol )

  // **************************************** 
  // *** Test random rectangular matrices ***
  // LAPACK's DGELS would attempt to find the minimum norm solution
  A = get_overd_test_matrix();
  B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);
  Pecos::util::qr_solve(A, B, result);
  gold_soln -= result;
  shifted_diff_norm = gold_soln.normFrobenius() + 1.0;
  TEST_FLOATING_EQUALITY( 1.0, shifted_diff_norm, tol )

}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, svd)
{
  // Might need to try a rank-deficient system for better testing
  RealMatrix A = get_test_matrix();

  RealMatrix gold_soln(NUMROWS,1);
  gold_soln.random();

  RealMatrix B(gold_soln);
  // Create the RHS needed to recover the gold soln via solve
  B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);

  int rank;
  RealMatrix result0;
  RealVector result1;
  Pecos::util::svd_solve(A, B, result0, result1, rank);

  // Could add check of singular values
  //A.print(std::cout);
  //std::cout << "\n\n";
  //result1.print(std::cout);

  // Test correctness
  gold_soln -= result0;
  Real shifted_diff_norm = gold_soln.normFrobenius() + 1.0;

  Real tol = 1.0e-8; // get_solver_solve_tol(psol);
  TEST_FLOATING_EQUALITY( 1.0, shifted_diff_norm, tol )

  // **************************************** 
  // *** Test a random rectangular matrix ***
  // Not really a corner case, every real matrix has an SVD with real entries
  A = get_overd_test_matrix();
  B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);
  Pecos::util::svd_solve(A, B, result0, result1, rank);
  gold_soln -= result0;
  shifted_diff_norm = gold_soln.normFrobenius() + 1.0;
  TEST_FLOATING_EQUALITY( 1.0, shifted_diff_norm, tol )

}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, conj_grad_solv)
{
  // Might need to try a rank-deficient system for better testing
  RealMatrix A = get_spd_test_matrix();

  RealMatrix gold_soln(NUMROWS,1);
  gold_soln.random();

  RealVector b(NUMROWS);
  // Create the RHS needed to recover the gold soln via solve
  b.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);

  int iters;
  Real r_norm;
  RealVector result;
  Pecos::util::conjugate_gradients_solve(A, b, result, r_norm, iters /* use defaults for tol, max_iter, verbosity */);

  // Test correctness
  TEST_EQUALITY( 5, iters )
  gold_soln -= result;
  Real shifted_diff_norm = gold_soln.normFrobenius() + 1.0;

  Real tol = 1.0e-14; // CG for a spd matrix should nail the solution
  TEST_ASSERT( r_norm < tol )
  TEST_FLOATING_EQUALITY( 1.0, shifted_diff_norm, tol )
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, lu_solv)
{
  // Might need to try a rank-deficient system for better testing
  RealMatrix A = get_test_matrix();

  RealMatrix gold_soln(NUMROWS,1);
  gold_soln.random();

  RealMatrix B(gold_soln);
  // Create the RHS needed to recover the gold soln via solve
  B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, gold_soln, 0.0);

  RealVector result1;
  RealVector result2;
  Pecos::util::lu_solve(A, B, result1, true /* copy A */, Teuchos::NO_TRANS);
  Pecos::util::lu_solve(A, B, result2, false, Teuchos::NO_TRANS);

  // Test correctness
  Real tol = 1.0e-14;
  for( int i=0; i<result1.length(); ++i )
    TEST_FLOATING_EQUALITY( result1(i), result2(i), tol )

  gold_soln -= result1;
  Real shifted_diff_norm = gold_soln.normFrobenius() + 1.0;
  TEST_FLOATING_EQUALITY( 1.0, shifted_diff_norm, tol )
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, substitution_solve)
{
  RealMatrix A = get_test_matrix();

  RealMatrix gold_soln(NUMROWS,1);
  gold_soln.random();

  RealMatrix B(gold_soln);
  // Create the RHS needed to recover the gold soln via solve
  for( int i=0; i<A.numRows(); ++i ) {
    B(i,0) = 1.0*gold_soln(i,0);
    for( int j=0; j<i; ++j )
      B(i,0) += A(i,j)*gold_soln(j,0);
  }

  RealMatrix result;
  Pecos::util::substitution_solve( A, B, result, Teuchos::NO_TRANS, Teuchos::LOWER_TRI, Teuchos::UNIT_DIAG );

  // Test correctness
  RealMatrix test_gold(gold_soln);
  test_gold -= result;
  Real shifted_diff_norm = test_gold.normFrobenius() + 1.0;

  Real tol = 1.0e-8; // get_solver_solve_tol(psol);
  TEST_FLOATING_EQUALITY( 1.0, shifted_diff_norm, tol )


  // Now test another variant using NON_UNIT_DIAG

  // Create the RHS needed to recover the gold soln via solve
  for( int i=0; i<A.numRows(); ++i ) {
    B(i,0) = A(i,i)*gold_soln(i,0);
    for( int j=0; j<i; ++j )
      B(i,0) += A(i,j)*gold_soln(j,0);
  }

  Pecos::util::substitution_solve( A, B, result, Teuchos::NO_TRANS, Teuchos::LOWER_TRI, Teuchos::NON_UNIT_DIAG );

  // Test correctness
  gold_soln -= result;
  shifted_diff_norm = gold_soln.normFrobenius() + 1.0;

  TEST_FLOATING_EQUALITY( 1.0, shifted_diff_norm, tol )
}

//----------------------------------------------------------------

namespace
{ 
  RealMatrix get_lu_test_matrix() {
    // Fix the seed so we get the correct matrix used in Matlab
    return get_test_matrix(31415);
  }

  void get_lu_gold_objects(RealMatrix & gold_L, RealMatrix & gold_U, IntVector & gold_pivots) {

    RealMatrix A = get_lu_test_matrix();

    // Set the gold test results based on Matlab output
    gold_L.shape(A.numRows(), A.numCols());
    gold_L(0,0) =  1.000000000000000;
    gold_L(1,0) =  0.688380967523950; gold_L(1,1) =  1.000000000000000;
    gold_L(2,0) = -0.085980045315419; gold_L(2,1) = -0.324496376307591; gold_L(2,2) =  1.000000000000000;
    gold_L(3,0) =  0.643071510911476; gold_L(3,1) =  0.186084843028333; gold_L(3,2) =  0.751463750262584; gold_L(3,3) =  1.000000000000000;
    gold_L(4,0) = -0.043127161426243; gold_L(4,1) =  0.699569598832772; gold_L(4,2) =  0.150356737437066; gold_L(4,3) =  0.237219736177574; gold_L(4,4) =  1.000000000000000;

    // Set the gold test results based on Matlab output
    gold_U.shape(A.numRows(), A.numCols());
    gold_U(0,0) =  0.503140000000000; gold_U(0,1) =  0.885923000000000; gold_U(0,2) = -0.371340000000000; gold_U(0,3) = -0.427508000000000; gold_U(0,4) = -0.368267000000000;
                                      gold_U(1,1) = -1.340183531891720; gold_U(1,2) =  0.596169388480343; gold_U(1,3) =  1.184995370664229; gold_U(1,4) =  1.179052993767142;
                                                                        gold_U(2,2) =  1.058995976199956; gold_U(2,3) =  1.118633546509109; gold_U(2,4) =  0.804724810603881;
                                                                                                          gold_U(3,3) = -0.342792022081881; gold_U(3,4) =  0.320098600735852;
                                                                                                                                            gold_U(4,4) = -0.988404242883647;
    gold_pivots.size(A.numRows());
    gold_pivots(0) = 1;
    gold_pivots(1) = 3;
    gold_pivots(2) = 0;
    gold_pivots(3) = 2;
    gold_pivots(4) = 4;
  }
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, piv_lu_factor)
{
  const RealMatrix A = get_lu_test_matrix();

  RealMatrix L;
  RealMatrix U;
  IntVector pivots;
  int max_iters = A.numRows();
  int num_initial_rows = 0;

  Pecos::util::pivoted_lu_factorization( A, L, U, pivots);

  // Test correctness
  RealMatrix gold_L;
  RealMatrix gold_U;
  IntVector gold_pivots;
  get_lu_gold_objects(gold_L, gold_U, gold_pivots);

  const Real tol = 2.e-5; // This tolerance was determined empirically via comparison with Matlab results
                          // Is this sufficient? - RWH
  for( int i=0; i<NUMROWS; ++i )
    for( int j=0; j<NUMCOLS; ++j )
      TEST_FLOATING_EQUALITY( L(i,j), gold_L(i,j), tol )

  for( int i=0; i<NUMROWS; ++i )
    for( int j=0; j<NUMCOLS; ++j )
      TEST_FLOATING_EQUALITY( U(i,j), gold_U(i,j), tol )

  for( int i=0; i<NUMROWS; ++i )
    TEST_EQUALITY( pivots(i), gold_pivots(i) )
}


//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, trunc_lu_piv_factor)
{
  const RealMatrix A = get_lu_test_matrix();

  RealMatrix L;
  RealMatrix U;
  IntVector pivots;
  int max_iters = A.numRows();
  int num_initial_rows = 0;

  Pecos::util::truncated_pivoted_lu_factorization( A, L, U, pivots, max_iters, num_initial_rows);

  // Test correctness
  RealMatrix gold_L;
  RealMatrix gold_U;
  IntVector gold_pivots;
  get_lu_gold_objects(gold_L, gold_U, gold_pivots);

  const Real tol = 2.e-5; // This tolerance was determined empirically via comparison with Matlab results
                          // Is this sufficient? - RWH
  for( int i=0; i<NUMROWS; ++i )
    for( int j=0; j<NUMCOLS; ++j )
      TEST_FLOATING_EQUALITY( L(i,j), gold_L(i,j), tol )

  for( int i=0; i<NUMROWS; ++i )
    for( int j=0; j<NUMCOLS; ++j )
      TEST_FLOATING_EQUALITY( U(i,j), gold_U(i,j), tol )

  for( int i=0; i<NUMROWS; ++i )
    TEST_EQUALITY( pivots(i), gold_pivots(i) )
}


//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, compl_piv_lu_factor)
{
  const RealMatrix A = get_lu_test_matrix();

  RealMatrix L;
  RealMatrix U;
  IntVector row_pivots;
  IntVector col_pivots;
  int max_iters = A.numRows();

  Pecos::util::complete_pivoted_lu_factorization( A, L, U, row_pivots, col_pivots, max_iters);

  RealMatrix chkA(L);
  chkA.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, L, U, 0.0);

  RealMatrix permutedA = get_lu_test_matrix();
  Pecos::util::permute_matrix_rows( permutedA, row_pivots);
  Pecos::util::permute_matrix_columns( permutedA, col_pivots);

  const Real tol = 2.e-16; // This tolerance was determined empirically via comparison with Matlab results
                           // Is this sufficient? - RWH
  for( int i=0; i<NUMROWS; ++i )
    for( int j=0; j<NUMCOLS; ++j )
      TEST_FLOATING_EQUALITY( chkA(i,j), permutedA(i,j), tol )
}

//----------------------------------------------------------------

namespace
{ 
  RealMatrix get_qr_test_matrix() {
    // Fix the seed so we get the correct matrix used in Matlab
    return get_test_matrix(31415);
  }

  void get_qr_gold_objects(RealMatrix & gold_Q, RealMatrix & gold_R ) {

    RealMatrix A = get_qr_test_matrix();

    // Set the gold test results based on Matlab output
    gold_Q.shape(A.numRows(), A.numCols());
    gold_Q(0,0) = -0.062431348529951; gold_Q(0,1) = -0.253837855859339; gold_Q(0,2) = -0.796816908284737; gold_Q(0,3) = -0.544370376124644; gold_Q(0,4) =  0.020307090459546;
    gold_Q(1,0) =  0.726114394344877; gold_Q(1,1) = -0.374486261339218; gold_Q(1,2) =  0.237587073830810; gold_Q(1,3) = -0.238963092071478; gold_Q(1,4) =  0.467939044502058;
    gold_Q(2,0) =  0.466943480665931; gold_Q(2,1) = -0.076791874049033; gold_Q(2,2) = -0.498515156520094; gold_Q(2,3) =  0.705513228319441; gold_Q(2,4) = -0.172628263414914;
    gold_Q(3,0) =  0.499843329312193; gold_Q(3,1) =  0.623688121865361; gold_Q(3,2) =  0.005682878696861; gold_Q(3,3) = -0.374013166908812; gold_Q(3,4) = -0.470373924292917;
    gold_Q(4,0) = -0.031315252698830; gold_Q(4,1) =  0.632805276309391; gold_Q(4,2) = -0.245123521167450; gold_Q(4,3) =  0.094460414305004; gold_Q(4,4) =  0.727714591528301;

    // Set the gold test results based on Matlab output
    gold_R.shape(A.numRows(), A.numCols());
    gold_R(0,0) =  0.692921120857201; gold_R(0,1) =  0.435964770875061; gold_R(0,2) =  0.137894805505129; gold_R(0,3) =  0.264458631508166; gold_R(0,4) =  0.589056375799687;
                                      gold_R(1,1) = -1.520383424492604; gold_R(1,2) =  0.447165579176710; gold_R(1,3) =  1.077124947661094; gold_R(1,4) =  0.561449734000936;
                                                                        gold_R(2,2) = -1.279573136431154; gold_R(2,3) = -1.160812817644980; gold_R(2,4) = -0.908246109658568;
                                                                                                          gold_R(3,3) = -0.249525546772342; gold_R(3,4) =  0.139641452194404;
                                                                                                                                            gold_R(4,4) = -0.719276189874912;
  }
}
TEUCHOS_UNIT_TEST(linear_algebra, qr_factor)
{
  const RealMatrix A = get_lu_test_matrix();

  RealMatrix Q;
  RealMatrix R;

  Pecos::util::qr_factorization( A, Q, R );

  // Test correctness
  RealMatrix gold_Q;
  RealMatrix gold_R;
  get_qr_gold_objects(gold_Q, gold_R);
  // Why are there sign differences in rows/cols?  - RWH
  const Real tol = 5.2e-5; // This tolerance was determined empirically via comparison with Matlab results
                          // Is this sufficient? - RWH
  for( int i=0; i<NUMROWS; ++i )
    for( int j=0; j<NUMCOLS; ++j )
      TEST_FLOATING_EQUALITY( std::fabs(Q(i,j)), std::fabs(gold_Q(i,j)), tol )

  for( int i=0; i<NUMROWS; ++i )
    for( int j=0; j<NUMCOLS; ++j )
      TEST_FLOATING_EQUALITY( std::fabs(R(i,j)), std::fabs(gold_R(i,j)), tol )
}

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(linear_algebra, permute_rows)
{
  RealMatrix L;
  RealMatrix U;
  IntVector pivots;
  get_lu_gold_objects(L, U, pivots);

  RealMatrix chkA(L);
  chkA.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, L, U, 0.0);

  RealMatrix permutedA = get_lu_test_matrix();
  Pecos::util::permute_matrix_rows( permutedA, pivots);

  const Real tol = 2.e-6; // This tolerance was determined empirically via comparison with Matlab results
                           // Is this sufficient? - RWH
  for( int i=0; i<NUMROWS; ++i )
    for( int j=0; j<NUMCOLS; ++j )
      TEST_FLOATING_EQUALITY( chkA(i,j), permutedA(i,j), tol )
}
