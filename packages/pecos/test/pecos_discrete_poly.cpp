/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */


#include <ctype.h>
#include <string>

#include <Teuchos_UnitTestHarness.hpp> 

#include "pecos_data_types.hpp"
#include "BasisPolynomial.hpp"
#include "CharlierOrthogPolynomial.hpp"
#include "KrawtchoukOrthogPolynomial.hpp"
#include "MeixnerOrthogPolynomial.hpp"
#include "HahnOrthogPolynomial.hpp"
#include "NumericGenOrthogPolynomial.hpp"

using namespace Pecos;

namespace {

  //------------------------------------
  // Compute known exact orthogonality value
  //------------------------------------
  Real charlier_exact_orthog(Real p, short order)
  {
    Real value =   std::pow(p,-order)*std::exp(p)*BasisPolynomial::factorial(order);

    return value;
  }

  //------------------------------------
  // Compute numerical inner product
  //------------------------------------
  Real charlier_inner_prod(unsigned nterms, Real p, short order1, short order2, BasisPolynomial& poly)
  {
    Real sum = 0.0;
    for( short i=0; i<nterms; ++i )
      sum += std::pow(p,i)/BasisPolynomial::factorial(i)*poly.type1_value(Real(i),order1)*poly.type1_value(Real(i),order2);

    return sum;
  }
}


//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(discrete_orthog_poly, charlier1)
{
  BasisPolynomial poly_basis = BasisPolynomial(CHARLIER_DISCRETE);
  auto ptr = std::dynamic_pointer_cast<CharlierOrthogPolynomial>
    (poly_basis.polynomial_rep());
  TEST_ASSERT( ptr != NULL );

  // Test deafult settings and accessors
  Real ap;  poly_basis.pull_parameter(P_LAMBDA, ap);
  TEST_EQUALITY( ap, 0.0 );

  const Real p = 0.1;
  const Real TEST_TOL = 2.e-5; // a relative tolerance based on the exact answers
  const unsigned NUM_TERMS_TO_SUM = 100; // the number of terms needed for the orthogonality sum to converge

  poly_basis.push_parameter(P_LAMBDA, p);

  // Test orthogonality of first 7 polynomials - covers hardcoded 1st and 2nd orders and recursion-based orders
  //    NOTE:  that there appears to be some degradation in numerical precision with increasing polynomial orders
  //    such that the tolerance check has been relaxed, and possibly the number of summation terms increased
  //    We also have to stop at order 7 for now.
  for( short i=0; i<8; ++i ) {
    Real exact_orth_val = charlier_exact_orthog(p, i);
    for( short j=0; j<8; ++j ) {
      Real numerical_orth_val = charlier_inner_prod(NUM_TERMS_TO_SUM, p, i, j, *ptr);
      if( i == j ) {
        TEST_FLOATING_EQUALITY( exact_orth_val, numerical_orth_val, TEST_TOL );
      }
      else {
        Real shifted_zero = 1.0 + numerical_orth_val;
        TEST_FLOATING_EQUALITY( shifted_zero, 1.0, TEST_TOL );
      }
    }
  }
}

//----------------------------------------------------------------


namespace {

  //------------------------------------
  // Compute known exact orthogonality value
  //------------------------------------
  Real krawtchouck_exact_orthog(short N, Real p, short order)
  {
    Real value =   std::pow(-1, order)*BasisPolynomial::factorial(order)/BasisPolynomial::pochhammer(-N,order)
                 * std::pow((1.0-p)/p,order);

    return value;
  }

  //------------------------------------
  // Compute numerical inner product
  //------------------------------------
  Real krawtchouck_inner_prod(short N, Real p, short order1, short order2, BasisPolynomial& poly)
  {
    Real sum = 0.0;
    for( short i=0; i<N+1; ++i )
      sum += BasisPolynomial::n_choose_k(N,i)*std::pow(p,i)*std::pow(1.0-p,N-i)*poly.type1_value(Real(i),order1)*poly.type1_value(Real(i),order2);

    return sum;
  }
}


//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(discrete_orthog_poly, krawtchouck1)
{
  BasisPolynomial poly_basis = BasisPolynomial(KRAWTCHOUK_DISCRETE);
  auto ptr = std::dynamic_pointer_cast<KrawtchoukOrthogPolynomial>
    (poly_basis.polynomial_rep());
  TEST_ASSERT( ptr != NULL );

  // Test deafult settings and accessors
  Real ap;  poly_basis.pull_parameter(BI_P_PER_TRIAL, ap);
  unsigned int  bp;  poly_basis.pull_parameter(BI_TRIALS,      bp);
  TEST_EQUALITY( ap, 0.0 ); TEST_EQUALITY( bp, 0 );

  const Real p = 0.1;
  const unsigned int N = 15;
  const Real TEST_TOL = 1.e-9; // a relative tolerance based on the exact answers

  poly_basis.push_parameter(BI_P_PER_TRIAL, p);
  poly_basis.push_parameter(BI_TRIALS,      N);

  // Test orthogonality of first 10 polynomials - covers hardcoded 1st and 2nd orders and recursion-based orders
  for( short i=0; i<11; ++i ) {
    Real exact_orth_val = krawtchouck_exact_orthog(N, p, i);
    for( short j=0; j<11; ++j ) {
      Real numerical_orth_val = krawtchouck_inner_prod(N, p, i, j, *ptr);
      if( i == j ) {
        TEST_FLOATING_EQUALITY( exact_orth_val, numerical_orth_val, TEST_TOL );
      }
      else {
        Real shifted_zero = 1.0 + numerical_orth_val;
        TEST_FLOATING_EQUALITY( shifted_zero, 1.0, TEST_TOL );
      }
    }
  }
}

//----------------------------------------------------------------


namespace {

  //------------------------------------
  // Compute known exact orthogonality value
  //------------------------------------
  Real meixner_exact_orthog(Real c, Real B, short order)
  {
    Real value =   BasisPolynomial::factorial(order)/(BasisPolynomial::pochhammer(B,order)
                 * std::pow(c,order)*std::pow((1.0-c),B));

    return value;
  }

  //------------------------------------
  // Compute numerical inner product
  //------------------------------------
  Real meixner_inner_prod(unsigned nterms, Real c, Real B, short order1, short order2, BasisPolynomial& poly)
  {
    Real sum = 0.0;
    for( short i=0; i<nterms; ++i )
      sum += BasisPolynomial::pochhammer(B,i)*std::pow(c,i)/BasisPolynomial::factorial(i)*poly.type1_value(Real(i),order1)*poly.type1_value(Real(i),order2);

    return sum;
  }
}


//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(discrete_orthog_poly, meixner1)
{
  BasisPolynomial poly_basis = BasisPolynomial(MEIXNER_DISCRETE);
  auto ptr = std::dynamic_pointer_cast<MeixnerOrthogPolynomial>
    (poly_basis.polynomial_rep());
  TEST_ASSERT( ptr != NULL );

  // Test deafult settings and accessors
  Real ap;  poly_basis.pull_parameter(NBI_P_PER_TRIAL, ap);
  unsigned int bp;  poly_basis.pull_parameter(NBI_TRIALS,      bp);
  TEST_EQUALITY( ap, 0.0 );  TEST_EQUALITY( bp,  1 );

  const Real c = 0.1;
  const unsigned int beta = 1;//const Real beta = 1.5; // *** MSE: Real was valid for old Meixner poly but not for NBI or GE distributions
  const Real TEST_TOL = 1.e-9; // a relative tolerance based on the exact answers
  const unsigned NUM_TERMS_TO_SUM = 40; // the number of terms needed for the orthogonality sum to converge

  poly_basis.push_parameter(NBI_P_PER_TRIAL, c);
  poly_basis.push_parameter(NBI_TRIALS,      beta);

  // Test orthogonality of first 10 polynomials - covers hardcoded 1st and 2nd orders and recursion-based orders
  for( short i=0; i<7; ++i ) {
    Real exact_orth_val = meixner_exact_orthog(c, beta, i);
    for( short j=0; j<7; ++j ) {
      Real numerical_orth_val = meixner_inner_prod(NUM_TERMS_TO_SUM, c, beta, i, j, *ptr);
      if( i == j ) {
        TEST_FLOATING_EQUALITY( exact_orth_val, numerical_orth_val, TEST_TOL );
      }
      else {
        Real shifted_zero = 1.0 + numerical_orth_val;
        TEST_FLOATING_EQUALITY( shifted_zero, 1.0, TEST_TOL );
      }
    }
  }
}

//----------------------------------------------------------------


namespace {

  //------------------------------------
  // Compute known exact orthogonality value
  //------------------------------------
  Real hahn_exact_orthog(unsigned int a, unsigned int b, int N, short order)
  {
    Real value =   std::pow(-1,order)*BasisPolynomial::pochhammer((order+a+b+1.0),(N+1))*BasisPolynomial::pochhammer((b+1.0),order)
                  *BasisPolynomial::factorial(order)/((2.0*order+a+b+1.0)*BasisPolynomial::pochhammer((a+1.0),order)
                  *BasisPolynomial::pochhammer(-N,order)*BasisPolynomial::factorial(N));

    return value;
  }

  //------------------------------------
  // Compute numerical inner product
  //------------------------------------
  Real hahn_inner_prod(unsigned int N, unsigned int a, unsigned int b, short order1, short order2, BasisPolynomial& poly)
  {
    Real sum = 0.0;
    for( short i=0; i<N+1; ++i ) {
      Real x = Real(i);
      short sa = short(a);
      short sb = short(b);
      Real term = BasisPolynomial::n_choose_k(sa+i,i)*BasisPolynomial::n_choose_k(sb+N-i,N-i)*poly.type1_value(x,order1)*poly.type1_value(x,order2);
      sum += term;
    }

    return sum;
  }
}


//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(discrete_orthog_poly, hahn1)
{
  BasisPolynomial poly_basis = BasisPolynomial(HAHN_DISCRETE);
  auto ptr = std::dynamic_pointer_cast<HahnOrthogPolynomial>
    (poly_basis.polynomial_rep());
  TEST_ASSERT( ptr != NULL );

  // Test deafult settings and accessors
  unsigned int ap;  poly_basis.pull_parameter(HGE_TOT_POP, ap);
  unsigned int bp;  poly_basis.pull_parameter(HGE_SEL_POP, bp);
  unsigned int gp;  poly_basis.pull_parameter(HGE_DRAWN  , gp);
  TEST_EQUALITY( ap, 0 );
  TEST_EQUALITY( bp, 0 );
  TEST_EQUALITY( gp, 0 );

  const unsigned int totalPop   = 4;
  const unsigned int selectPop  = 6;
  const unsigned int N          = 10;
  const Real TEST_TOL  = 5.e-8; // a relative tolerance based on the exact answers

  poly_basis.push_parameter(HGE_TOT_POP, totalPop);
  poly_basis.push_parameter(HGE_SEL_POP, selectPop);
  poly_basis.push_parameter(HGE_DRAWN,   N);

  // Test orthogonality of first 10 polynomials - covers hardcoded 1st and 2nd orders and recursion-based orders
  for( short i=0; i<11; ++i ) {
    Real exact_orth_val = hahn_exact_orthog(totalPop, selectPop, N, i);
    for( short j=0; j<11; ++j ) {
      Real numerical_orth_val = hahn_inner_prod(N, totalPop, selectPop, i, j, *ptr);
      if( i == j ) {
        TEST_FLOATING_EQUALITY( exact_orth_val, numerical_orth_val, TEST_TOL );
      }
      else {
        Real shifted_zero = 1.0 + numerical_orth_val;
        TEST_FLOATING_EQUALITY( shifted_zero, 1.0, TEST_TOL );
      }
    }
  }
}

//----------------------------------------------------------------

// Test numerically-generated distributions for histogram point
// variables

namespace {

/// populate a map of int, str, real to real
template<typename ScalarType>
void array_to_map(size_t num_vals, ScalarType pt_vals[], double pt_mass[],
		  std::map<ScalarType, double>& output_map)
{
  for (size_t i=0; i<num_vals; ++i)
    output_map[pt_vals[i]] = pt_mass[i];
}


/// check pairwise orthogonality of all polynomials in the basis
template<typename ScalarType>
void histpt_check_orthog(size_t num_vals, ScalarType pt_vals[],
			 double pt_mass[], BasisPolynomial& poly_basis,
			 const unsigned short max_order, const double tol,
			 Teuchos::FancyOStream &out, bool &success)
{
  // check orthogonality to self and other orders without using the
  // inner product, since that's what we're checking...
  //
  // int_x { p_i(x) * p_j(x) * d(x) } = sum_k { p_i(x_k) * p_j(x_k) * m_k }

  //std::cerr << "i\t" << "j\t" << "integral\t" << "norm_sq" << '\n';
  for (unsigned short i = 0; i<=max_order; ++i) {
    for (unsigned short j = 0; j<=i; ++j) {
      double integral = 0.0;
      for (size_t k=0; k<num_vals; ++k) {
	double p_i_k = poly_basis.type1_value(pt_vals[k], i);
	double p_j_k = poly_basis.type1_value(pt_vals[k], j);
	integral += p_i_k * p_j_k * pt_mass[k];
      }
      if (i == j) {
	double norm_sq = poly_basis.norm_squared(i);
	TEUCHOS_TEST_FLOATING_EQUALITY( integral, norm_sq, tol, out, success );
	//std::cerr << i << "\t" << j << "\t" << integral << "\t" << norm_sq << '\n';
      }
      else {
	// shift from 0.0 to avoid numerical issues
	TEUCHOS_TEST_FLOATING_EQUALITY( 1.0 + integral, 1.0, tol, out, success );
	//std::cerr << i << "\t" << j << "\t" << integral << "\t" << 0.0 << '\n';
      }
    }
  }
}

} // namespace


TEUCHOS_UNIT_TEST(discrete_orthog_poly, hist_pt_int)
{
  BasisPolynomial poly_basis = BasisPolynomial(NUM_GEN_ORTHOG);
  auto ptr = std::dynamic_pointer_cast<NumericGenOrthogPolynomial>
    (poly_basis.polynomial_rep());
  TEST_ASSERT( ptr != NULL );

  // Test orthogonality to discrete data
  size_t num_vals = 20;
    // 20 not-so-randomly generated points (must be sorted)
  int pt_vals[] = {
    1, 2, 3, 5, 6, 8, 10, 13, 14, 18,
    21, 22, 26, 30, 34, 38, 42, 46, 50, 54
  };
  //  masses must sum to 1.0
  double pt_mass[] = {
    7.421166048564262e-02,  2.519880495959160e-02,  4.997944801838294e-02,
    6.905619480130368e-02,  8.800520252258573e-02,  9.476072279800409e-02,
    5.405504290268991e-02,  1.369359931170338e-02,  1.474756002656411e-02,
    2.543717961618541e-02,  8.304772955591133e-02,  2.511850146980441e-02,
    8.043668134029533e-02,  2.405588275825480e-02,  9.179451654564216e-02,
    3.457209536214227e-02,  1.942007146110596e-02,  2.480256493445103e-02,
    6.085412342241643e-02,  4.675241770732283e-02
  };

  IntRealMap pt_pairs;
  array_to_map(num_vals, pt_vals, pt_mass, pt_pairs);

  ptr->histogram_pt_distribution(pt_pairs);
  ptr->coefficients_norms_flag(true);

  // discrete data is more challenging
  const unsigned short max_order = 4;
  const double tol = 1.0e-6;
  histpt_check_orthog(num_vals, pt_vals, pt_mass, poly_basis, max_order, tol,
		      out, success);
}


TEUCHOS_UNIT_TEST(discrete_orthog_poly, hist_pt_str)
{
  BasisPolynomial poly_basis = BasisPolynomial(NUM_GEN_ORTHOG);
  auto ptr = std::dynamic_pointer_cast<NumericGenOrthogPolynomial>
    (poly_basis.polynomial_rep());
  TEST_ASSERT( ptr != NULL );

  // Test orthogonality to discrete data
  size_t num_vals = 20;
  // 20 not-so-randomly generated points (must be sorted)
  String pt_vals[] = {
    "aa", "bb",  "cc",  "dd",  "ee",  "ff",  "gg",  "hh",  "ii",  "jj",
    "kk",  "ll",  "mm",  "nn",  "oo",  "pp",  "qq",  "rr",  "ss",  "tt"
  };
  // string data are mapped to 0-based indices; use indices to evaluate
  int pt_inds[] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 11, 12, 13, 14, 15, 16, 17, 18, 19
  };

  //  masses must sum to 1.0
  double pt_mass[] = {
    7.421166048564262e-02,  2.519880495959160e-02,  4.997944801838294e-02,
    6.905619480130368e-02,  8.800520252258573e-02,  9.476072279800409e-02,
    5.405504290268991e-02,  1.369359931170338e-02,  1.474756002656411e-02,
    2.543717961618541e-02,  8.304772955591133e-02,  2.511850146980441e-02,
    8.043668134029533e-02,  2.405588275825480e-02,  9.179451654564216e-02,
    3.457209536214227e-02,  1.942007146110596e-02,  2.480256493445103e-02,
    6.085412342241643e-02,  4.675241770732283e-02
  };

  StringRealMap pt_pairs;
  array_to_map(num_vals, pt_vals, pt_mass, pt_pairs);

  ptr->histogram_pt_distribution(pt_pairs);
  ptr->coefficients_norms_flag(true);

  // discrete data is more challenging, since well-distributed can
  // have slightly tighter tolerance
  const unsigned short max_order = 4;
  const double tol = 1.0e-6;
  // use indices for pointwise evaluation
  histpt_check_orthog(num_vals, pt_inds, pt_mass, poly_basis, max_order, tol,
		      out, success);
}


TEUCHOS_UNIT_TEST(discrete_orthog_poly, hist_pt_real)
{
  BasisPolynomial poly_basis = BasisPolynomial(NUM_GEN_ORTHOG);
  auto ptr = std::dynamic_pointer_cast<NumericGenOrthogPolynomial>
    (poly_basis.polynomial_rep());
  TEST_ASSERT( ptr != NULL );

  // Test orthogonality to discrete data
  size_t num_vals = 20;
  // 20 randomly generated points on [0,2] (must be sorted)
  double pt_vals[] = {
    1.563510575063674e-01,  1.676427559938651e-01,  3.047560379384460e-01,
    3.243646163864855e-01,  3.312974589995619e-01,  4.579539374336377e-01,
    5.259425690802886e-01,  6.224300840896098e-01,  8.853565395508927e-01,
    9.010831970049955e-01,  1.057066271012425e+00,  1.076684870520114e+00,
    1.203963882803273e+00,  1.308158196953565e+00,  1.378429006280016e+00,
    1.496303185647419e+00,  1.588569081367814e+00,  1.651633954979095e+00,
    1.826674723003339e+00,  1.992269433253771e+00
  };
  // masses must sum to 1.0
  double pt_mass[] = {
    7.421166048564262e-02,  2.519880495959160e-02,  4.997944801838294e-02,
    6.905619480130368e-02,  8.800520252258573e-02,  9.476072279800409e-02,
    5.405504290268991e-02,  1.369359931170338e-02,  1.474756002656411e-02,
    2.543717961618541e-02,  8.304772955591133e-02,  2.511850146980441e-02,
    8.043668134029533e-02,  2.405588275825480e-02,  9.179451654564216e-02,
    3.457209536214227e-02,  1.942007146110596e-02,  2.480256493445103e-02,
    6.085412342241643e-02,  4.675241770732283e-02
  };

  RealRealMap pt_pairs;
  array_to_map(num_vals, pt_vals, pt_mass, pt_pairs);

  ptr->histogram_pt_distribution(pt_pairs);
  ptr->coefficients_norms_flag(true);

  const unsigned short max_order = 10;
  const double tol = 1.0e-12;
  histpt_check_orthog(num_vals, pt_vals, pt_mass, poly_basis, max_order, tol,
		      out, success);
}
