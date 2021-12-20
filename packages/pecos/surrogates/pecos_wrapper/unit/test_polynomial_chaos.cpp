/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "Teuchos_UnitTestHarness.hpp"
//#include "PolynomialChaosExpansionWrap.hpp"
#include "PolynomialChaosExpansion.hpp"
#include "RegressionBuilder.hpp"

using namespace Pecos;
using namespace Pecos::util;
using namespace Pecos::surrogates;

namespace {

bool test_polynomial_approximation(short basis_type, Real tol){
//PolynomialChaosExpansionWrap poly;
  PolynomialChaosExpansion poly;
  RealMatrix samples, basis_matrix;
  size_t nvars=2, nsamples=5, seed=1;
  get_canonical_uniform_samples(nvars, nsamples, seed, samples);
  samples.print(std::cout);
  Pecos::ShortArray basis_types(nvars,Pecos::LEGENDRE_ORTHOG);
  poly.initialize_polynomial_basis_from_basis_types(basis_types);
  poly.generate_basis_matrix(samples, basis_matrix);
  basis_matrix.print(std::cout);
  return true;
}

TEUCHOS_UNIT_TEST(pecos_polynomial_approximation, pce)
{
  TEST_ASSERT(test_polynomial_approximation(Pecos::LEGENDRE_ORTHOG, 1e-8));
}

} // namespace
