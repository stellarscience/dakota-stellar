/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        GenLaguerreOrthogPolynomial
//- Description:  Implementation code for GenLaguerreOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "GenLaguerreOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.hpp"
#endif
#include "pecos_stat_util.hpp"

//#define DEBUG


namespace Pecos {

Real GenLaguerreOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  // employ Horner's rule for improved efficiency and precision
  Real t1_val;
  switch (order) {
  case 0:
    t1_val = 1.;                                                          break;
  case 1:
    t1_val = -x + alphaPoly + 1.;                                         break;
  case 2: {
    Real ap2 = alphaPoly+2.;
    t1_val = (x*(x - 2.*ap2) + (alphaPoly+1)*ap2)/2.;
    break;
  }
  case 3: {
    Real ap3 = alphaPoly+3., ap2 = alphaPoly+2.;
    t1_val = (x*(x*(-x + 3.*ap3) - 3.*ap2*ap3) + (alphaPoly+1.)*ap2*ap3)/6.;
    break;
  }
  default: {
    // Support higher order polynomials using the 3 point recursion formula
    Real ap3 = alphaPoly+3., ap2 = alphaPoly+2., ap1 = alphaPoly+1.,
      La_n   = (x*(x*(-x + 3.*ap3) - 3.*ap2*ap3) + ap1*ap2*ap3)/6., // La_3
      La_nm1 = (x*(x - 2.*ap2) + ap1*ap2)/2.;                       // La_2
    for (size_t i=3; i<order; i++) {
      t1_val = ( (2.*i+1.+alphaPoly-x)*La_n - (i+alphaPoly)*La_nm1 )
	/ (i+1.); // La_np1
      if (i != order-1) {
	La_nm1 = La_n;
	La_n   = t1_val;
      }
    }
    break;
  }
  }

  return t1_val;
}


Real GenLaguerreOrthogPolynomial::type1_gradient(Real x, unsigned short order)
{
  Real t1_grad;
#ifdef DEBUG
  // See Abramowitz & Stegun, Section 22.8, p.783
  //t1_grad = (order) ? (order*type1_value(x, order)
  //                  - (order+alphaPoly)*type1_value(x, order-1))/x : 0.;
  if (order) {
    Real La_n = type1_value(x, order), La_nminus1 = type1_value(x, order-1);
    t1_grad = (order*La_n - (order+alphaPoly)*La_nminus1)/x;
  }
  else
    t1_grad = 0.;
  PCout << "Gen Laguerre gradient approach 1: " << t1_grad << '\n';
#endif // DEBUG

  // The previous approach, while very compact, produces 0/0 = NaN at x = 0.
  // To avoid NaN issue at lower bound, differentiate the 3 pt value recursion
  // to get a 3 point gradient recursion
  switch (order) {
  case 0:
    t1_grad = 0.;                                                         break;
  case 1:
    t1_grad = -1.;                                                        break;
  case 2:
    t1_grad = x - (alphaPoly + 2.);                                       break;
  case 3: {
    Real ap3 = alphaPoly+3.;
    t1_grad = (x*(-x + 2.*ap3) - (alphaPoly+2.)*ap3)/2.;
    break;
  }
  default: {
    // Support higher order polynomials using the 3 point recursion formula
    Real ap3 = alphaPoly+3., ap2 = alphaPoly+2.,
      dLadx_n   = (x*(-x + 2.*ap3) - ap2*ap3)/2., // L'a_3
      dLadx_nm1 = x - ap2;                        // L'a_2
    for (size_t i=3; i<order; i++) {
      t1_grad = ( (2.*i+1.+alphaPoly-x)*dLadx_n - type1_value(x,i) -
		  (i+alphaPoly)*dLadx_nm1 ) / (i+1.); // dLadx_np1
      if (i != order-1) {
	dLadx_nm1 = dLadx_n;
	dLadx_n   = t1_grad;
      }
    }
    break;
  }
  }
#ifdef DEBUG
  PCout << "Gen Laguerre gradient approach 2: " << t1_grad << '\n';
#endif // DEBUG

  return t1_grad;
}


Real GenLaguerreOrthogPolynomial::type1_hessian(Real x, unsigned short order)
{
  Real t1_hess;
  switch (order) {
  case 0: case 1:
    t1_hess = 0.;                                                         break;
  case 2:
    t1_hess = 1.;                                                         break;
  case 3:
    t1_hess = alphaPoly + 3. - x;                                         break;
  default: {
    // Support higher order polynomials using the 3 point recursion formula
    Real d2Ladx2_n = alphaPoly + 3. - x, d2Ladx2_nm1 = 1.; // L''a_3, L''a_2
    for (size_t i=3; i<order; i++) {
      t1_hess = ( (2.*i+1.+alphaPoly-x)*d2Ladx2_n - type1_gradient(x,i) -
		  (i+alphaPoly)*d2Ladx2_nm1 ) / (i+1.); // dLadx_np1
      if (i != order-1) {
	d2Ladx2_nm1 = d2Ladx2_n;
	d2Ladx2_n   = t1_hess;
      }
    }
    break;
  }
  }

  return t1_hess;
}


Real GenLaguerreOrthogPolynomial::norm_squared(unsigned short order)
{
  // For integer alphaPoly, Gamma(alphaPoly+n+1)/n!/Gamma(alphaPoly+1)
  // = (alphaPoly+n)!/n!/alphaPoly!
  //return n_choose_k(alphaPoly+n,n);

  // For real alphaPoly: Gamma(alphaPoly+n+1)/Gamma(alphaPoly+1)/n!
  // = (alphaPoly+1)(alphaPoly+2)...(alphaPoly+1+(n-1))/n!
  // = pochhammer(alphaPoly+1,n)/n!
  return pochhammer(alphaPoly+1,order) / factorial(order);
}


const RealArray& GenLaguerreOrthogPolynomial::
collocation_points(unsigned short order)
{
  // pull this out from default below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "GenLaguerreOrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  UShortRealArrayMap::iterator it = collocPointsMap.find(order);
  if (it != collocPointsMap.end())
    return it->second;

  RealArray& colloc_pts = collocPointsMap[order]; // create new array
  colloc_pts.resize(order);
  switch (order) {
  case 1: // zeros of L^(alphaPoly)_1(x) for one gen Gauss-Laguerre pt:
    colloc_pts[0] =  1. + alphaPoly;
    break;
  case 2: { // zeros of L^(alphaPoly)_2(x) for two gen Gauss-Laguerre pts:
    Real srap2 = sqrt(alphaPoly + 2.);
    colloc_pts[0] = alphaPoly + 2. - srap2;
    colloc_pts[1] = alphaPoly + 2. + srap2;
    break;
  }
  default:
#ifdef HAVE_SPARSE_GRID
    // sandia_rules.C calculates points/weights together
    RealArray& colloc_wts = collocWeightsMap[order];
    if (colloc_wts.size() != order)
      colloc_wts.resize(order);
    webbur::gen_laguerre_compute(order, alphaPoly, &colloc_pts[0],
				 &colloc_wts[0]);
    Real wt_factor = weight_factor();
    for (size_t i=0; i<order; i++)
      colloc_wts[i] *= wt_factor; // polynomial weight fn -> PDF
#else
    PCerr << "Error: overflow in maximum quadrature order limit (2) in "
	  << "GenLaguerreOrthogPolynomial::collocation_points().  Configure "
	  << "with VPISparseGrid to extend range." << std::endl;
    abort_handler(-1);
#endif
    break;
  }
  return colloc_pts;
}


const RealArray& GenLaguerreOrthogPolynomial::
type1_collocation_weights(unsigned short order)
{
  // Derived from -(A_{n+1} gamma_n)/(A_n Phi_n'(x_i) Phi_{n+1}(x_i)),
  // which for L^(alphaPoly)(x), is Gamma(n+alphaPoly) x_i /
  // (n! (n+alphaPoly) Gamma(1+alphaPoly) (L^(alphaPoly)_{n-1}(x_i))^2).

  // The sums of the weights = 1, which is the integral of the density function
  // x^alphaPoly exp(-x)/Gamma(alpha+1) over the support range of [0,+infinity].

  UShortRealArrayMap::iterator it = collocWeightsMap.find(order);
  if (it != collocWeightsMap.end())
    return it->second;

  RealArray& colloc_wts = collocWeightsMap[order]; // create new array
  colloc_wts.resize(order);
  switch (order) {
  case 1: // weight for one generalized Gauss-Laguerre point:
    colloc_wts[0] = 1.;
    break;
  default:
#ifdef HAVE_SPARSE_GRID
    // sandia_rules.C calculates points/weights together
    RealArray& colloc_pts = collocPointsMap[order];
    if (colloc_pts.size() != order)
      colloc_pts.resize(order);
    webbur::gen_laguerre_compute(order, alphaPoly, &colloc_pts[0],
				 &colloc_wts[0]);
    Real wt_factor = weight_factor();
    for (size_t i=0; i<order; i++)
      colloc_wts[i] *= wt_factor; // polynomial weight fn -> PDF
#else
    // define Gauss wts from Gauss pts using formula above
    const RealArray& colloc_pts = collocation_points(order);
    for (size_t i=0; i<order; i++) {
      Real x_i = colloc_pts[i];
      // For integer alphaPoly:
      //colloc_wts[i] = factorial_ratio(order+alphaPoly-1, order) * x_i /
      //  (order+alphaPoly) / factorial(alphaPoly) /
      //  std::pow(type1_value(x_i,order-1),2);
      // For real alphaPoly:
      colloc_wts[i] = pochhammer(alphaPoly+1., order)*x_i/factorial(order)
	/ std::pow((order+alphaPoly) * type1_value(x_i, order-1), 2);
    }
#endif
    break;
  }
  return colloc_wts;
}


Real GenLaguerreOrthogPolynomial::weight_factor()
{
//#ifdef HAVE_BOOST
  wtFactor = 1./bmth::tgamma(alphaPoly + 1.);
/*
#elif HAVE_GSL
  wtFactor = 1./gsl_sf_gamma(alphaPoly + 1.);
#else
  PCerr << "Error: BOOST or GSL required in GenLaguerreOrthogPolynomial::"
        << "weight_factor()." << std::endl;
  abort_handler(-1);
#endif
*/
  return wtFactor;
}

} // namespace Pecos
