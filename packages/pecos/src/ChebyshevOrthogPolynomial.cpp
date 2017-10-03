/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        ChebyshevOrthogPolynomial
//- Description:  Implementation code for ChebyshevOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "ChebyshevOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.hpp"
#endif

//#define DEBUG


namespace Pecos {

Real ChebyshevOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  Real t1_val;
  switch (order) {
  case 0:
    t1_val = 1.;                                                          break;
  case 1:
    t1_val = x;                                                           break;
  case 2:
    t1_val = 2.*x*x - 1.;                                                 break;
  case 3:
    t1_val = x*(4.*x*x - 3.);                                             break;
  case 4: {
    Real x2 = x*x;
    t1_val = 8.*x2*(x2 - 1.) + 1.;                                        break;
  }
  case 5: {
    Real x2 = x*x;
    t1_val = x*(x2*(16.*x2 - 20.) + 5.);                                  break;
  }
  case 6: {
    Real x2 = x*x;
    t1_val = x2*(x2*(32.*x2 - 48.) + 18.) - 1.;                           break;
  }
  case 7: {
    Real x2 = x*x;
    t1_val = x*(x2*(x2*(64.*x2 - 112.) + 56.) - 7.);                      break;
  }
  case 8: {
    Real x2 = x*x;
    t1_val = x2*(x2*(x2*(128.*x2 - 256.) + 160.) - 32.) + 1.;             break;
  }
  case 9: {
    Real x2 = x*x;
    t1_val = x*(x2*(x2*(x2*(256.*x2 - 576.) + 432.) - 120.) + 9.);        break;
  }
  default:
    // Support higher order polynomials using the 3 point recursion formula:
    Real x2 = x*x,
      T_n       = x*(x2*(x2*(x2*(256.*x2 - 576.) + 432.) - 120.) + 9.), // T_9
      T_nminus1 = x2*(x2*(x2*(128.*x2 - 256.) + 160.) - 32.) + 1.;      // T_8
    for (size_t i=9; i<order; i++) {
      t1_val = 2.*x*T_n - T_nminus1; // T_nplus1
      if (i != order-1) {
	T_nminus1 = T_n;
	T_n       = t1_val;
      }
    }
    break;
  }

  return t1_val;
}


Real ChebyshevOrthogPolynomial::type1_gradient(Real x, unsigned short order)
{
  Real t1_grad;
  switch (order) {
  case 0:
    t1_grad = 0.;                                                         break;
  case 1:
    t1_grad = 1;                                                          break;
  case 2:
    t1_grad = 4.*x;                                                       break;
  case 3:
    t1_grad = 12.*x*x - 3.;                                               break;
  case 4:
    t1_grad = x*(32.*x*x - 16.);                                          break;
  case 5: {
    Real x2 = x*x;
    t1_grad = x2*(80.*x2 - 60.) + 5.;                                     break;
  }
  case 6: {
    Real x2 = x*x;
    t1_grad = x*(x2*(192.*x2 - 192.) + 36.);                              break;
  }
  case 7: {
    Real x2 = x*x;
    t1_grad = x2*(x2*(448.*x2 - 560.) + 168.) - 7.;                       break;
  }
  case 8: {
    Real x2 = x*x;
    t1_grad = x*(x2*(x2*(1024.*x2 - 1536.) + 640.) - 64.);                break;
  }
  case 9: {
    Real x2 = x*x;
    t1_grad = x2*(x2*(x2*(2304.*x2 - 4032.) + 2160.) - 360.) + 9.;        break;
  }
  default:
    // Support higher order polynomials using a 3 point recursion formula:
    Real x2 = x*x,
      dTdx_n = x2*(x2*(x2*(2304.*x2 - 4032.) + 2160.) - 360.) + 9., // P'_9
      dTdx_nminus1 = x*(x2*(x2*(1024.*x2 - 1536.) + 640.) - 64.);   // P'_8
    for (size_t i=9; i<order; i++) {
      // dTdx_nplus1:
      t1_grad = 2.*x*dTdx_n + 2.*type1_value(x,i) - dTdx_nminus1;
      if (i != order-1) {
	dTdx_nminus1 = dTdx_n;
	dTdx_n       = t1_grad;
      }
    }
    break;
  }

  return t1_grad;
}


Real ChebyshevOrthogPolynomial::type1_hessian(Real x, unsigned short order)
{
  Real t1_hess;
  switch (order) {
  case 0: case 1:
    t1_hess = 0.;                                                         break;
  case 2:
    t1_hess = 4.;                                                         break;
  case 3:
    t1_hess = 24.*x;                                                      break;
  case 4:
    t1_hess = 96.*x*x - 16.;                                              break;
  case 5:
    t1_hess = x*(320.*x*x - 120.);                                        break;
  case 6: {
    Real x2 = x*x;
    t1_hess = x2*(960.*x2 - 576.) + 36.;                                  break;
  }
  case 7: {
    Real x2 = x*x;
    t1_hess = x*(x2*(2688.*x2 - 2240.) + 336.);                           break;
  }
  case 8: {
    Real x2 = x*x;
    t1_hess = x2*(x2*(7168.*x2 - 7680.) + 1920.) - 64.;                   break;
  }
  case 9: {
    Real x2 = x*x;
    t1_hess = x*(x2*(x2*(18432.*x2 - 24192.) + 8640.) - 720.);        break;
  }
  default:
    // Support higher order polynomials using a 3 point recursion formula:
    Real x2 = x*x,
      d2Tdx2_n = x*(x2*(x2*(18432.*x2 - 24192.) + 8640.) - 720.), // P'_9
      d2Tdx2_nminus1 = x2*(x2*(7168.*x2 - 7680.) + 1920.) - 64.;   // P'_8
    for (size_t i=9; i<order; i++) {
      // d2Tdx2_nplus1:
      t1_hess = 2.*x*d2Tdx2_n + 4.*type1_gradient(x,i) - d2Tdx2_nminus1;
      if (i != order-1) {
	d2Tdx2_nminus1 = d2Tdx2_n;
	d2Tdx2_n       = t1_hess;
      }
    }
    break;
  }

  return t1_hess;
}


Real ChebyshevOrthogPolynomial::norm_squared(unsigned short order)
{ return (order) ? PI/2. : PI; }


const RealArray& ChebyshevOrthogPolynomial::
collocation_points(unsigned short order)
{
  // pull this out from default below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in Chebyshev"
	  << "OrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  if (collocPoints.size() != order) { // if not already computed
    collocPoints.resize(order);

#ifdef HAVE_SPARSE_GRID
    // separable calculation of points/weights in sandia_rules.C
    if (collocRule == CLENSHAW_CURTIS)
      webbur::clenshaw_curtis_compute_points(order, &collocPoints[0]);
    else if (collocRule == FEJER2)
      webbur::fejer2_compute_points(order, &collocPoints[0]);
    else {
      PCerr << "Error: unsupported collocation point type in ChebyshevOrthog"
	    << "Polynomial::collocation_points()." << std::endl;
      abort_handler(-1);
    }
#else
    PCerr << "Error: configuration with VPISparseGrid package required in "
	  << "ChebyshevOrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
#endif
  }

  return collocPoints;
}


const RealArray& ChebyshevOrthogPolynomial::
type1_collocation_weights(unsigned short order)
{
  // The sums of the weights = 1, which is the integral of the density
  // function 1/2 over the support range of [-1,+1].  These differ from
  // VPISparseGrid by a constant factor of 1/2.

  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in Chebyshev"
	  << "OrthogPolynomial::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  if (collocWeights.size() != order) { // if not already computed
    collocWeights.resize(order);

#ifdef HAVE_SPARSE_GRID
    // separable calculation of points/weights in sandia_rules.C
    if (collocRule == CLENSHAW_CURTIS)
      webbur::clenshaw_curtis_compute_weights(order, &collocWeights[0]);
    else if (collocRule == FEJER2)
      webbur::fejer2_compute_weights(order, &collocWeights[0]);
    else {
      PCerr << "Error: unsupported collocation weight type in ChebyshevOrthog"
	    << "Polynomial::type1_collocation_weights()." << std::endl;
      abort_handler(-1);
    }
    for (size_t i=0; i<order; i++)
      collocWeights[i] *= wtFactor;
#else
    PCerr << "Error: configuration with VPISparseGrid package required in "
	  << "ChebyshevOrthogPolynomial::type1_collocation_weights()."
	  << std::endl;
    abort_handler(-1);
#endif
  }

  return collocWeights;
}

} // namespace Pecos
