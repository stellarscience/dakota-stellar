/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LegendreOrthogPolynomial
//- Description:  Implementation code for LegendreOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "LegendreOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.hpp"
#endif

//#define DEBUG


namespace Pecos {

Real LegendreOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  // employ Horner's rule for improved efficiency and precision
  Real t1_val;
  switch (order) {
  case 0:
    t1_val = 1.;                                                          break;
  case 1:
    t1_val = x;                                                           break;
  case 2:
    t1_val = (3.*x*x - 1.)/2.;                                            break;
  case 3:
    t1_val = x*(5.*x*x - 3.)/2.;                                          break;
  case 4: {
    Real x2 = x*x;
    t1_val = (x2*(35.*x2 - 30.) + 3.)/8.;                                 break;
  }
  case 5: {
    Real x2 = x*x;
    t1_val = x*(x2*(63.*x2 - 70.) + 15.)/8.;                              break;
  }
  case 6: {
    Real x2 = x*x;
    t1_val = (x2*(x2*(231.*x2 - 315.) + 105.) - 5.)/16.;                  break;
  }
  case 7: {
    Real x2 = x*x;
    t1_val = x*(x2*(x2*(429.*x2 - 693.) + 315.) - 35.)/16.;               break;
  }
  case 8: {
    Real x2 = x*x;
    t1_val = (x2*(x2*(x2*(6435.*x2 - 12012.) + 6930.) - 1260.) + 35.)/128.;
    break;
  }
  case 9: {
    Real x2 = x*x;
    t1_val = x*(x2*(x2*(x2*(12155.*x2 - 25740.) + 18018.) - 4620.) + 315.)/128.;
    break;
  }
  case 10: {
    Real x2 = x*x;
    t1_val  = (x2*(x2*(x2*(x2*(46189.*x2 - 109395.) + 90090.) - 30030.)
	    + 3465.) - 63.)/256.;
    break;
  }
  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real x2 = x*x, P_n = (x2*(x2*(x2*(x2*(46189.*x2 - 109395.) + 90090.)
		       - 30030.) + 3465.) - 63.)/256., // P_10
      P_nminus1 = x*(x2*(x2*(x2*(12155.*x2 - 25740.) + 18018.) - 4620.) + 315.)
		/ 128.;                                // P_9
    for (size_t i=10; i<order; i++) {
      t1_val = ( (2.*i+1.)*x*P_n - i*P_nminus1 ) / (i+1.); // P_nplus1
      if (i != order-1) {
	P_nminus1 = P_n;
	P_n       = t1_val;
      }
    }
    break;
  }
  }

  return t1_val;
}


Real LegendreOrthogPolynomial::type1_gradient(Real x, unsigned short order)
{
  Real t1_grad;
#ifdef DEBUG
  // See Abramowitz & Stegun, Section 22.8, p.783
  if (order) {
    Real P_n = type1_value(x, order), P_nminus1 = type1_value(x, order-1);
    t1_grad = order*(x*P_n - P_nminus1)/(x*x - 1.);
  }
  else
    t1_grad = 0.;
  PCout << "Legendre gradient approach 1: " << t1_grad << '\n';
#endif // DEBUG

  // The previous approach, while very compact, produces 0/0 = NaN at x = +/-1.
  // It is therefore only used for cross-validation purposes.  To avoid NaN
  // issues at bounds, differentiate the 3 pt value recursion to get a 3 point
  // gradient recursion
  switch (order) {
  case 0:
    t1_grad = 0.;                                                         break;
  case 1:
    t1_grad = 1;                                                          break;
  case 2:
    t1_grad = 3.*x;                                                       break;
  case 3:
    t1_grad = (15.*x*x - 3.)/2.;                                          break;
  case 4:
    t1_grad = x*(35.*x*x - 15.)/2.;                                       break;
  case 5: {
    Real x2 = x*x;
    t1_grad = (x2*(315.*x2 - 210.) + 15.)/8.;                             break;
  }
  case 6: {
    Real x2 = x*x;
    t1_grad = x*(x2*(693.*x2 - 630.) + 105.)/8.;                          break;
  }
  default:
    // Support higher order polynomials using a 3 point recursion formula:
    Real x2 = x*x, dPdx_n = x*(x2*(693.*x2 - 630.) + 105.)/8., // P'_6
      dPdx_nminus1 = (x2*(315.*x2 - 210.) + 15.)/8.;           // P'_5
    for (size_t i=6; i<order; i++) {
      t1_grad // dPdx_nplus1
	= ( (2.*i+1.)*(x*dPdx_n + type1_value(x,i)) - i*dPdx_nminus1 ) / (i+1.);
      if (i != order-1) {
	dPdx_nminus1 = dPdx_n;
	dPdx_n       = t1_grad;
      }
    }
    break;
  }
#ifdef DEBUG
  PCout << "Legendre gradient approach 2: " << t1_grad << '\n';
#endif // DEBUG

  return t1_grad;
}


Real LegendreOrthogPolynomial::type1_hessian(Real x, unsigned short order)
{
  Real t1_hess;
  switch (order) {
  case 0: case 1:
    t1_hess = 0.;                                                         break;
  case 2:
    t1_hess = 3.;                                                         break;
  case 3:
    t1_hess = 15.*x;                                                      break;
  case 4:
    t1_hess = (105.*x*x - 15.)/2.;                                        break;
  case 5:
    t1_hess = x*(315.*x*x - 105.)/2.;                                     break;
  case 6: {
    Real x2 = x*x;
    t1_hess = (x2*(3465.*x2 - 1890.) + 105.)/8.;                          break;
  }
  default:
    // Support higher order polynomials using a 3 point recursion formula:
    Real x2 = x*x, d2Pdx2_n = (x2*(3465.*x2 - 1890.) + 105.)/8., // P''_6
      d2Pdx2_nminus1 = x*(315.*x*x - 105.)/2.;                   // P''_5
    for (size_t i=6; i<order; i++) {
      t1_hess = ((2.*i+1.)*(x*d2Pdx2_n + 2.*type1_gradient(x,i)) -
		 i*d2Pdx2_nminus1) / (i+1.); // d2Pdx2_nplus1
      if (i != order-1) {
	d2Pdx2_nminus1 = d2Pdx2_n;
	d2Pdx2_n       = t1_hess;
      }
    }
    break;
  }

  return t1_hess;
}


Real LegendreOrthogPolynomial::norm_squared(unsigned short order)
{
  // Abramowitz & Stegun: w(x) = 1
  //return 2./(2.*order + 1.);

  // sampling density f(x) = 1/(U-L) = 1/2 for [L,U] = [-1,1]
  return 1./(2.*order + 1.);
}


const RealArray& LegendreOrthogPolynomial::
collocation_points(unsigned short order)
{
  // pull this outside block below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "LegendreOrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  bool rule_err = false;
  if (collocPoints.size() != order) { // if not already computed
    collocPoints.resize(order);
    switch (collocRule) {
    case GAUSS_PATTERSON:
#ifdef HAVE_SPARSE_GRID
      webbur::patterson_lookup_points(order, &collocPoints[0]);
#else
      rule_err = true;
#endif
      break;
    case CLENSHAW_CURTIS: 
#ifdef HAVE_SPARSE_GRID
      webbur::clenshaw_curtis_compute_points(order, &collocPoints[0]);
#else
      rule_err = true;
#endif
      break;
    case FEJER2: 
#ifdef HAVE_SPARSE_GRID
      webbur::fejer2_compute_points(order, &collocPoints[0]);
#else
      rule_err = true;
#endif
      break;
    case GAUSS_LEGENDRE:
#ifdef HAVE_SPARSE_GRID
      if (order <= 33) // retrieve full precision tabulated values
	webbur::legendre_lookup_points(order, &collocPoints[0]);
      else { // sandia_rules calculates points/weights together
	if (collocWeights.size() != order)
	  collocWeights.resize(order);
	webbur::legendre_compute(order, &collocPoints[0], &collocWeights[0]);
	for (size_t i=0; i<order; i++)
	  collocWeights[i] *= wtFactor; // polynomial weight fn -> PDF
      }
#else
      switch (order) {
      case 1: // zeros of P_1(x) for one Gauss-Legendre point:
	collocPoints[0] = 0.0;
	break;
      case 2: { // zeros of P_2(x) for two Gauss-Legendre points:
	Real z1osr3 = 1./std::sqrt(3.);
	collocPoints[0] = -z1osr3;
	collocPoints[1] =  z1osr3;
	break;
      }
      case 3: { // zeros of P_3(x) for three Gauss-Legendre points:
	Real sr3o5 = std::sqrt(3./5.);
	collocPoints[0] = -sr3o5;
	collocPoints[1] =  0.0;
	collocPoints[2] =  sr3o5;
	break;
      }
      case 4: { // zeros of P_4(x) for four Gauss-Legendre points:
	Real sr30 = std::sqrt(30.), sr525p70sr30 = std::sqrt(525.+70.*sr30)/35.,
	  sr525m70sr30 = std::sqrt(525.-70.*sr30)/35.;
	collocPoints[0] = -sr525p70sr30;
	collocPoints[1] = -sr525m70sr30;
	collocPoints[2] =  sr525m70sr30;
	collocPoints[3] =  sr525p70sr30;
	break;
      }
      case 5: { // zeros of P_5(x) for five Gauss-Legendre points:
	Real sr70 = std::sqrt(70.), sr245p14sr70 = std::sqrt(245.+14.*sr70)/21.,
	  sr245m14sr70 = std::sqrt(245.-14.*sr70)/21.;
	collocPoints[0] = -sr245p14sr70;
	collocPoints[1] = -sr245m14sr70;
	collocPoints[2] =  0.0;
	collocPoints[3] =  sr245m14sr70;
	collocPoints[4] =  sr245p14sr70;
	break;
      }
      // tabulated values from Abramowitz & Stegun have limited precision
      case 6:
	collocPoints[0] = -0.932469514203152;
	collocPoints[1] = -0.661209386466265;
	collocPoints[2] = -0.238619186083197;
	collocPoints[3] = -collocPoints[2];
	collocPoints[4] = -collocPoints[1];
	collocPoints[5] = -collocPoints[0]; break;
      case 7:
	collocPoints[0] = -0.949107912342759;
	collocPoints[1] = -0.741531185599394;
	collocPoints[2] = -0.405845151377397;
	collocPoints[3] =  0.0;
	collocPoints[4] = -collocPoints[2];
	collocPoints[5] = -collocPoints[1];
	collocPoints[6] = -collocPoints[0]; break;
      case 8:
	collocPoints[0] = -0.960289856497536;
	collocPoints[1] = -0.796666477413627;
	collocPoints[2] = -0.525532409916329;
	collocPoints[3] = -0.183434642495650;
	collocPoints[4] = -collocPoints[3];
	collocPoints[5] = -collocPoints[2];
	collocPoints[6] = -collocPoints[1];
	collocPoints[7] = -collocPoints[0]; break;
      case 9:
	collocPoints[0] = -0.968160239507626;
	collocPoints[1] = -0.836031107326636;
	collocPoints[2] = -0.613371432700590;
	collocPoints[3] = -0.324253423403809;
	collocPoints[4] =  0.0;
	collocPoints[5] = -collocPoints[3];
	collocPoints[6] = -collocPoints[2];
	collocPoints[7] = -collocPoints[1];
	collocPoints[8] = -collocPoints[0]; break;
      case 10:
	collocPoints[0] = -0.973906528517172;
	collocPoints[1] = -0.865063366688985;
	collocPoints[2] = -0.679409568299024;
	collocPoints[3] = -0.433395394129247;
	collocPoints[4] = -0.148874338981631;
	collocPoints[5] = -collocPoints[4];
	collocPoints[6] = -collocPoints[3];
	collocPoints[7] = -collocPoints[2];
	collocPoints[8] = -collocPoints[1];
	collocPoints[9] = -collocPoints[0]; break;
      default:
	PCerr << "Error: overflow in maximum quadrature order limit (10) in "
	      << "LegendreOrthogPolynomial::collocation_points().  Configure "
	      << "with VPISparseGrid to extend range." << std::endl;
	abort_handler(-1); break;
      }
#endif
      break;
    default:
      rule_err = true; break;
    }
  }

  if (rule_err) {
    PCerr << "Error: unsupported collocation rule in LegendreOrthogPolynomial"
	  << "::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  return collocPoints;
}


const RealArray& LegendreOrthogPolynomial::
type1_collocation_weights(unsigned short order)
{
  // pull this outside block below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in Legendre"
	  << "OrthogPolynomial::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  // The sums of the weights = 1, which is the integral of the density
  // function 1/2 over the support range of [-1,+1].  These differ from
  // Abramowitz & Stegun by a constant factor of 1/2.

  bool rule_err = false;
  if (collocWeights.size() != order) { // if not already computed
    collocWeights.resize(order);
    switch (collocRule) {
    case GAUSS_PATTERSON:
#ifdef HAVE_SPARSE_GRID
      webbur::patterson_lookup_weights(order, &collocWeights[0]);
#else
      rule_err = true;
#endif
      break;
    case CLENSHAW_CURTIS: 
#ifdef HAVE_SPARSE_GRID
      webbur::clenshaw_curtis_compute_weights(order, &collocWeights[0]);
#else
      rule_err = true;
#endif
      break;
    case FEJER2: 
#ifdef HAVE_SPARSE_GRID
      webbur::fejer2_compute_weights(order, &collocWeights[0]);
#else
      rule_err = true;
#endif
      break;
    case GAUSS_LEGENDRE:
#ifdef HAVE_SPARSE_GRID
      if (order <= 33) // retrieve full precision tabulated values
	webbur::legendre_lookup_weights(order, &collocWeights[0]);
      else { // sandia_rules calculates points/weights together
	if (collocPoints.size() != order)
	  collocPoints.resize(order);
	webbur::legendre_compute(order, &collocPoints[0], &collocWeights[0]);
      }
#else
      switch (order) {
      case 1: // weights for one Gauss-Legendre point:
	collocWeights[0] = 1.0; break;
      case 2: // weights for two Gauss-Legendre points:
	collocWeights[0] = collocWeights[1] = 0.5; break;
      case 3: // weights for three Gauss-Legendre points:
	collocWeights[0] = collocWeights[2] = 5./18.;
	collocWeights[1] = 4./9.; break;
      case 4: { // weights for four Gauss-Legendre points:
	Real sr30 = std::sqrt(30.);
	collocWeights[0] = collocWeights[3] = (18.-sr30)/72.;
	collocWeights[1] = collocWeights[2] = (18.+sr30)/72.; break;
      }
      case 5: { // weights for five Gauss-Legendre points:
	Real sr70 = std::sqrt(70.);
	collocWeights[0] = collocWeights[4] = (322.-13.*sr70)/1800.;
	collocWeights[1] = collocWeights[3] = (322.+13.*sr70)/1800.;
	collocWeights[2] = 64./225.; break;
      }
      // tabulated values from Abramowitz & Stegun have limited precision
      case 6:
	collocWeights[0] = collocWeights[5] = 0.171324492379170 * wtFactor;
	collocWeights[1] = collocWeights[4] = 0.360761573048139 * wtFactor;
	collocWeights[2] = collocWeights[3] = 0.467913934572691 * wtFactor;
	break;
      case 7:
	collocWeights[0] = collocWeights[6] = 0.129484966168870 * wtFactor;
	collocWeights[1] = collocWeights[5] = 0.279705391489277 * wtFactor;
	collocWeights[2] = collocWeights[4] = 0.381830050505119 * wtFactor;
	collocWeights[3] = 0.417959183673469 * wtFactor;	break;
      case 8:
	collocWeights[0] = collocWeights[7] = 0.101228536290376 * wtFactor;
	collocWeights[1] = collocWeights[6] = 0.222381034453374 * wtFactor;
	collocWeights[2] = collocWeights[5] = 0.313706645877887 * wtFactor;
	collocWeights[3] = collocWeights[4] = 0.362683783378362 * wtFactor;
	break;
      case 9:
	collocWeights[0] = collocWeights[8] = 0.081274388361574 * wtFactor;
	collocWeights[1] = collocWeights[7] = 0.180648160694857 * wtFactor;
	collocWeights[2] = collocWeights[6] = 0.260610696402935 * wtFactor;
	collocWeights[3] = collocWeights[5] = 0.312347077040003 * wtFactor;
	collocWeights[4] = 0.330239355001260 * wtFactor;	break;
      case 10:
	collocWeights[0] = collocWeights[9] = 0.066671344308688 * wtFactor;
	collocWeights[1] = collocWeights[8] = 0.149451349150581 * wtFactor;
	collocWeights[2] = collocWeights[7] = 0.219086362515982 * wtFactor;
	collocWeights[3] = collocWeights[6] = 0.269266719309996 * wtFactor;
	collocWeights[4] = collocWeights[5] = 0.295524224714753 * wtFactor;
	break;
      default:
	// define Gauss wts from Gauss pts using
	// -(A_{n+1} gamma_n)/(A_n Phi_n'(x_i) Phi_{n+1}(x_i)),
	// which for P(x) with w(x) = 1/2 is (1-x_i^2)/(n P_{n-1}(x_i))^2.
	const RealArray& colloc_pts = collocation_points(order); Real x_i;
	for (size_t i=0; i<order; i++) {
	  x_i = colloc_pts[i];
	  collocWeights[i] = (1.-x_i*x_i)
	                  / std::pow(order*type1_value(x_i,order-1),2);
	}
	break;
      }
#endif
      break;
    default:
      rule_err = true; break;
    }

#ifdef HAVE_SPARSE_GRID
    for (size_t i=0; i<order; i++)
      collocWeights[i] *= wtFactor;
#endif
  }

  if (rule_err) {
    PCerr << "Error: unsupported collocation rule in LegendreOrthogPolynomial"
	  << "::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  return collocWeights;
}

} // namespace Pecos
