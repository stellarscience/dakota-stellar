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

  UShortRealArrayMap::iterator it = collocPointsMap.find(order);
  if (it != collocPointsMap.end())
    return it->second;

  RealArray& colloc_pts = collocPointsMap[order]; // create new array
  colloc_pts.resize(order);
  bool rule_err = false;
  switch (collocRule) {
  case GAUSS_PATTERSON:
#ifdef HAVE_SPARSE_GRID
    webbur::patterson_lookup_points(order, &colloc_pts[0]);
#else
    rule_err = true;
#endif
    break;
  case CLENSHAW_CURTIS: 
#ifdef HAVE_SPARSE_GRID
    webbur::clenshaw_curtis_compute_points(order, &colloc_pts[0]);
#else
    rule_err = true;
#endif
    break;
  case FEJER2: 
#ifdef HAVE_SPARSE_GRID
    webbur::fejer2_compute_points(order, &colloc_pts[0]);
#else
    rule_err = true;
#endif
    break;
  case GAUSS_LEGENDRE:
#ifdef HAVE_SPARSE_GRID
    if (order <= 33) // retrieve full precision tabulated values
      webbur::legendre_lookup_points(order, &colloc_pts[0]);
    else { // sandia_rules calculates points/weights together
      RealArray& colloc_wts = collocWeightsMap[order];
      if (colloc_wts.size() != order)
	colloc_wts.resize(order);
      webbur::legendre_compute(order, &colloc_pts[0], &colloc_wts[0]);
      for (size_t i=0; i<order; i++)
	colloc_wts[i] *= wtFactor; // polynomial weight fn -> PDF
    }
#else
    switch (order) {
    case 1: // zeros of P_1(x) for one Gauss-Legendre point:
      colloc_pts[0] = 0.0;
      break;
    case 2: { // zeros of P_2(x) for two Gauss-Legendre points:
      Real z1osr3 = 1./std::sqrt(3.);
      colloc_pts[0] = -z1osr3;
      colloc_pts[1] =  z1osr3;
      break;
    }
    case 3: { // zeros of P_3(x) for three Gauss-Legendre points:
      Real sr3o5 = std::sqrt(3./5.);
      colloc_pts[0] = -sr3o5;
      colloc_pts[1] =  0.0;
      colloc_pts[2] =  sr3o5;
      break;
    }
    case 4: { // zeros of P_4(x) for four Gauss-Legendre points:
      Real sr30 = std::sqrt(30.), sr525p70sr30 = std::sqrt(525.+70.*sr30)/35.,
	sr525m70sr30 = std::sqrt(525.-70.*sr30)/35.;
      colloc_pts[0] = -sr525p70sr30;
      colloc_pts[1] = -sr525m70sr30;
      colloc_pts[2] =  sr525m70sr30;
      colloc_pts[3] =  sr525p70sr30;
      break;
    }
    case 5: { // zeros of P_5(x) for five Gauss-Legendre points:
      Real sr70 = std::sqrt(70.), sr245p14sr70 = std::sqrt(245.+14.*sr70)/21.,
	sr245m14sr70 = std::sqrt(245.-14.*sr70)/21.;
      colloc_pts[0] = -sr245p14sr70;
      colloc_pts[1] = -sr245m14sr70;
      colloc_pts[2] =  0.0;
      colloc_pts[3] =  sr245m14sr70;
      colloc_pts[4] =  sr245p14sr70;
      break;
    }
    // tabulated values from Abramowitz & Stegun have limited precision
    case 6:
      colloc_pts[0] = -0.932469514203152;
      colloc_pts[1] = -0.661209386466265;
      colloc_pts[2] = -0.238619186083197;
      colloc_pts[3] = -colloc_pts[2];
      colloc_pts[4] = -colloc_pts[1];
      colloc_pts[5] = -colloc_pts[0]; break;
    case 7:
      colloc_pts[0] = -0.949107912342759;
      colloc_pts[1] = -0.741531185599394;
      colloc_pts[2] = -0.405845151377397;
      colloc_pts[3] =  0.0;
      colloc_pts[4] = -colloc_pts[2];
      colloc_pts[5] = -colloc_pts[1];
      colloc_pts[6] = -colloc_pts[0]; break;
    case 8:
      colloc_pts[0] = -0.960289856497536;
      colloc_pts[1] = -0.796666477413627;
      colloc_pts[2] = -0.525532409916329;
      colloc_pts[3] = -0.183434642495650;
      colloc_pts[4] = -colloc_pts[3];
      colloc_pts[5] = -colloc_pts[2];
      colloc_pts[6] = -colloc_pts[1];
      colloc_pts[7] = -colloc_pts[0]; break;
    case 9:
      colloc_pts[0] = -0.968160239507626;
      colloc_pts[1] = -0.836031107326636;
      colloc_pts[2] = -0.613371432700590;
      colloc_pts[3] = -0.324253423403809;
      colloc_pts[4] =  0.0;
      colloc_pts[5] = -colloc_pts[3];
      colloc_pts[6] = -colloc_pts[2];
      colloc_pts[7] = -colloc_pts[1];
      colloc_pts[8] = -colloc_pts[0]; break;
    case 10:
      colloc_pts[0] = -0.973906528517172;
      colloc_pts[1] = -0.865063366688985;
      colloc_pts[2] = -0.679409568299024;
      colloc_pts[3] = -0.433395394129247;
      colloc_pts[4] = -0.148874338981631;
      colloc_pts[5] = -colloc_pts[4];
      colloc_pts[6] = -colloc_pts[3];
      colloc_pts[7] = -colloc_pts[2];
      colloc_pts[8] = -colloc_pts[1];
      colloc_pts[9] = -colloc_pts[0]; break;
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

  if (rule_err) {
    PCerr << "Error: unsupported collocation rule in LegendreOrthogPolynomial"
	  << "::collocation_points()." << std::endl;
    abort_handler(-1);
  }
  return colloc_pts;
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

  UShortRealArrayMap::iterator it = collocWeightsMap.find(order);
  if (it != collocWeightsMap.end())
    return it->second;
  
  RealArray& colloc_wts = collocWeightsMap[order]; // create new array
  colloc_wts.resize(order);
  bool rule_err = false;
  switch (collocRule) {
  case GAUSS_PATTERSON:
#ifdef HAVE_SPARSE_GRID
    webbur::patterson_lookup_weights(order, &colloc_wts[0]);
#else
    rule_err = true;
#endif
    break;
  case CLENSHAW_CURTIS: 
#ifdef HAVE_SPARSE_GRID
    webbur::clenshaw_curtis_compute_weights(order, &colloc_wts[0]);
#else
    rule_err = true;
#endif
    break;
  case FEJER2: 
#ifdef HAVE_SPARSE_GRID
    webbur::fejer2_compute_weights(order, &colloc_wts[0]);
#else
    rule_err = true;
#endif
    break;
  case GAUSS_LEGENDRE:
#ifdef HAVE_SPARSE_GRID
    if (order <= 33) // retrieve full precision tabulated values
      webbur::legendre_lookup_weights(order, &colloc_wts[0]);
    else { // sandia_rules calculates points/weights together
      RealArray& colloc_pts = collocPointsMap[order];
      if (colloc_pts.size() != order)
	colloc_pts.resize(order);
      webbur::legendre_compute(order, &colloc_pts[0], &colloc_wts[0]);
    }
#else
    switch (order) {
    case 1: // weights for one Gauss-Legendre point:
      colloc_wts[0] = 1.0; break;
    case 2: // weights for two Gauss-Legendre points:
      colloc_wts[0] = colloc_wts[1] = 0.5; break;
    case 3: // weights for three Gauss-Legendre points:
      colloc_wts[0] = colloc_wts[2] = 5./18.;
      colloc_wts[1] = 4./9.; break;
    case 4: { // weights for four Gauss-Legendre points:
      Real sr30 = std::sqrt(30.);
      colloc_wts[0] = colloc_wts[3] = (18.-sr30)/72.;
      colloc_wts[1] = colloc_wts[2] = (18.+sr30)/72.; break;
    }
    case 5: { // weights for five Gauss-Legendre points:
      Real sr70 = std::sqrt(70.);
      colloc_wts[0] = colloc_wts[4] = (322.-13.*sr70)/1800.;
      colloc_wts[1] = colloc_wts[3] = (322.+13.*sr70)/1800.;
      colloc_wts[2] = 64./225.; break;
    }
    // tabulated values from Abramowitz & Stegun have limited precision
    case 6:
      colloc_wts[0] = colloc_wts[5] = 0.171324492379170 * wtFactor;
      colloc_wts[1] = colloc_wts[4] = 0.360761573048139 * wtFactor;
      colloc_wts[2] = colloc_wts[3] = 0.467913934572691 * wtFactor;
      break;
    case 7:
      colloc_wts[0] = colloc_wts[6] = 0.129484966168870 * wtFactor;
      colloc_wts[1] = colloc_wts[5] = 0.279705391489277 * wtFactor;
      colloc_wts[2] = colloc_wts[4] = 0.381830050505119 * wtFactor;
      colloc_wts[3] = 0.417959183673469 * wtFactor;	break;
    case 8:
      colloc_wts[0] = colloc_wts[7] = 0.101228536290376 * wtFactor;
      colloc_wts[1] = colloc_wts[6] = 0.222381034453374 * wtFactor;
      colloc_wts[2] = colloc_wts[5] = 0.313706645877887 * wtFactor;
      colloc_wts[3] = colloc_wts[4] = 0.362683783378362 * wtFactor;
      break;
    case 9:
      colloc_wts[0] = colloc_wts[8] = 0.081274388361574 * wtFactor;
      colloc_wts[1] = colloc_wts[7] = 0.180648160694857 * wtFactor;
      colloc_wts[2] = colloc_wts[6] = 0.260610696402935 * wtFactor;
      colloc_wts[3] = colloc_wts[5] = 0.312347077040003 * wtFactor;
      colloc_wts[4] = 0.330239355001260 * wtFactor;	break;
    case 10:
      colloc_wts[0] = colloc_wts[9] = 0.066671344308688 * wtFactor;
      colloc_wts[1] = colloc_wts[8] = 0.149451349150581 * wtFactor;
      colloc_wts[2] = colloc_wts[7] = 0.219086362515982 * wtFactor;
      colloc_wts[3] = colloc_wts[6] = 0.269266719309996 * wtFactor;
      colloc_wts[4] = colloc_wts[5] = 0.295524224714753 * wtFactor;
      break;
    default:
      // define Gauss wts from Gauss pts using
      // -(A_{n+1} gamma_n)/(A_n Phi_n'(x_i) Phi_{n+1}(x_i)),
      // which for P(x) with w(x) = 1/2 is (1-x_i^2)/(n P_{n-1}(x_i))^2.
      const RealArray& colloc_pts = collocation_points(order); Real x_i;
      for (size_t i=0; i<order; i++) {
	x_i = colloc_pts[i];
	colloc_wts[i] = (1.-x_i*x_i)
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
    colloc_wts[i] *= wtFactor;
#endif

  if (rule_err) {
    PCerr << "Error: unsupported collocation rule in LegendreOrthogPolynomial::"
	  << "type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }
  return colloc_wts;
}

} // namespace Pecos
