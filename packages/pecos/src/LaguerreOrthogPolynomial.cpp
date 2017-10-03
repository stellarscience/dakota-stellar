/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LaguerreOrthogPolynomial
//- Description:  Implementation code for LaguerreOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "LaguerreOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.hpp"
#endif

//#define DEBUG


namespace Pecos {

Real LaguerreOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  // employ Horner's rule for improved efficiency and precision
  Real t1_val;
  switch (order) {
  case 0:
    t1_val = 1.;                                                          break;
  case 1:
    t1_val = -x + 1.;                                                     break;
  case 2:
    t1_val = (x*(x - 4.) + 2.) / 2.;                                      break;
  case 3:
    t1_val = (x*(x*(-x + 9.) - 18.) + 6.) / 6.;                           break;
  case 4:
    t1_val = (x*(x*(x*(x - 16.) + 72.) - 96.) + 24.) / 24.;               break;
  case 5:
    t1_val = (x*(x*(x*(x*(-x + 25.) - 200.) + 600.) - 600.) + 120.) / 120.;
    break;
  case 6:
    t1_val = (x*(x*(x*(x*(x*(x - 36.) + 450.) - 2400.) + 5400.) - 4320.) + 720.)
           / 720.;
    break;
  case 7:
    t1_val = (x*(x*(x*(x*(x*(x*(-x + 49.) - 882.) + 7350.) - 29400.) + 52920.)
	   - 35280.) + 5040.) / 5040.;
    break;
  case 8:
    t1_val = (x*(x*(x*(x*(x*(x*(x*(x - 64.) + 1568.) - 18816.) + 117600.)
	   - 376320.) + 564480.) - 322560.) + 40320.) / 40320.;
    break;
  case 9:
    t1_val = (x*(x*(x*(x*(x*(x*(x*(x*(-x + 81.) - 2592.) + 42336.) - 381024.)
	   + 1905120.) - 5080320.) + 6531840.) - 3265920.) + 362880.) / 362880.;
    break;
  case 10:
    t1_val = (x*(x*(x*(x*(x*(x*(x*(x*(x*(x - 100.) + 4050.) - 86400.)
	   + 1058400.) - 7620480.) + 31752000.) - 72576000.) + 81648000.)
	   - 36288000.) + 3628800.) / 3628800.;
    break;
  default:
    // Support higher order polynomials using the 3 point recursion formula
    Real L_n = (x*(x*(x*(x*(x*(x*(x*(x*(x*(x - 100.) + 4050.) - 86400.)
	     + 1058400.) - 7620480.) + 31752000.) - 72576000.) + 81648000.)
	     - 36288000.) + 3628800.) / 3628800., // L_10
      L_nminus1 = (x*(x*(x*(x*(x*(x*(x*(x*(-x + 81.) - 2592.) + 42336.)
		- 381024.) + 1905120.) - 5080320.) + 6531840.) - 3265920.)
		+ 362880.) / 362880.;             // L_9
    for (size_t i=10; i<order; i++) {
      t1_val = ( (2.*i+1.-x)*L_n - i*L_nminus1 ) / (i+1.); // L_nplus1
      if (i != order-1) {
	L_nminus1 = L_n;
	L_n       = t1_val;
      }
    }
    break;
  }

  return t1_val;
}


Real LaguerreOrthogPolynomial::
type1_gradient(Real x, unsigned short order)
{ 
  Real t1_grad;
#ifdef DEBUG
  // See Abramowitz & Stegun, Section 22.8, p.783
  //t1_grad = (order) ?
  //  order*(type1_value(x, order) - type1_value(x, order-1))/x : 0.;
  if (order) {
    Real L_n = type1_value(x, order), L_nminus1 = type1_value(x, order-1);
    t1_grad = order*(L_n - L_nminus1)/x;
  }
  else
    t1_grad = 0.;
  PCout << "Laguerre gradient approach 1: " << t1_grad << '\n';
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
    t1_grad = x - 2.;                                                     break;
  case 3:
    t1_grad = (x*(-x + 6.) - 6.)/2.;                                      break;
  case 4:
    t1_grad = (x*(x*(x - 12.) + 36.) - 24.)/6.;                           break;
  case 5:
    t1_grad = (x*(x*(x*(-x + 20.) - 120.) + 240.) - 120.) / 24.;          break;
  case 6:
    t1_grad = (x*(x*(x*(x*(x - 30.) + 300.) - 1200.) + 1800.) - 720.) / 120.;
    break;
  default:
    // Support higher order polynomials using the 3 point recursion formula
    Real dLdx_n
      = (x*(x*(x*(x*(x - 30.) + 300.) - 1200.) + 1800.) - 720.)/120., // L'_6
      dLdx_nminus1 = (x*(x*(x*(-x + 20.) - 120.) + 240.) - 120.)/24.; // L'_5
    for (size_t i=6; i<order; i++) {
      t1_grad // dLdx_nplus1
	= ( (2.*i+1.-x)*dLdx_n - type1_value(x,i) - i*dLdx_nminus1 ) / (i+1.);
      if (i != order-1) {
	dLdx_nminus1 = dLdx_n;
	dLdx_n       = t1_grad;
      }
    }
    break;
  }
#ifdef DEBUG
  PCout << "Laguerre gradient approach 2: " << t1_grad << '\n';
#endif // DEBUG

  return t1_grad;
}


Real LaguerreOrthogPolynomial::type1_hessian(Real x, unsigned short order)
{
  Real t1_hess;
  switch (order) {
  case 0: case 1:
    t1_hess = 0.;                                                         break;
  case 2:
    t1_hess = 1.;                                                         break;
  case 3:
    t1_hess = 3. - x;                                                     break;
  case 4:
    t1_hess = (x*(x - 8.) + 12.)/2.;                                      break;
  case 5:
    t1_hess = (x*(x*(-x + 15.) - 60.) + 60.) / 6.;          break;
  case 6:
    t1_hess = (x*(x*(x*(x - 24.) + 180.) - 480.) + 360.) / 24.;
    break;
  default:
    // Support higher order polynomials using the 3 point recursion formula
    Real d2Ldx2_n = (x*(x*(x*(x - 24.) + 180.) - 480.) + 360.) / 24., // L''_6
      d2Ldx2_nminus1 = (x*(x*(-x + 15.) - 60.) + 60.) / 6.;           // L''_5
    for (size_t i=6; i<order; i++) {
      t1_hess = ( (2.*i+1.-x)*d2Ldx2_n - 2.*type1_gradient(x,i) -
		  i*d2Ldx2_nminus1 ) / (i+1.); // d2Ldx2_nplus1
      if (i != order-1) {
	d2Ldx2_nminus1 = d2Ldx2_n;
	d2Ldx2_n       = t1_hess;
      }
    }
    break;
  }

  return t1_hess;
}


Real LaguerreOrthogPolynomial::norm_squared(unsigned short order)
{ return 1.; }


const RealArray& LaguerreOrthogPolynomial::
collocation_points(unsigned short order)
{
  // pull this outside block below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "LaguerreOrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  if (collocPoints.size() != order) { // if not already computed
    collocPoints.resize(order);
#ifdef HAVE_SPARSE_GRID
    if (order <= 20) // retrieve full precision tabulated values
      webbur::laguerre_lookup_points(order, &collocPoints[0]);
    else { // calculates points/weights together
      if (collocWeights.size() != order)
	collocWeights.resize(order);
      webbur::laguerre_compute(order, &collocPoints[0], &collocWeights[0]);
    }
#else
    switch (order) {
    case 1: // zeros of L_1(x) for one Gauss-Laguerre point:
      collocPoints[0] =  1.0; break;
    case 2: { // zeros of L_2(x) for two Gauss-Laguerre points:
      Real sr2 = std::sqrt(2.);
      collocPoints[0] =  2. - sr2;
      collocPoints[1] =  2. + sr2; break;
    }
    // Only ~12 digits of precision in Abramowitz & Stegun tabulated values
    case 3:
      collocPoints[0] =  0.415774556783;
      collocPoints[1] =  2.294280360279;
      collocPoints[2] =  6.289945082937; break;
    case 4:
      collocPoints[0] =  0.322547689619;
      collocPoints[1] =  1.745761101158;
      collocPoints[2] =  4.536620296921;
      collocPoints[3] =  9.395070912301; break;
    case 5:
      collocPoints[0] =  0.263560319718;
      collocPoints[1] =  1.413403059107;
      collocPoints[2] =  3.596425771041;
      collocPoints[3] =  7.085810005859;
      collocPoints[4] = 12.640800844276; break;
    case 6:
      collocPoints[0] =  0.222846604179;
      collocPoints[1] =  1.188932101673;
      collocPoints[2] =  2.992736326059;
      collocPoints[3] =  5.775143569105;
      collocPoints[4] =  9.837467418383;
      collocPoints[5] = 15.982873980602; break;
    case 7:
      collocPoints[0] =  0.193043676560;
      collocPoints[1] =  1.026664895339;
      collocPoints[2] =  2.567876744951;
      collocPoints[3] =  4.900353084526;
      collocPoints[4] =  8.182153444563;
      collocPoints[5] = 12.734180291798;
      collocPoints[6] = 19.395727862263; break;
    case 8:
      collocPoints[0] =  0.170279632305;
      collocPoints[1] =  0.903701776799;
      collocPoints[2] =  2.251086629866;
      collocPoints[3] =  4.266700170288;
      collocPoints[4] =  7.045905402393;
      collocPoints[5] = 10.758516010181;
      collocPoints[6] = 15.740678641278;
      collocPoints[7] = 22.863131736889; break;
    case 9:
      collocPoints[0] =  0.152322227732;
      collocPoints[1] =  0.807220022742;
      collocPoints[2] =  2.005135155619;
      collocPoints[3] =  3.783473973331;
      collocPoints[4] =  6.204956777877;
      collocPoints[5] =  9.372985251688;
      collocPoints[6] = 13.466236911092;
      collocPoints[7] = 18.833597788992;
      collocPoints[8] = 26.374071890927; break;
    case 10:
      collocPoints[0] =  0.137793470540;
      collocPoints[1] =  0.729454549503;
      collocPoints[2] =  1.808342901740;
      collocPoints[3] =  3.401433697855;
      collocPoints[4] =  5.552496140064;
      collocPoints[5] =  8.330152746764;
      collocPoints[6] = 11.843785837900;
      collocPoints[7] = 16.279257831378;
      collocPoints[8] = 21.996585811981;
      collocPoints[9] = 29.920697012274; break;
    default:
      PCerr << "Error: overflow in maximum quadrature order limit (10) in "
	    << "LaguerreOrthogPolynomial::collocation_points().  Configure "
	    << "with VPISparseGrid to extend range." << std::endl;
      abort_handler(-1); break;
    }
#endif
  }

  return collocPoints;
}


const RealArray& LaguerreOrthogPolynomial::
type1_collocation_weights(unsigned short order)
{
  // pull this outside block below since order=0 is initial colloc wts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in Laguerre"
	  << "OrthogPolynomial::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  // The sums of the weights = 1, which is the integral of the density function
  // exp(-x) over the support range of [0,+infinity].

  if (collocWeights.size() != order) { // if not already computed
    collocWeights.resize(order);
#ifdef HAVE_SPARSE_GRID
    if (order <= 20) // tabulated values from sandia_rules have full precision
      webbur::laguerre_lookup_weights(order, &collocWeights[0]);
    else { // sandia_rules calculates points/weights together
      if (collocPoints.size() != order)
	collocPoints.resize(order);
      webbur::laguerre_compute(order, &collocPoints[0], &collocWeights[0]);
    }
#else
    switch (order) {
    case 1: // weights for one Gauss-Laguerre point:
      collocWeights[0] = 1.0; break;
    case 2: { // weights for two Gauss-Laguerre points:
      Real sr2 = std::sqrt(2.);
      collocWeights[0] = (2. + sr2)/4.;
      collocWeights[1] = (2. - sr2)/4.; break;
    }
    // Only ~12 digits of precision in Abramowitz & Stegun tabulated values
    case 3:
      collocWeights[0] = 0.711093009929;
      collocWeights[1] = 0.278517733569;
      collocWeights[2] = 0.0103892565016; break;
    case 4:
      collocWeights[0] = 0.603154104342;
      collocWeights[1] = 0.357418692438;
      collocWeights[2] = 0.0388879085150;
      collocWeights[3] = 0.000539294705561; break;
    case 5:
      collocWeights[0] = 0.521755610583;
      collocWeights[1] = 0.398666811083;
      collocWeights[2] = 0.0759424496817;
      collocWeights[3] = 0.00361175867992;
      collocWeights[4] = 2.33699723858e-5; break;
    case 6:
      collocWeights[0] = 0.458964673950;
      collocWeights[1] = 0.417000830772;
      collocWeights[2] = 0.113373382074;
      collocWeights[3] = 0.0103991974531;
      collocWeights[4] = 0.000261017202815;
      collocWeights[5] = 8.98547906430e-7; break;
    case 7:
      collocWeights[0] = 0.409318951701;
      collocWeights[1] = 0.421831277862;
      collocWeights[2] = 0.147126348658;
      collocWeights[3] = 0.0206335144687;
      collocWeights[4] = 0.00107401014328;
      collocWeights[5] = 1.58654643486e-5;
      collocWeights[6] = 3.17031547900e-8; break;
    case 8:
      collocWeights[0] = 0.369188589342;
      collocWeights[1] = 0.418786780814;
      collocWeights[2] = 0.175794986637;
      collocWeights[3] = 0.0333434922612;
      collocWeights[4] = 0.00279453623523;
      collocWeights[5] = 9.07650877336e-5;
      collocWeights[6] = 8.48574671627e-7;
      collocWeights[7] = 1.04800117487e-9; break;
    case 9:
      collocWeights[0] = 0.336126421798;
      collocWeights[1] = 0.411213980424;
      collocWeights[2] = 0.199287525371;
      collocWeights[3] = 0.0474605627657;
      collocWeights[4] = 0.00559962661079;
      collocWeights[5] = 0.000305249767093;
      collocWeights[6] = 6.59212302608e-6;
      collocWeights[7] = 4.11076933035e-8;
      collocWeights[8] = 3.29087403035e-11; break;
    case 10:
      collocWeights[0] = 0.308441115765;
      collocWeights[1] = 0.401119929155;
      collocWeights[2] = 0.218068287612;
      collocWeights[3] = 0.0620874560987;
      collocWeights[4] = 0.00950151697518;
      collocWeights[5] = 0.000753008388588;
      collocWeights[6] = 2.82592334960e-5;
      collocWeights[7] = 4.24931398496e-7;
      collocWeights[8] = 1.83956482398e-9;
      collocWeights[9] = 9.91182721961e-13; break;
    default:
      // define Gauss wts from Gauss pts using
      // -(A_{n+1} gamma_n)/(A_n Phi_n'(x_i) Phi_{n+1}(x_i)),
      // which for L(x), is x_i/(n L_{n-1}(x_i))^2.
      const RealArray& colloc_pts = collocation_points(order);
      for (size_t i=0; i<order; i++)
	collocWeights[i] = colloc_pts[i] / std::pow(order *
			   type1_value(colloc_pts[i], order-1), 2);
      break;
    }
#endif
  }

  return collocWeights;
}

} // namespace Pecos
