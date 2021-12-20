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

  UShortRealArrayMap::iterator it = collocPointsMap.find(order);
  if (it != collocPointsMap.end())
    return it->second;
    
  RealArray& colloc_pts = collocPointsMap[order]; // create new array
  colloc_pts.resize(order);
#ifdef HAVE_SPARSE_GRID
  if (order <= 20) // retrieve full precision tabulated values
    webbur::laguerre_lookup_points(order, &colloc_pts[0]);
  else { // calculates points/weights together
    RealArray& colloc_wts = collocWeightsMap[order];
    if (colloc_wts.size() != order)
      colloc_wts.resize(order);
    webbur::laguerre_compute(order, &colloc_pts[0], &colloc_wts[0]);
  }
#else
  switch (order) {
  case 1: // zeros of L_1(x) for one Gauss-Laguerre point:
    colloc_pts[0] =  1.0; break;
  case 2: { // zeros of L_2(x) for two Gauss-Laguerre points:
    Real sr2 = std::sqrt(2.);
    colloc_pts[0] =  2. - sr2;
    colloc_pts[1] =  2. + sr2; break;
  }
  // Only ~12 digits of precision in Abramowitz & Stegun tabulated values
  case 3:
    colloc_pts[0] =  0.415774556783;
    colloc_pts[1] =  2.294280360279;
    colloc_pts[2] =  6.289945082937; break;
  case 4:
    colloc_pts[0] =  0.322547689619;
    colloc_pts[1] =  1.745761101158;
    colloc_pts[2] =  4.536620296921;
    colloc_pts[3] =  9.395070912301; break;
  case 5:
    colloc_pts[0] =  0.263560319718;
    colloc_pts[1] =  1.413403059107;
    colloc_pts[2] =  3.596425771041;
    colloc_pts[3] =  7.085810005859;
    colloc_pts[4] = 12.640800844276; break;
  case 6:
    colloc_pts[0] =  0.222846604179;
    colloc_pts[1] =  1.188932101673;
    colloc_pts[2] =  2.992736326059;
    colloc_pts[3] =  5.775143569105;
    colloc_pts[4] =  9.837467418383;
    colloc_pts[5] = 15.982873980602; break;
  case 7:
    colloc_pts[0] =  0.193043676560;
    colloc_pts[1] =  1.026664895339;
    colloc_pts[2] =  2.567876744951;
    colloc_pts[3] =  4.900353084526;
    colloc_pts[4] =  8.182153444563;
    colloc_pts[5] = 12.734180291798;
    colloc_pts[6] = 19.395727862263; break;
  case 8:
    colloc_pts[0] =  0.170279632305;
    colloc_pts[1] =  0.903701776799;
    colloc_pts[2] =  2.251086629866;
    colloc_pts[3] =  4.266700170288;
    colloc_pts[4] =  7.045905402393;
    colloc_pts[5] = 10.758516010181;
    colloc_pts[6] = 15.740678641278;
    colloc_pts[7] = 22.863131736889; break;
  case 9:
    colloc_pts[0] =  0.152322227732;
    colloc_pts[1] =  0.807220022742;
    colloc_pts[2] =  2.005135155619;
    colloc_pts[3] =  3.783473973331;
    colloc_pts[4] =  6.204956777877;
    colloc_pts[5] =  9.372985251688;
    colloc_pts[6] = 13.466236911092;
    colloc_pts[7] = 18.833597788992;
    colloc_pts[8] = 26.374071890927; break;
  case 10:
    colloc_pts[0] =  0.137793470540;
    colloc_pts[1] =  0.729454549503;
    colloc_pts[2] =  1.808342901740;
    colloc_pts[3] =  3.401433697855;
    colloc_pts[4] =  5.552496140064;
    colloc_pts[5] =  8.330152746764;
    colloc_pts[6] = 11.843785837900;
    colloc_pts[7] = 16.279257831378;
    colloc_pts[8] = 21.996585811981;
    colloc_pts[9] = 29.920697012274; break;
  default:
    PCerr << "Error: overflow in maximum quadrature order limit (10) in "
	  << "LaguerreOrthogPolynomial::collocation_points().  Configure "
	  << "with VPISparseGrid to extend range." << std::endl;
    abort_handler(-1); break;
  }
#endif

  return colloc_pts;
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

  UShortRealArrayMap::iterator it = collocWeightsMap.find(order);
  if (it != collocWeightsMap.end())
    return it->second;

  RealArray& colloc_wts = collocWeightsMap[order]; // create new array
  colloc_wts.resize(order);
#ifdef HAVE_SPARSE_GRID
  if (order <= 20) // tabulated values from sandia_rules have full precision
    webbur::laguerre_lookup_weights(order, &colloc_wts[0]);
  else { // sandia_rules calculates points/weights together
    RealArray& colloc_pts = collocPointsMap[order];
    if (colloc_pts.size() != order)
      colloc_pts.resize(order);
    webbur::laguerre_compute(order, &colloc_pts[0], &colloc_wts[0]);
  }
#else
  switch (order) {
  case 1: // weights for one Gauss-Laguerre point:
    colloc_wts[0] = 1.0; break;
  case 2: { // weights for two Gauss-Laguerre points:
    Real sr2 = std::sqrt(2.);
    colloc_wts[0] = (2. + sr2)/4.;
    colloc_wts[1] = (2. - sr2)/4.; break;
  }
  // Only ~12 digits of precision in Abramowitz & Stegun tabulated values
  case 3:
    colloc_wts[0] = 0.711093009929;
    colloc_wts[1] = 0.278517733569;
    colloc_wts[2] = 0.0103892565016; break;
  case 4:
    colloc_wts[0] = 0.603154104342;
    colloc_wts[1] = 0.357418692438;
    colloc_wts[2] = 0.0388879085150;
    colloc_wts[3] = 0.000539294705561; break;
  case 5:
    colloc_wts[0] = 0.521755610583;
    colloc_wts[1] = 0.398666811083;
    colloc_wts[2] = 0.0759424496817;
    colloc_wts[3] = 0.00361175867992;
    colloc_wts[4] = 2.33699723858e-5; break;
  case 6:
    colloc_wts[0] = 0.458964673950;
    colloc_wts[1] = 0.417000830772;
    colloc_wts[2] = 0.113373382074;
    colloc_wts[3] = 0.0103991974531;
    colloc_wts[4] = 0.000261017202815;
    colloc_wts[5] = 8.98547906430e-7; break;
  case 7:
    colloc_wts[0] = 0.409318951701;
    colloc_wts[1] = 0.421831277862;
    colloc_wts[2] = 0.147126348658;
    colloc_wts[3] = 0.0206335144687;
    colloc_wts[4] = 0.00107401014328;
    colloc_wts[5] = 1.58654643486e-5;
    colloc_wts[6] = 3.17031547900e-8; break;
  case 8:
    colloc_wts[0] = 0.369188589342;
    colloc_wts[1] = 0.418786780814;
    colloc_wts[2] = 0.175794986637;
    colloc_wts[3] = 0.0333434922612;
    colloc_wts[4] = 0.00279453623523;
    colloc_wts[5] = 9.07650877336e-5;
    colloc_wts[6] = 8.48574671627e-7;
    colloc_wts[7] = 1.04800117487e-9; break;
  case 9:
    colloc_wts[0] = 0.336126421798;
    colloc_wts[1] = 0.411213980424;
    colloc_wts[2] = 0.199287525371;
    colloc_wts[3] = 0.0474605627657;
    colloc_wts[4] = 0.00559962661079;
    colloc_wts[5] = 0.000305249767093;
    colloc_wts[6] = 6.59212302608e-6;
    colloc_wts[7] = 4.11076933035e-8;
    colloc_wts[8] = 3.29087403035e-11; break;
  case 10:
    colloc_wts[0] = 0.308441115765;
    colloc_wts[1] = 0.401119929155;
    colloc_wts[2] = 0.218068287612;
    colloc_wts[3] = 0.0620874560987;
    colloc_wts[4] = 0.00950151697518;
    colloc_wts[5] = 0.000753008388588;
    colloc_wts[6] = 2.82592334960e-5;
    colloc_wts[7] = 4.24931398496e-7;
    colloc_wts[8] = 1.83956482398e-9;
    colloc_wts[9] = 9.91182721961e-13; break;
  default:
    // define Gauss wts from Gauss pts using
    // -(A_{n+1} gamma_n)/(A_n Phi_n'(x_i) Phi_{n+1}(x_i)),
    // which for L(x), is x_i/(n L_{n-1}(x_i))^2.
    const RealArray& colloc_pts = collocation_points(order);
    for (size_t i=0; i<order; i++)
      colloc_wts[i] = colloc_pts[i]
	            / std::pow(order * type1_value(colloc_pts[i], order-1), 2);
    break;
  }
#endif

  return colloc_wts;
}

} // namespace Pecos
