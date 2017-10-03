/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HermiteOrthogPolynomial
//- Description:  Implementation code for HermiteOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "HermiteOrthogPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.hpp"
#endif


namespace Pecos {


Real HermiteOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  // employ Horner's rule for improved efficiency and precision
  Real t1_val;
  switch (order) {
  case 0:
    t1_val = 1.;
    break;
  case 1:
    t1_val = x;
    break;
  case 2:
    t1_val = x*x - 1.;
    break;
  case 3:
    t1_val = x*(x*x - 3.);
    break;
  case 4: {
    Real x2 = x*x;
    t1_val = x2*(x2 - 6.) + 3.;
    break;
  }
  case 5: {
    Real x2 = x*x;
    t1_val = x*(x2*(x2 - 10.) + 15.);
    break;
  }
  case 6: {
    Real x2 = x*x;
    t1_val = x2*(x2*(x2 - 15.) + 45.) - 15.;
    break;
  }
  case 7: {
    Real x2 = x*x;
    t1_val = x*(x2*(x2*(x2 - 21.) + 105.) - 105.);
    break;
  }
  case 8: {
    Real x2 = x*x;
    t1_val = x2*(x2*(x2*(x2 - 28.) + 210.) - 420.) + 105.;
    break;
  }
  case 9: {
    Real x2 = x*x;
    t1_val = x*(x2*(x2*(x2*(x2 - 36.) + 378.) - 1260.) + 945.);
    break;
  }
  case 10: {
    Real x2 = x*x;
    t1_val = x2*(x2*(x2*(x2*(x2 - 45.) + 630.) - 3150.) + 4725.) - 945.;
    break;
  }
  default:
    // Support higher order polynomials using the 3 point recursion formula:
    Real x2 = x*x,
      He_n = x2*(x2*(x2*(x2*(x2 - 45.) + 630.) - 3150.) + 4725.) - 945.,// He_10
      He_nminus1 = x*(x2*(x2*(x2*(x2 - 36.) + 378.) - 1260.) + 945.);   // He_9
    for (size_t i=10; i<order; i++) {
      t1_val = x*He_n - i*He_nminus1; // He_nplus1
      if (i != order-1) {
	He_nminus1 = He_n;
	He_n       = t1_val;
      }
    }
    break;
  }

  return t1_val;
}


Real HermiteOrthogPolynomial::type1_gradient(Real x, unsigned short order)
{ return (order) ? order*type1_value(x, order-1) : 0.; }


Real HermiteOrthogPolynomial::type1_hessian(Real x, unsigned short order)
{ return (order>1) ? order*(order-1)*type1_value(x, order-2) : 0.; }


Real HermiteOrthogPolynomial::norm_squared(unsigned short order)
{ return factorial(order); }


const RealArray& HermiteOrthogPolynomial::
collocation_points(unsigned short order)
{
  // pull this outside block below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "HermiteOrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  bool rule_err = false;
  if (collocPoints.size() != order) { // if not already computed
    collocPoints.resize(order);
#ifdef HAVE_SPARSE_GRID
    if (collocRule == GENZ_KEISTER) {
      webbur::hermite_genz_keister_lookup_points(order, &collocPoints[0]);
      for (size_t i=0; i<order; i++)
	collocPoints[i] *= ptFactor; // scale H_n by sr2 to get He_n
    }
    else if (collocRule == GAUSS_HERMITE) {
      if (order <= 20) { // retrieve full precision tabulated values
	webbur::hermite_lookup_points(order, &collocPoints[0]);
	for (size_t i=0; i<order; i++)
	  collocPoints[i] *= ptFactor; // scale H_n by sr2 to get He_n
      }
      else { // sandia_rules calculates points/weights together
	if (collocWeights.size() != order)
	  collocWeights.resize(order);
	webbur::hermite_compute(order, &collocPoints[0], &collocWeights[0]);
	for (size_t i=0; i<order; i++) {
	  collocPoints[i]  *= ptFactor; // scale H_n by sr2 to get He_n
	  collocWeights[i] *= wtFactor; // polynomial weight fn -> PDF
	}
      }
    }
    else
      rule_err = true;
#else
    if (collocRule == GENZ_KEISTER) {
      PCerr << "Error: VPISparseGrid required for Genz-Keister points in "
	    << "HermiteOrthogPolynomial::collocation_points()." << std::endl;
      abort_handler(-1);
    }
    else if (collocRule == GAUSS_HERMITE) {
      switch (order) {
      case 1: // zeros of He_1(x) for one Gauss-Hermite point:
	collocPoints[0] = 0.0;  break;
      case 2: // zeros of He_2(x) for two Gauss-Hermite points:
	collocPoints[0] = -1.0;
	collocPoints[1] =  1.0; break;
      case 3: { // zeros of He_3(x) for three Gauss-Hermite points:
	Real sr3 = std::sqrt(3.);
	collocPoints[0] = -sr3;
	collocPoints[1] =  0.0;
	collocPoints[2] =  sr3; break;
      }
      case 4: { // zeros of He_4(x) for four Gauss-Hermite points:
	Real sr3 = std::sqrt(3.), sr6 = std::sqrt(6.),
	  sr3psr6 = std::sqrt(3.+sr6), sr3msr6 = std::sqrt(3.-sr6);
	collocPoints[0] = -sr3psr6;
	collocPoints[1] = -sr3msr6;
	collocPoints[2] =  sr3msr6;
	collocPoints[3] =  sr3psr6; break;
      }
      case 5: { // zeros of He_5(x) for five Gauss-Hermite points:
	Real sr10 = std::sqrt(10.), sr5psr10 = std::sqrt(5.+sr10),
	  sr5msr10 = std::sqrt(5.-sr10);
	collocPoints[0] = -sr5psr10;
	collocPoints[1] = -sr5msr10;
	collocPoints[2] =  0.0;
	collocPoints[3] =  sr5msr10;
	collocPoints[4] =  sr5psr10; break;
      }
	// tabulated values from Abramowitz & Stegun have limited precision
      case 6:
	collocPoints[0] = -2.350604973674492 * ptFactor;
	collocPoints[1] = -1.335849074013697 * ptFactor;
	collocPoints[2] = -0.436077411927617 * ptFactor;
	collocPoints[3] = -collocPoints[2];
	collocPoints[4] = -collocPoints[1];
	collocPoints[5] = -collocPoints[0]; break;
      case 7:
	collocPoints[0] = -2.651961356835233 * ptFactor;
	collocPoints[1] = -1.673551628767471 * ptFactor;
	collocPoints[2] = -0.816287882858965 * ptFactor;
	collocPoints[3] =  0.0;
	collocPoints[4] = -collocPoints[2];
	collocPoints[5] = -collocPoints[1];
	collocPoints[6] = -collocPoints[0]; break;
      case 8:
	collocPoints[0] = -2.930637420257244 * ptFactor;
	collocPoints[1] = -1.981656756695843 * ptFactor;
	collocPoints[2] = -1.157193712446780 * ptFactor;
	collocPoints[3] = -0.381186990207322 * ptFactor;
	collocPoints[4] = -collocPoints[3];
	collocPoints[5] = -collocPoints[2];
	collocPoints[6] = -collocPoints[1];
	collocPoints[7] = -collocPoints[0]; break;
      case 9:
	collocPoints[0] = -3.190993201781528 * ptFactor;
	collocPoints[1] = -2.266580584531843 * ptFactor;
	collocPoints[2] = -1.468553289216668 * ptFactor;
	collocPoints[3] = -0.723551018752838 * ptFactor;
	collocPoints[4] =  0.0;
	collocPoints[5] = -collocPoints[3];
	collocPoints[6] = -collocPoints[2];
	collocPoints[7] = -collocPoints[1];
	collocPoints[8] = -collocPoints[0]; break;
      case 10:
	collocPoints[0] = -3.436159118837738 * ptFactor;
	collocPoints[1] = -2.532731674232790 * ptFactor;
	collocPoints[2] = -1.756683649299882 * ptFactor;
	collocPoints[3] = -1.036610829789514 * ptFactor;
	collocPoints[4] = -0.342901327223705 * ptFactor;
	collocPoints[5] = -collocPoints[4];
	collocPoints[6] = -collocPoints[3];
	collocPoints[7] = -collocPoints[2];
	collocPoints[8] = -collocPoints[1];
	collocPoints[9] = -collocPoints[0]; break;
      default:
	PCerr << "Error: overflow in maximum quadrature order limit (10) in "
	      << "HermiteOrthogPolynomial::collocation_points().  Configure "
	      << "with VPISparseGrid to extend range." << std::endl;
	abort_handler(-1); break;
      }
    }
    else
      rule_err = true;
#endif
  }

  if (rule_err) {
    PCerr << "Error: unsupported collocation rule in "
	  << "HermiteOrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  return collocPoints;
}


const RealArray& HermiteOrthogPolynomial::
type1_collocation_weights(unsigned short order)
{
  // pull this outside block below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in Hermite"
	  << "OrthogPolynomial::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  // The sums of the weights = 1, which is the integral of the density function
  // 1/sqrt(2*PI) exp(-x^2/2) over the support range of [-infinity,+infinity]
  // (the std normal CDF for +infinity).

  bool rule_err = false;
  if (collocWeights.size() != order) { // if not already computed
    collocWeights.resize(order);
#ifdef HAVE_SPARSE_GRID
    if (collocRule == GENZ_KEISTER) {
      webbur::hermite_genz_keister_lookup_weights(order, &collocWeights[0]);
      for (size_t i=0; i<order; i++)
	collocWeights[i] *= wtFactor; // polynomial weight fn -> PDF
    }
    else if (collocRule == GAUSS_HERMITE) {
      if (order <= 20) { // retrieve full precision tabulated values
	webbur::hermite_lookup_weights(order, &collocWeights[0]);
	for (size_t i=0; i<order; i++)
	  collocWeights[i] *= wtFactor; // polynomial weight fn -> PDF
      }
      else { // sandia_rules calculates points/weights together
	if (collocPoints.size() != order)
	  collocPoints.resize(order);
	webbur::hermite_compute(order, &collocPoints[0], &collocWeights[0]);
	for (size_t i=0; i<order; i++) {
	  collocPoints[i]  *= ptFactor; // scale H_n pts by sr2 to get He_n pts
	  collocWeights[i] *= wtFactor; // polynomial weight fn -> PDF
	}
      }
    }
    else
      rule_err = true;
#else
    if (collocRule == GENZ_KEISTER) {
      PCerr << "Error: VPISparseGrid required for Genz-Keister points in "
	    << "HermiteOrthogPolynomial::type1_collocation_weights()."
	    << std::endl;
      abort_handler(-1);
    }
    else if (collocRule == GAUSS_HERMITE) {
      switch (order) {
      case 1: // weights for one Gauss-Hermite point:
	collocWeights[0] = 1.0; break;
      case 2: // weights for two Gauss-Hermite points:
	collocWeights[0] = collocWeights[1] = 0.5; break;
      case 3: // weights for three Gauss-Hermite points:
	collocWeights[0] = collocWeights[2] = 1./6.;
	collocWeights[1] = 2./3.; break;
      case 4: { // weights for four Gauss-Hermite points:
	Real sr6 = std::sqrt(6.);
	collocWeights[0] = collocWeights[3] = 1./4./(3.+sr6);
	collocWeights[1] = collocWeights[2] = 1./4./(3.-sr6); break;
      }
      case 5: { // weights for five Gauss-Hermite points:
	Real w2sr10 = 2.*std::sqrt(10.);
	collocWeights[0] = collocWeights[4] = 3./20./(7.+w2sr10);
	collocWeights[1] = collocWeights[3] = 3./20./(7.-w2sr10);
	collocWeights[2] = 8./15.; break;
      }
	// tabulated values from Abramowitz & Stegun have limited precision
      case 6:
	collocWeights[0] = collocWeights[5] = 4.530009905509e-3 * wtFactor;
	collocWeights[1] = collocWeights[4] = 0.1570673203229 * wtFactor;
	collocWeights[2] = collocWeights[3] = 0.7246295952244 * wtFactor; break;
      case 7:
	collocWeights[0] = collocWeights[6] = 9.717812450995e-4 * wtFactor;
	collocWeights[1] = collocWeights[5] = 5.451558281913e-2 * wtFactor;
	collocWeights[2] = collocWeights[4] = 0.4256072526101 * wtFactor;
	collocWeights[3] = 0.8102646175568 * wtFactor; break;
      case 8:
	collocWeights[0] = collocWeights[7] = 1.996040722114e-4 * wtFactor;
	collocWeights[1] = collocWeights[6] = 1.707798300741e-2 * wtFactor;
	collocWeights[2] = collocWeights[5] = 0.2078023258149 * wtFactor;
	collocWeights[3] = collocWeights[4] = 0.6611470125582 * wtFactor; break;
      case 9:
	collocWeights[0] = collocWeights[8] = 3.960697726326e-5 * wtFactor;
	collocWeights[1] = collocWeights[7] = 4.943624275537e-3 * wtFactor;
	collocWeights[2] = collocWeights[6] = 8.847452739438e-2 * wtFactor;
	collocWeights[3] = collocWeights[5] = 0.4326515590026 * wtFactor;
	collocWeights[4] = 0.7202352156061 * wtFactor; break;
      case 10:
	collocWeights[0] = collocWeights[9] = 7.640432855233e-6 * wtFactor;
	collocWeights[1] = collocWeights[8] = 1.343645746781e-3 * wtFactor;
	collocWeights[2] = collocWeights[7] = 3.387439445548e-2 * wtFactor;
	collocWeights[3] = collocWeights[6] = 0.2401386110823 * wtFactor;
	collocWeights[4] = collocWeights[5] = 0.6108626337353 * wtFactor; break;
      default:
	// define Gauss wts from Gauss pts using
	// -(A_{n+1} gamma_n)/(A_n Phi_n'(x_i) Phi_{n+1}(x_i)),
	// which for He(x), is n!/(n He_{n-1}(x_i))^2.
	const RealArray& colloc_pts = collocation_points(order);
	for (size_t i=0; i<order; i++)
	  collocWeights[i] = factorial(order)/std::pow(order*
	    type1_value(colloc_pts[i], order-1), 2);
	break;
      }
    }
    else
      rule_err = true;
#endif
  }

  if (rule_err) {
    PCerr << "Error: unsupported collocation rule in HermiteOrthogPolynomial::"
	  << "type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  return collocWeights;
}

} // namespace Pecos
