/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "CharlierOrthogPolynomial.hpp"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace Pecos {

Real CharlierOrthogPolynomial::type1_value( Real x, unsigned short order )
{
  Real result;
  switch ( order ) {
  case 0:
    result = 1.; break;

  case 1:
    result = 1. - x/lambdaStat; break;

  case 2:
    result = 1. + x*(x-1.-2.*lambdaStat)/(lambdaStat*lambdaStat); break;

  case 3: {
    Real lam2 = lambdaStat*lambdaStat, lam3 = lambdaStat*lam2;
    result = 1. + x*(-3.*lam2+(2.+3.*lambdaStat-x)*(-1.+x))/lam3;
    break;
  }

  // This version is not consistent with recursion below (breaks the unit test)
  //case 4: {
  //  Real lam2 = lambdaStat*lambdaStat, lam3 = lambdaStat*lam2, x2 = x*x;
  //  result = (lam3*(2.+lambdaStat)-2.*(3.+2.*lambdaStat)*
  //            (1.+lambdaStat+lam2)*x+(11.+2.*lambdaStat*(7.+3.*lambdaStat))
  //            *x2-2.*(3.+2.*lambdaStat)*x2*x + x2*x2)/(lam2*lam2);
  //  break;
  //}

  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real fm2 = type1_value(x, order-2);
    Real fm1 = type1_value(x, order-1);
    result = ( (Real(order)-1.0+lambdaStat-x)*fm1 - (Real(order)-1.0)*fm2 ) / lambdaStat;
    // Unroll recursion for performance
    //Real lam2 = lambdaStat*lambdaStat, lam3 = lambdaStat*lam2,x2 = x*x,
    //  Ch_nm1 = 1.+x*(-3.*lam2+(2.+3.*lambdaStat-x)*(-1.+x))/lam3,  //3
    //  Ch_n = (lam3*(2.+lambdaStat)-2.*(3.+2.*lambdaStat)*
    //          (1.+lambdaStat+lam2)*x+(11.+2.*lambdaStat*(7.+3.*lambdaStat))
    //          *x2-2.*(3.+2.*lambdaStat)*x2*x + x2*x2)/(lam2*lam2); //4
    //for (size_t i=5; i<=order; i++) {
    //  Real om1 = (Real)i - 1.;
    //  result = ((om1+lambdaStat-x)*Ch_n-om1*Ch_nm1)/lambdaStat; // Ch_nplus1
    //  if (i < order)
    //    { Ch_nm1 = Ch_n;  Ch_n = result; }
    //}
    break;
    // The recusion above produces a different result to the following recursion
    // Ch_n(x,a)=1/a*x*Ch_n(x-1,a) - Ch_n(x,a)
    // Specifically every odd polynomial is different by a factor of -1.
  }
  }
  return result;
}


Real CharlierOrthogPolynomial::type1_gradient( Real x, unsigned short order )
{
  Real result;
  switch ( order ) {
  case 0:
    result = 0.; break;

  case 1:
    result = -1./lambdaStat; break;

  case 2:
    result = (2.*(-lambdaStat+x)-1.)/(lambdaStat*lambdaStat); break;

  case 3:
    result = (-2.+(6.-3.*x)*x+lambdaStat*(-3.-3.*lambdaStat+6.*x))
           / (lambdaStat*lambdaStat*lambdaStat);
    break;

  case 4: {
    Real lam2 = lambdaStat*lambdaStat;
    result = (-6.+lambdaStat*(-10.+(-10.-4.*lambdaStat)*lambdaStat)+
	      x*(22.+lambdaStat*(28.+12.*lambdaStat)+
		 x*(-18.-12.*lambdaStat+4.*x)))/(lam2*lam2);
    break;
  }

  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real lam2 = lambdaStat*lambdaStat;
    Real dChdx_nm1 = (-2.+(6.-3.*x)*x+lambdaStat*(-3.-3.*lambdaStat+6.*x))
                   / (lambdaStat*lam2),
         dChdx_n   = (-6.+lambdaStat*(-10.+(-10.-4.*lambdaStat)*lambdaStat)+
		      x*(22.+lambdaStat*(28.+12.*lambdaStat)+
			 x*(-18.-12.*lambdaStat+4.*x)))/(lam2*lam2);
    for ( size_t i=5; i<=order; i++ ) {
      Real om1 = (Real)i - 1.;
      result = ((om1+lambdaStat-x)*dChdx_n-type1_value(x,order)-om1*dChdx_nm1)
	     / lambdaStat;
      if (i < order)
	{ dChdx_nm1 = dChdx_n;	dChdx_n = result; }
    }
    break;
  }
  }
  return result;
};


Real CharlierOrthogPolynomial::type1_hessian( Real x, unsigned short order )
{
  Real result;
  switch ( order ) {
  case 0: case 1:
    result = 0.; break;

  case 2:
    result = 2./(lambdaStat*lambdaStat); break;

  case 3:
    result = 6.*(lambdaStat-x+1.)/(lambdaStat*lambdaStat*lambdaStat); break;

  case 4: {
    Real lam2 = lambdaStat*lambdaStat;
    result
      = (2.*(11.+6.*lam2+2.*lambdaStat*(7.-6.*x)+6.*(-3.+x)*x))/(lam2*lam2);
    break;
  }

  default:{
    // Support higher order polynomials using the 3 point recursion formula:
    Real lam2 = lambdaStat*lambdaStat;
    Real d2Chdx2_nm1 = 6.*(lambdaStat-x+1.)/(lam2*lambdaStat),
         d2Chdx2_n   = (2.*(11.+6.*lam2+2.*lambdaStat*(7.-6.*x)+6.*(-3.+x)*x))
                     / (lam2*lam2);
    for ( size_t i=5; i<=order; i++ ) {
      Real om1 = (Real)i - 1.;
      result = ((om1+lambdaStat-x)*d2Chdx2_n-2.*type1_gradient(x,order)
		-om1*d2Chdx2_nm1)/lambdaStat;
      if (i < order)
	{ d2Chdx2_nm1 = d2Chdx2_n;  d2Chdx2_n = result; }
    }
    break;
  }
  }
  return result;
};


Real CharlierOrthogPolynomial::norm_squared( unsigned short order )
{
  return std::pow( lambdaStat, order ) * boost::math::factorial<Real>( order );
}

} // namespace Pecos
