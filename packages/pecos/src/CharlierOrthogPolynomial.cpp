#include "CharlierOrthogPolynomial.hpp"
#include "MathTools.hpp"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace Pecos {

Real CharlierOrthogPolynomial::type1_value( Real x, unsigned short order )
{
  Real result = 1.;
  switch ( order ) {
  case 0:{
    result = 1.;
    break;
  }
  case 1:{
    result = (alphaPoly-x)/alphaPoly;
    break;
  }
  case 2:{
    Real alpha2 = alphaPoly*alphaPoly;
    result = (alpha2+x*(-2.*alphaPoly+x-1.))/alpha2;
    break;
  }
  case 3:{
    Real alpha2 = alphaPoly*alphaPoly, alpha3 = alphaPoly*alpha2;
    result=(alpha3+(-3.*alpha2+(2.+3.*alphaPoly-x)*(-1.+x))*x)/alpha3;
    break;
  }
  case 4:{
    Real alpha2 = alphaPoly*alphaPoly, alpha3 = alphaPoly*alpha2, x2 = x*x;
    result = (alpha3*(2.+alphaPoly)-2.*(3.+2.*alphaPoly)*(1.+alphaPoly+alpha2)*x+(11.+2.*alphaPoly*(7.+3.*alphaPoly))*x2-2.*(3.+2.*alphaPoly)*x2*x + x2*x2)/(alpha2*alpha2);
    break;
  }
  default:{
    // Support higher order polynomials using the 3 point recursion formula:
    Real Ch_nminus1 = type1_value(x,3), //Ch_3
      Ch_n =  type1_value(x,4); //Ch_4
    for ( size_t i=4; i<order; i++ ) {
      result = ((i+alphaPoly-x)*Ch_n-i*Ch_nminus1)/alphaPoly; // Ch_nplus1
      if (i != order-1) {
	Ch_nminus1 = Ch_n;
	Ch_n       = result;
      }
    }
    break;
    // This above recusion produces a different result to the following recursion
    // Ch_n(x,a)=1/a*x*Ch_n(x-1,a) - Ch_n(x,a)
    // Specifically every odd polynomial is different by a factor of -1.
  }
  }
  return result;
}

Real CharlierOrthogPolynomial::type1_gradient( Real x, unsigned short order )
{
  Real result = 0.;
  switch ( order ) {
  case 0:{
    result = 0.;
    break;
  }
  case 1:{
    result = -1/alphaPoly;
    break;
  }
  case 2:{
    Real alpha2 = alphaPoly*alphaPoly;
    result = (2.*(-alphaPoly+x)-1.)/alpha2;
    break;
  }
  case 3:{
    Real alpha3 = alphaPoly*alphaPoly*alphaPoly;
    result=(-2.+(6.-3.*x)*x+alphaPoly*(-3.-3.*alphaPoly+6.*x))/alpha3;
    break;
  }
  case 4:{
    Real alpha2 = alphaPoly*alphaPoly;
    result = (-6.+alphaPoly*(-10.+(-10.-4.*alphaPoly)*alphaPoly)+x*(22.+alphaPoly*(28.+12.*alphaPoly)+x*(-18.-12.*alphaPoly+4.*x)))/(alpha2*alpha2);
    break;
  }
  default:{
    // Support higher order polynomials using the 3 point recursion formula:
    Real dChdx_nminus1 = type1_gradient(x,3), dChdx_n =  type1_gradient(x,4);
    for ( size_t i=4; i<order; i++ ) {
      result = ((i+alphaPoly-x)*dChdx_n-type1_value(x,order)-i*dChdx_nminus1)/alphaPoly;
      if (i != order-1) {
	dChdx_nminus1 = dChdx_n;
	dChdx_n       = result;
      }
    }
    break;
  }
  }
  return result;
};

Real CharlierOrthogPolynomial::type1_hessian( Real x, unsigned short order )
{
  Real result = 0.;
  switch ( order ) {
  case 0:{
    result = 0.;
    break;
  }
  case 1:{
    result = 0.;
    break;
  }
  case 2:{
    Real alpha2 = alphaPoly*alphaPoly;
    result = 2./alpha2;
    break;
  }
  case 3:{
    Real alpha3 = alphaPoly*alphaPoly*alphaPoly;
    result=6.*(alphaPoly-x+1.)/alpha3;
    break;
  }
  case 4:{
    Real alpha2 = alphaPoly*alphaPoly;
    result = (2.*(11.+6.*alpha2+2.*alphaPoly*(7.-6.*x)+6.*(-3.+x)*x))/alpha2*alpha2;
    break;
  }
  default:{
    // Support higher order polynomials using the 3 point recursion formula:
    Real d2Chdx2_nminus1 = type1_hessian(x,3), d2Chdx2_n = type1_hessian(x,4);
    for ( size_t i=4; i<order; i++ ) {
      result = ((i+alphaPoly-x)*d2Chdx2_n-2.*type1_gradient(x,order)-i*d2Chdx2_nminus1)/alphaPoly;
      if (i != order-1) {
	d2Chdx2_nminus1 = d2Chdx2_n;
	d2Chdx2_n       = result;
      }
    }
    break;
  }
  }
  return result;
};


Real CharlierOrthogPolynomial::norm_squared( unsigned short order )
{
  return std::pow( alphaPoly, order ) * boost::math::factorial<Real>( (Real)order );
}

} // namespace Pecos
