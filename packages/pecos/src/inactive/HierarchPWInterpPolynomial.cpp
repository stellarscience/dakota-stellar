/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchPWInterpPolynomial
//- Description:  Implementation code for HierarchPWInterpPolynomial class
//-               
//- Owner:        Christopher Miller, University of Maryland
//- Contact:      cmiller@math.umd.edu

#include "HierarchPWInterpPolynomial.hpp"

namespace Pecos {

/// Point evaluation of a basis element.
Real HierarchPWInterpPolynomial::
type1_value(const Real& x, const unsigned int i) 
{
  const Real    nodalPoint = pointSet.get_interp_point(i);
  const Real  leftEndPoint = pointSet.get_left_neighbor(i);
  const Real rightEndPoint = pointSet.get_right_neighbor(i);

  Real t1_val;
  switch (basisPolyType) {
  case PIECEWISE_QUADRATIC_INTERP:
    PCerr << "Warning: Quadratic interpolation not currently implemented in "
	  << "HierarchPWInterpPolynomial. Defaulting to linear interpolation.";
    // no break, fall through
  case PIECEWISE_LINEAR_INTERP:
    if ( i == 0 ) {
      if ( ( x < leftEndPoint ) || ( x > rightEndPoint ) ) t1_val = 0;
      else t1_val = 1;
    } else {
      if ( ( x < nodalPoint ) && ( x >= leftEndPoint ) ) {
	t1_val = (x - leftEndPoint)/(nodalPoint - leftEndPoint);
      } else if ( ( x >= nodalPoint ) && ( x <= rightEndPoint ) ) {
	if ( rightEndPoint == nodalPoint ) t1_val = 1;
	else t1_val = (rightEndPoint - x)/(rightEndPoint - nodalPoint);
      } else t1_val = 0;
    }
    break;
  case PIECEWISE_CUBIC_INTERP:
    if ( i == 0 ) {
      if ( ( x < leftEndPoint ) || ( x > rightEndPoint ) ) t1_val = 0;
      else t1_val = 1;
    } else {
      if ( ( x < nodalPoint ) && ( x >= leftEndPoint ) ) {
	Real t = ( x - leftEndPoint )/( nodalPoint - leftEndPoint );
	t1_val = t*t*(3.-2.*t);
      }
      else if ( ( x >= nodalPoint ) && ( x < rightEndPoint ) ) {
	Real t = ( x - nodalPoint )/( rightEndPoint - nodalPoint );
	Real tm1 = t - 1.;
	t1_val = tm1*tm1*(1.+2.*t);
      } else if ( ( x == nodalPoint ) && ( x == rightEndPoint ) ) {
	t1_val = 1;
      }
      else t1_val = 0.;
    }
    break;
  }
  return t1_val;
}

/// Point evaluation of the gradient of a basis element.
Real HierarchPWInterpPolynomial::
type1_gradient(const Real& x, const unsigned int i) 
{
  const Real nodalPoint = pointSet.get_interp_point(i);
  const Real leftEndPoint = pointSet.get_left_neighbor(i);
  const Real rightEndPoint = pointSet.get_right_neighbor(i);
  Real t1_grad;
  switch (basisPolyType) {
  case PIECEWISE_QUADRATIC_INTERP:
    PCerr << "Warning: Quadratic interpolation not currently implemented in "
	  << "HierarchPWInterpPolynomial. Defaulting to linear interpolation.";
    // no break, fall through
  case PIECEWISE_LINEAR_INTERP:
    if ( i == 0 ) {
      t1_grad = 0;
    } else {
      if ( ( x < nodalPoint ) && ( x > leftEndPoint ) ) {
	t1_grad = 1/(nodalPoint - leftEndPoint);
      } else if ( ( x > nodalPoint ) && ( x < rightEndPoint ) ) {
	t1_grad = -1/(rightEndPoint - nodalPoint);
      } else t1_grad = 0;
    }
    break;
  case PIECEWISE_CUBIC_INTERP:
    if ( i == 0 ) {
      t1_grad = 0;
    } else { 
      if ( ( x < nodalPoint ) && ( x > leftEndPoint ) ) {
	Real interval = nodalPoint - leftEndPoint;
	Real t = ( x - leftEndPoint )/interval;
	Real dt_dx = 1./interval;
	t1_grad = 6.*t*(1.-t)*dt_dx; // dh01/dt * dt/dx
      }
      else if ( ( x > nodalPoint ) && ( x < rightEndPoint ) ) {
	Real interval = rightEndPoint - nodalPoint;
	Real t = ( x - nodalPoint )/interval;
	Real dt_dx = 1./interval;
	t1_grad = 6.*t*(t-1.)*dt_dx; // dh00/dt * dt/dx
      }
      else t1_grad = 0;
    }
    break;
  }
  return t1_grad;
}

Real HierarchPWInterpPolynomial::
type2_value(const Real& x, const unsigned int i)
{
  const Real nodalPoint = pointSet.get_interp_point(i);
  const Real leftEndPoint = pointSet.get_left_neighbor(i);
  const Real rightEndPoint = pointSet.get_right_neighbor(i);

  Real t2_val;
  switch (basisPolyType) {
  case PIECEWISE_CUBIC_INTERP:
    if ( i == 0 ) {
      t2_val = 0;
    }
    else if (x <= nodalPoint && x > leftEndPoint) {
      // left half interval: m_k+1=1, m_k=p_k=p_k+1=0
      Real interval = nodalPoint - leftEndPoint;
      Real t = (x-leftEndPoint)/interval;
      t2_val = interval*t*t*(t-1.); // interval*h11(t) -> h11(\xi)
    }
    // final condition here prevents double evaluations in the general case
    else if (x >= nodalPoint && x < rightEndPoint ) { 
      // right half interval: m_k=1, m_k+1=p_k=p_k+1=0
      Real interval = rightEndPoint-nodalPoint;
      Real t = (x-nodalPoint)/interval;
      Real tm1 = t-1.;
      t2_val = interval*tm1*tm1*t;  // interval*h10(t) -> h10(\xi)
    }
    else
      t2_val = 0.;
    break;
  default:
    t2_val = 0;
    break;
  }
  return t2_val;
}

Real HierarchPWInterpPolynomial::
type2_gradient(const Real& x, const unsigned int i)
{
  const Real nodalPoint = pointSet.get_interp_point(i);
  const Real leftEndPoint = pointSet.get_left_neighbor(i);
  const Real rightEndPoint = pointSet.get_right_neighbor(i);

  Real t2_grad;
  switch (basisPolyType) {
  case PIECEWISE_CUBIC_INTERP:
    if ( i == 0 ) {
      t2_grad = 0;
    } else {
      if (x < nodalPoint && x > leftEndPoint) {
	// left half interval: m_k+1=1, m_k=p_k=p_k+1=0
	Real interval = nodalPoint - leftEndPoint;
	Real t = (x-leftEndPoint)/interval;
	t2_grad = t*(3.*t-2.); 
      }
      else if (x > nodalPoint && x < rightEndPoint) {
	// right half interval: m_k=1, m_k+1=p_k=p_k+1=0
	Real interval = rightEndPoint-nodalPoint;
	Real t = (x-nodalPoint)/interval;
	t2_grad = t*(3.*t-4.)+1;
      }
      else if ( x == nodalPoint ) {
	t2_grad = 1;
      }
      else
	t2_grad = 0.;
    }
    break;
  default:
    t2_grad = 0;
    break;
  }
  return t2_grad;
}

} // End namespace Pecos
