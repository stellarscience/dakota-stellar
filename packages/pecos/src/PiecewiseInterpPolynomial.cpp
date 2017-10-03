/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        PiecewiseInterpPolynomial
//- Description:  Implementation code for PiecewiseInterpPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "PiecewiseInterpPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.hpp"
#endif

namespace Pecos {


void PiecewiseInterpPolynomial::precompute_data()
{
  // For num_interp_pts == 1, we unconditionally return value = 1, grad = 0.
  size_t num_interp_pts = interpPts.size();
  if (num_interp_pts < 1) {
    PCerr << "Error: PiecewiseInterpPolynomial requires at least one point."
	  << std::endl;
    abort_handler(-1);
  }

  // for equally spaced pts (requires at least 2 points at interval bounds):
  if (num_interp_pts > 1 && collocRule == NEWTON_COTES) {
    size_t num_intervals = num_interp_pts - 1;
    interpInterval = (interpPts[num_intervals] - interpPts[0])/num_intervals;
  }

  // for non-equidistant points
  //if (num_interp_pts > 1 && collocRule != NEWTON_COTES) {
  //  size_t i, num_intervals = num_interp_pts - 1;
  //  interpIntervals.size(num_intervals);
  //  for (i=0; i<num_intervals; ++i)
  //    interpIntervals[i] = interpPts[k+1] - interpPts[k];
  //}
}


/** Compute value of the piecewise interpolation polynomial corresponding
    to interpolation point i. */
Real PiecewiseInterpPolynomial::type1_value(Real x, unsigned short i)
{
  // handle special case of a single interpolation point
  size_t num_interp_pts = interpPts.size();
  if (num_interp_pts == 1)
    return 1.;

  // does x lie within interval corresponding to interpolation point i
  Real pt_i = interpPts[i], t1_val;
  switch (basisPolyType) {
  case PIECEWISE_LINEAR_INTERP:
    switch (collocRule) {
    case NEWTON_COTES: {
      // linear spline interpolant with equidistant pts on [a,b]
      Real abs_dist = std::abs(x - pt_i);
      t1_val = (abs_dist < interpInterval) ? 1. - abs_dist/interpInterval : 0.;
      break;
    }
    default:
      // linear spline interpolant with general point spacing on [a,b];
      // forward/backward looking indices protected by x{<,>}pt_i (closed rules)
      if (x == pt_i) // for x=pt_i: val=1, not 0
	t1_val = 1.;
      else if (x < pt_i && x > interpPts[i-1])
	t1_val = 1. - (x - pt_i)/(interpPts[i-1] - pt_i);
      else if (x > pt_i && x < interpPts[i+1])
	t1_val = 1. - (x - pt_i)/(interpPts[i+1] - pt_i);
      else
	t1_val = 0.;
      break;
    }
    break;
  case PIECEWISE_QUADRATIC_INTERP:
    switch (collocRule) {
    case NEWTON_COTES: {
      // quadratic spline interpolant with equidistant pts on [a,b]
      Real dist = x - pt_i;
      if (std::abs(dist) < interpInterval) {
	Real ratio = dist/interpInterval;
	t1_val = 1.-ratio*ratio;
      }
      else
	t1_val = 0.;
      break;
    }
    default:
      // quadratic spline interpolant with general point spacing on [a,b]
      if (i == 0) { // 1-sided as equidistant 2-sided
	Real pt_ip1 = interpPts[1];
	if (x < pt_ip1) {
	  Real ratio = (x - pt_i)/(pt_ip1 - pt_i);
	  t1_val = 1.-ratio*ratio;
	}
	else t1_val = 0.;
      }
      else if (i == num_interp_pts-1) { // 1-sided as equidistant 2-sided
	Real pt_im1 = interpPts[i-1];
	if (x > pt_im1) {
	  Real ratio = (x - pt_i)/(pt_i - pt_im1);
	  t1_val = 1.-ratio*ratio;
	}
	else t1_val = 0.;
      }
      else { // 2-sided non-equidistant
	Real pt_im1 = interpPts[i-1], pt_ip1 = interpPts[i+1];
	// Note: this interpolant does not have zero derivative at x=pt_i
	t1_val = (x > pt_im1 && x < pt_ip1) ?
	  (x-pt_im1)*(pt_ip1-x)/(pt_i-pt_im1)/(pt_ip1-pt_i) : 0.;
      }
      break;
    }
    break;
  case PIECEWISE_CUBIC_INTERP:
    // cubic Hermite spline interpolant with general point spacing
    // glue together shape fn h01 from [i-1,i] with h00 from [i,i+1]
    if (x < pt_i) {
      Real pt_im1 = interpPts[i-1];
      if (x > pt_im1) { // left half interval: p_k+1=1, p_k=m_k=m_k+1=0
	Real t = (x-pt_im1)/(pt_i-pt_im1);
	t1_val = t*t*(3.-2.*t);     // h01(t)
      }
      else t1_val = 0.;
    }
    else if (x > pt_i) {
      Real pt_ip1 = interpPts[i+1];
      if (x < pt_ip1) { // right half interval: p_k=1, p_k+1=m_k=m_k+1=0
	Real t = (x - pt_i)/(pt_ip1-pt_i), tm1 = t-1.;
	t1_val = tm1*tm1*(1.+2.*t); // h00(t)
      }
      else t1_val = 0.;
    }
    else // x == pt_i
      t1_val = 1.;
    break;
  }
  return t1_val;
}


Real PiecewiseInterpPolynomial::type2_value(Real x, unsigned short i)
{
  // handle special case of a single interpolation point
  if (interpPts.size() == 1)
    switch (basisPolyType) {
    case PIECEWISE_LINEAR_INTERP: case PIECEWISE_QUADRATIC_INTERP:
      return 0.; break;
    case PIECEWISE_CUBIC_INTERP:
      return x; break; // integral of grad = 1 condition
    }

  Real t2_val;
  switch (basisPolyType) {
  case PIECEWISE_LINEAR_INTERP: case PIECEWISE_QUADRATIC_INTERP:
    t2_val = 0.;
    break;
  case PIECEWISE_CUBIC_INTERP: {
    // cubic Hermite spline interpolant with equidistant pts
    // glue together shape fn h11 from [i-1,i] with h10 from [i,i+1]
    Real pt_i = interpPts[i];
    if (x < pt_i) {
      Real pt_im1 = interpPts[i-1];
      if (x > pt_im1) { // left half interval: m_k+1=1, m_k=p_k=p_k+1=0
	Real interval = pt_i-pt_im1, t = (x-pt_im1)/interval;
	t2_val = interval*t*t*(t-1.); // interval*h11(t) -> h11(\xi)
      }
      else t2_val = 0.;
    }
    else if (x > pt_i) {
      Real pt_ip1 = interpPts[i+1];
      if (x < pt_ip1) { // right half interval: m_k=1, m_k+1=p_k=p_k+1=0
	Real interval = pt_ip1-pt_i, t = (x-pt_i)/interval, tm1 = t-1.;
	t2_val = interval*tm1*tm1*t;  // interval*h10(t) -> h10(\xi)
      }
      else t2_val = 0.;
    }
    else // x == pt_i
      t2_val = 0.;
    break;
  }
  }
  return t2_val;
}


/** Compute derivative with respect to x of the piecewise interpolation
    polynomial corresponding to interpolation point i. */
Real PiecewiseInterpPolynomial::type1_gradient(Real x, unsigned short i)
{ 
  // handle special case of a single interpolation point
  size_t num_interp_pts = interpPts.size();
  if (num_interp_pts == 1)
    return 0.; // derivative of value = 1 condition

  // does x lie within interval corresponding to interpolation point i
  Real pt_i = interpPts[i], t1_grad;
  switch (basisPolyType) {
  case PIECEWISE_LINEAR_INTERP:
    switch (collocRule) {
    case NEWTON_COTES: {
      // linear spline interpolant with equidistant pts on [a,b]
      Real dist = x - pt_i, abs_dist = std::abs(dist);
      // zero gradient at transition points and outside local support
      if (abs_dist == 0. || abs_dist >= interpInterval)
	t1_grad = 0.;
      else
	t1_grad = (dist < 0.) ? 1./interpInterval : -1./interpInterval;
      break;
    }
    default:
      if (x < pt_i && x > interpPts[i-1])
	t1_grad =  1./(pt_i - interpPts[i-1]);
      else if (x > pt_i && x < interpPts[i+1])
	t1_grad = -1./(interpPts[i+1] - pt_i);
      else // includes x == pt_i
	t1_grad = 0.;
      break;
    }
    break;
  case PIECEWISE_QUADRATIC_INTERP:
    switch (collocRule) {
    case NEWTON_COTES: {
      Real dist = x - pt_i;
      t1_grad = (std::abs(dist) < interpInterval) ?
	-2.*dist/(interpInterval*interpInterval) : 0.; // includes dist=0 case
      break;
    }
    default:
      // quadratic spline interpolant with general point spacing on [a,b]
      if (i == 0) { // 1-sided as equidistant 2-sided
	Real pt_ip1 = interpPts[1];
	if (x < pt_ip1) {
	  Real interval = pt_ip1 - pt_i; // 1-sided as equidistant 2-sided
	  t1_grad = -2.*(x - pt_i)/(interval*interval);
	}
	else t1_grad = 0.;
      }
      else if (i == num_interp_pts-1) { // 1-sided as equidistant 2-sided
	Real pt_im1 = interpPts[i-1];
	if (x > pt_im1) {
	  Real interval = pt_i - pt_im1; // 1-sided as equidistant 2-sided
	  t1_grad = -2.*(x - pt_i)/(interval*interval);
	}
	else t1_grad = 0.;
      }
      else { // 2-sided non-equidistant
	Real pt_im1 = interpPts[i-1], pt_ip1 = interpPts[i+1];
	// Note: this interpolant does not have zero derivative at x=pt_i
	if (x > pt_im1 && x < pt_ip1) {
	  Real interval1 = pt_i - pt_im1, interval2 = pt_ip1 - pt_i,
	    ratio1 = (x - pt_im1)/interval1, ratio2 = (pt_ip1 - x)/interval2;
	  // ratio1 * dratio2/dx + ratio2 * dratio1/dx
	  t1_grad = ratio2/interval1 - ratio1/interval2;
	}
	else t1_grad = 0.;
      }
      break;
    }
    break;
  case PIECEWISE_CUBIC_INTERP:
    // cubic Hermite spline interpolant with equidistant pts
    // glue together shape fn h01 from [i-1,i] with h00 from [i,i+1]
    if (x < pt_i) {
      Real pt_im1 = interpPts[i-1];
      if (x > pt_im1) { // left half interval: m_k+1=1, m_k=p_k=p_k+1=0
	Real interval = pt_i-pt_im1, t = (x-pt_im1)/interval, dt_dx=1./interval;
	t1_grad = 6.*t*(1.-t)*dt_dx; // dh01/dt * dt/dx
      }
      else t1_grad = 0.;
    }
    else if (x > pt_i) {
      Real pt_ip1 = interpPts[i+1];
      if (x < pt_ip1) {	// right half interval: m_k=1, m_k+1=p_k=p_k+1=0
	Real interval = pt_ip1-pt_i, t = (x-pt_i)/interval, dt_dx = 1./interval;
	t1_grad = 6.*t*(t-1.)*dt_dx; // dh00/dt * dt/dx
      }
      else t1_grad = 0.;
    }
    else // x == pt_i
      t1_grad = 0.;
    break;
  }
  return t1_grad;
}


Real PiecewiseInterpPolynomial::type2_gradient(Real x, unsigned short i)
{
  // handle special case of a single interpolation point
  if (interpPts.size() == 1)
    switch (basisPolyType) {
    case PIECEWISE_LINEAR_INTERP: case PIECEWISE_QUADRATIC_INTERP:
      return 0.; break;
    case PIECEWISE_CUBIC_INTERP:
      return 1.; break;
    }

  Real t2_grad;
  switch (basisPolyType) {
  case PIECEWISE_LINEAR_INTERP: case PIECEWISE_QUADRATIC_INTERP:
    t2_grad = 0.;
    break;
  case PIECEWISE_CUBIC_INTERP: {
    // cubic Hermite spline interpolant with equidistant pts
    // glue together shape fn h11 from [i-1,i] with h10 from [i,i+1]
    Real pt_i = interpPts[i];
    if (x < pt_i) {
      Real pt_im1 = interpPts[i-1];
      if (x > pt_im1) { // left half interval: m_k+1=1, m_k=p_k=p_k+1=0
	Real interval = pt_i-pt_im1, t = (x-pt_im1)/interval;//dt_dx=1./interval
	t2_grad = t*(3.*t-2.);//*interval*dt_dx (terms cancel)
      }
      else t2_grad = 0.;
    }
    else if (x > pt_i) {
      Real pt_ip1 = interpPts[i+1];
      if (x < pt_ip1) { // right half interval: m_k=1, m_k+1=p_k=p_k+1=0
	Real interval = pt_ip1-pt_i, t = (x-pt_i)/interval;//,dt_dx=1./interval;
	t2_grad = t*(3.*t-4.)+1.;//*interval*dt_dx (terms cancel)
      }
      else t2_grad = 0.;
    }
    else // x == pt_i
      t2_grad = 1.;
    break;
  }
  }
  return t2_grad;
}


const RealArray& PiecewiseInterpPolynomial::
collocation_points(unsigned short order)
{
  // pull this outside block below since order=0 is initial interpPts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum order (1) in PiecewiseInterp"
	  << "Polynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  // Computation of interpPts depends only on collocRule (not on basisPolyType).
  // Points are defined on [-1,1] (unlike Ma & Zabaras, Jakeman, etc.).
  // Bypass webbur::{hce,hcc}_compute_points() since it replicates points
  // in an array of size 2*order to match type1/2 weight aggregation.

  bool mode_err = false;
  if (interpPts.size() != order) { // if not already computed
    interpPts.resize(order);
    if (order == 1)
      interpPts[0] = 0.;
    else if (collocRule == NEWTON_COTES) {
      Real val = 2./((Real)(order - 1));
      for (unsigned short i=0; i<order; ++i)
	interpPts[i] = val*i - 1.;
      // update interpInterval (could also call precompute_data())
      size_t num_intervals = order - 1;
      interpInterval = (interpPts[num_intervals] - interpPts[0])/num_intervals;
    }
    else if (collocRule == CLENSHAW_CURTIS) {
#ifdef HAVE_SPARSE_GRID
      webbur::clenshaw_curtis_compute_points(order, &interpPts[0]);
#else
      PCerr << "Error: configuration with VPISparseGrid package required in "
	    << "PiecewiseInterpPolynomial::collocation_points()."<< std::endl;
      abort_handler(-1);
#endif
    }
    else
      mode_err = true;
  }

  if (mode_err) {
    PCerr << "Error: unsupported interpolation mode in "
	  << "PiecewiseInterpPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  return interpPts;
}


const RealArray& PiecewiseInterpPolynomial::
type1_collocation_weights(unsigned short order)
{
  // pull this outside block below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum order (1) in PiecewiseInterp"
	  << "Polynomial::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  bool mode_err = false;
  if (interpPts.size() != order)
    collocation_points(order);
  if (type1InterpWts.size() != order) { // if not already computed
    type1InterpWts.resize(order);
    if (order == 1)
      type1InterpWts[0] = 1.;
    else
      switch (basisPolyType) {
      case PIECEWISE_LINEAR_INTERP: case PIECEWISE_CUBIC_INTERP:
	//   Left end:  (x_1 - a)/2/(b-a)           = (x_1 - a)/4
	//   Interior:  (x_{i+1} - x_{i-1})/2/(b-a) = (x_{i+1} - x_{i-1})/4
	//   Right end: (b - x_{order-2})/2/(b-a)   = (b - x_{order-2})/4
	// Bypass webbur::{hce,hcc}_compute_weights() since it aggregates
	// type1/2 weights in a single array of size 2*order
	if (collocRule == NEWTON_COTES) {
	  Real val = interpInterval/4.;
	  type1InterpWts[0] = type1InterpWts[order-1] = val; // left/right ends
	  val *= 2.;
	  for (unsigned short i=1; i<order-1; ++i)
	    type1InterpWts[i] = val; // interior
	}
	else if (collocRule == CLENSHAW_CURTIS) {
	  // bases for left and right ends span one interp interval
	  type1InterpWts[0]       = (interpPts[1]      -interpPts[0])/4.;
	  type1InterpWts[order-1] = (interpPts[order-1]-interpPts[order-2])/4.;
	  // interior bases span two interp intervals
	  for (unsigned short i=1; i<order-1; ++i)
	    type1InterpWts[i] = (interpPts[i+1]-interpPts[i-1])/4.;
	}
	else
	  mode_err = true;
	break;
      case PIECEWISE_QUADRATIC_INTERP:
	mode_err = true; break;
      }
  }

  if (mode_err) {
    PCerr << "Error: unsupported interpolation mode in PiecewiseInterp"
	  << "Polynomial::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  return type1InterpWts;
}


const RealArray& PiecewiseInterpPolynomial::
type2_collocation_weights(unsigned short order)
{
  // pull this outside block below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum order (1) in PiecewiseInterp"
	  << "Polynomial::type2_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  bool mode_err = false;
  switch (basisPolyType) {
  case PIECEWISE_LINEAR_INTERP: case PIECEWISE_QUADRATIC_INTERP:
    if (!type2InterpWts.empty())
      type2InterpWts.clear();
    break;
  case PIECEWISE_CUBIC_INTERP:
    //   Left end:   (x_1 - a)^2/12/(b-a)         =  (x_1 - a)/24
    //   Interior:   (x_{i+1} - 2x_i + x_{i-1})(x_{i+1} - x_{i-1})/12/(b-a)
    //            =  (x_{i+1} - 2x_i + x_{i-1})(x_{i+1} - x_{i-1})/24
    //   Right end: -(b - x_{order-2})^2/12/(b-a) = -(b - x_{order-2})/24
    // Bypass webbur::{hce,hcc}_compute_weights() since it aggregates
    // type1/2 weights in a single array of size 2*order
    if (interpPts.size() != order)
      collocation_points(order);
    if (type2InterpWts.size() != order) { // if not already computed
      type2InterpWts.resize(order);
      if (order == 1)
	type2InterpWts[0] = 0.; // TO DO: verify
      else if (collocRule == NEWTON_COTES) {
	Real val = interpInterval*interpInterval/24.;
	type2InterpWts[0]       =  val; //  left end
	type2InterpWts[order-1] = -val; // right end
	for (unsigned short i=1; i<order-1; ++i)
	  type2InterpWts[i] = 0.;       //  interior
      }
      else if (collocRule == CLENSHAW_CURTIS) {
	Real val = interpPts[1] - interpPts[0];
	type2InterpWts[0]       =  val*val/24.;  //  left end
	val = interpPts[order-1] - interpPts[order-2];
	type2InterpWts[order-1] = -val*val/24.;  // right end
	for (unsigned short i=1; i<order-1; ++i) //  interior
	  type2InterpWts[i] = (interpPts[i+1] - interpPts[i-1])*
	    (interpPts[i+1] - 2.*interpPts[i] + interpPts[i-1])/24.;
      }
      else
	mode_err = true;
    }
    break;
  }

  if (mode_err) {
    PCerr << "Error: unsupported interpolation mode in PiecewiseInterp"
	  << "Polynomial::type2_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  return type2InterpWts;
}

} // namespace Pecos
