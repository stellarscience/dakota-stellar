/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LagrangeInterpPolynomial
//- Description:  Implementation code for LagrangeInterpPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "LagrangeInterpPolynomial.hpp"


namespace Pecos {

/** Pre-compute denominator products that are only a function of the
    interpolation points. */
void LagrangeInterpPolynomial::precompute_data()
{
  // precompute w_j for all forms of Lagrange interpolation
  size_t i, j, num_interp_pts = interpPts.size(); Real prod;
  if (bcWeights.empty())
    bcWeights.sizeUninitialized(num_interp_pts);
  for (i=0; i<num_interp_pts; ++i) {
    prod = 1.;
    Real interp_pt_i = interpPts[i];
    for (j=0; j<num_interp_pts; ++j)
      if (i != j)
	prod *= interp_pt_i - interpPts[j];
    bcWeights[i] = 1. / prod; // prefer repeated multiplication to division
  }
}


/** Shared initialization code. */
void LagrangeInterpPolynomial::
init_new_point(Real x, short request_order, short& compute_order)
{
  // define compute_order from requested and available data.  If grad factors
  // are requested, value factors will also be needed (multidimensional grad
  // components involve values in non-differentiated dimensions).
  if (x == newPoint) {
    compute_order = request_order - (newPtOrder & request_order);
    if (request_order == 2 && !(newPtOrder & 1))
      compute_order |= 1;                           // augment request
    if (compute_order) newPtOrder |= compute_order; // augment total available
    else               return;
  }
  else {
    newPtOrder = compute_order = (request_order & 2) ? 3 : request_order;
    newPoint   = x; exactIndex = exactDeltaIndex = _NPOS;
  }
}


/** Define the bcValueFactors (and exactIndex if needed) corresponding to x. */
void LagrangeInterpPolynomial::set_new_point(Real x, short request_order)
{
  short compute_order;
  init_new_point(x, request_order, compute_order);
  allocate_factors(compute_order);

  // second form of barycentric interpolation: precompute value factors and
  // grad factor terms or identify exactIndex
  size_t j, num_interp_pts = interpPts.size();
  Real diff_inv_sum; RealVector diffs;
  if (exactIndex == _NPOS) { // exact match may have been previously detected
    // first, compute diffs vector or exactIndex
    diffs.sizeUninitialized(num_interp_pts);
    for (j=0; j<num_interp_pts; ++j) {
      diffs[j] = newPoint - interpPts[j];
      if (diffs[j] == 0.) // no tolerance needed due to favorable stability
	{ exactIndex = exactDeltaIndex = j; break; }
    }

    // now compute value factors and grad factor terms based on exactIndex
    if (exactIndex == _NPOS) {
      if (compute_order & 1) bcValueFactorSum = 0.;
      if (compute_order & 2) { diffProduct = 1.; diff_inv_sum = 0.; }
      for (j=0; j<num_interp_pts; ++j) {
	if (compute_order & 1)
	  bcValueFactorSum += bcValueFactors[j] = bcWeights[j] / diffs[j];
	if (compute_order & 2)
	  { diffProduct *= diffs[j]; diff_inv_sum += 1. / diffs[j]; }
      }
    }
  }
  // catch previous or current identification of exactIndex
  if (exactIndex != _NPOS && (compute_order & 1)) {
    bcValueFactors = 0.; bcValueFactors[exactIndex] = 1.;
    //bcValueFactorSum = 1.; // not currently used if exactIndex
  }

  // second form of barycentric interpolation: precompute gradient factors
  if (compute_order & 2) {
    // bcGradFactors (like bcValueFactors) differ from the actual gradient
    // values by diffProduct, which gets applied after a tensor summation
    // using the barycentric gradient scaling.
    if (exactIndex == _NPOS)
      for (j=0; j<num_interp_pts; ++j) // bcValueFactors must be available
	bcGradFactors[j] = bcValueFactors[j] // * diffProduct
	  * (diff_inv_sum - 1. / diffs[j]); // back out jth inverse diff
    else { // Berrut and Trefethen, 2004
      // for this case, bcGradFactors are the actual gradient values
      // and no diffProduct scaling needs to be subsequently applied
      bcGradFactors[exactIndex] = 0.;
      for (j=0; j<num_interp_pts; ++j)
	if (j != exactIndex)
	  bcGradFactors[exactIndex] -= bcGradFactors[j] = bcWeights[j]
	    / bcWeights[exactIndex] / (interpPts[exactIndex] - interpPts[j]);
    }
  }
}


/** Define the bcValueFactors (and exactIndex if needed) corresponding to x. */
void LagrangeInterpPolynomial::
set_new_point(Real x, short request_order, const UShortArray& delta_key)
{
  short compute_order;
  // TO DO: capture change in delta_key...
  // Note: no longer any reuse --> just move code cases here and encapsulate
  //       any lower level shared logic
  init_new_point(x, request_order, compute_order);//, delta_key);
  allocate_factors(compute_order);

  // second form of barycentric interpolation: precompute value factors and
  // grad factor terms or identify exactIndex
  size_t j, vj, num_interp_pts = interpPts.size(),
    num_delta_pts = delta_key.size();
  Real diff_inv_sum; RealVector diffs;

  if (exactIndex == _NPOS) { // exact match may have been previously detected
    // compute diffs vector and detect exactIndex within all of interpPts
    diffs.sizeUninitialized(num_interp_pts);
    for (j=0; j<num_interp_pts; ++j) {
      diffs[j] = newPoint - interpPts[j];
      if (diffs[j] == 0.) // no tol reqd due to favorable stability
	{ exactIndex = j; break; }
    }
    // detect exactDeltaIndex within delta points only
    // (see, for example, InterpPolyApproximation::barycentric_exact_index())
    exactDeltaIndex = (exactIndex == _NPOS) ? _NPOS :
      find_index(delta_key, exactIndex);

    // now compute value factors and grad factor terms based on exactIndex
    if (exactIndex == _NPOS) {
      if (compute_order & 1) bcValueFactorSum = 0.;
      if (compute_order & 2) { diffProduct = 1.; diff_inv_sum = 0.; }
      for (j=0; j<num_interp_pts; ++j) {
	if (compute_order & 1)
	  bcValueFactorSum += bcValueFactors[j] = bcWeights[j] / diffs[j];
	if (compute_order & 2)
	  { diffProduct *= diffs[j]; diff_inv_sum += 1. / diffs[j]; }
      }
    }
  }
  // catch previous or current identification of exactIndex
  if (exactIndex != _NPOS && (compute_order & 1)) {
    bcValueFactors = 0.; bcValueFactors[exactIndex] = 1.;
    //bcValueFactorSum = 1.; // not currently used if exactIndex
  }

  // second form of barycentric interpolation: precompute gradient factors
  if (compute_order & 2) {
    // bcGradFactors (like bcValueFactors) differ from the actual gradient
    // values by diffProduct, which gets applied after a tensor summation
    // using the barycentric gradient scaling.
    if (exactIndex == _NPOS)
      for (j=0; j<num_delta_pts; ++j) { // bcValueFactors must be available
	vj = delta_key[j];
	bcGradFactors[vj] = bcValueFactors[vj] // * diffProduct
	  * (diff_inv_sum - 1. / diffs[vj]); // back out vj-th inverse diff
      }
    else { // Berrut and Trefethen, 2004
      // for this case, bcGradFactors are the actual gradient values
      // and no diffProduct scaling needs to be subsequently applied
      bcGradFactors[exactIndex] = 0.;
      for (j=0; j<num_interp_pts; ++j)
	if (j != exactIndex)
	  bcGradFactors[exactIndex] -= bcGradFactors[j] = bcWeights[j]
	    / bcWeights[exactIndex] / (interpPts[exactIndex] - interpPts[j]);
    }
  }
}


/** Compute value of the Lagrange polynomial (1st barycentric form)
    corresponding to interpolation point i using data from previous
    call to set_new_point(). */
Real LagrangeInterpPolynomial::type1_value(unsigned short i)
{
  // second form of the barycentric interpolation formula
  if (exactIndex == _NPOS) return bcValueFactors[i] * diffProduct;
  else                     return (exactIndex == i) ? 1. : 0.;

  /*
  // first form of the barycentric interpolation formula
  if (exactIndex == _NPOS)
    return bcWeights[i] * diffProduct / (newPoint - interpPts[i]);
  else
    return (exactIndex == i) ? 1. : 0.;
  */
}


/** Compute derivative with respect to x of the Lagrange polynomial
    (1st barycentric form) corresponding to interpolation point i
    using data from previous call to set_new_point(). */
Real LagrangeInterpPolynomial::type1_gradient(unsigned short i)
{
  // second form of the barycentric interpolation formula
  return (exactIndex == _NPOS) ?
    bcGradFactors[i] * diffProduct : bcGradFactors[i];

  /*
  // first form of the barycentric interpolation formula
  if (exactIndex == _NPOS) {
    Real sum = 0.,
      t1_i = bcWeights[i] * diffProduct / (newPoint - interpPts[i]);
    size_t j, num_interp_pts = interpPts.size();
    for (j=0; j<num_interp_pts; ++j)
      if (j != i)
	sum += t1_i / (newPoint - interpPts[j]);
    return sum;
  }
  else // double loop fallback
    return type1_gradient(newPoint, i);
  */
}


/** Compute value of the Lagrange polynomial (traditional characteristic
    polynomial form) corresponding to interpolation point i. */
Real LagrangeInterpPolynomial::type1_value(Real x, unsigned short i)
{
  size_t j, num_interp_pts = interpPts.size();
  Real t1_val = bcWeights[i];
  for (j=0; j<num_interp_pts; ++j)
    if (i != j)
      t1_val *= x - interpPts[j];
  return t1_val;
}


/** Compute derivative with respect to x of the Lagrange polynomial (traditional
    characteristic polynomial form) corresponding to interpolation point i. */
Real LagrangeInterpPolynomial::type1_gradient(Real x, unsigned short i)
{ 
  size_t j, k, num_interp_pts = interpPts.size();
  Real sum = 0., prod;
  for (j=0; j<num_interp_pts; ++j) {
    if (j != i) {
      prod = 1.;
      for (k=0; k<num_interp_pts; ++k)
	if (k != j && k != i)
	  prod *= x - interpPts[k];
      sum += prod;
    }
  }
  return sum * bcWeights[i];
}

} // namespace Pecos
