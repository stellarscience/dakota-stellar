/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedNodalInterpPolyApproxData
//- Description:  Implementation code for SharedNodalInterpPolyApproxData class
//-               
//- Owner:        Mike Eldred

#include "SharedNodalInterpPolyApproxData.hpp"
#include "TensorProductDriver.hpp"
#include "CombinedSparseGridDriver.hpp"
#include "InterpolationPolynomial.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG
//#define VBD_DEBUG

namespace Pecos {


void SharedNodalInterpPolyApproxData::allocate_data()
{
  SharedInterpPolyApproxData::allocate_data();

  // We are migrating towards consistent usage of INTERPOLATION_OF_PRODUCTS,
  // but its usage of higher-order reinterpolation of covariance is currently
  // too slow for production usage.  Thus, we only activate it when needed to
  // support new capability, such as gradient-enhanced interpolation.
  momentInterpType = (basisConfigOptions.useDerivs) ?
    REINTERPOLATION_OF_PRODUCTS : PRODUCT_OF_INTERPOLANTS_FAST;

  // Map out current flow to integration driver:
  // > NonDStochCollocation::initialize_u_space_model() constructs driver_basis
  //   which is passed along to Pecos::*Driver::initialize_grid(poly_basis) at
  //   construct time
  // > this function is invoked at run time as part of ApproximationInterface::
  //   build_approximation()

  // need full integrand precision 2m-1 for interpolant order m-1
  // TO DO: discrepancies in growth will necessitate management of alternate
  // quad orders.  For example, an interpolatory nested rule has grown to 9
  // points at level 3 --> order 8 interpolant, but the non-nested level 3
  // Gauss rule only has 7 pts with prec = 13 --> insufficient for variance.

  // Determine if alt grid is necessary (interpolatory rules, nested Gauss,
  // etc.) and appropriate (smooth global basis; not piecewise).
  // basisConfigOptions.piecewiseBasis is strictly enforced, but 
  // basisConfigOptions.{useDerivs,nestedRules} are requests that may not
  // be acted upon depending on driverMode and u_types (resort to checking
  // collocation rules for these cases).
  bool alt_moment_grid = false;
  const std::vector<BasisPolynomial>& driver_basis
    = driverRep->polynomial_basis();
  /* TO DO: re-activate this block
  if (!basisConfigOptions.piecewiseBasis) {
    // check for interpolatory / nested rules or grad-enhanced interpolants
    for (size_t i=0; i<numVars; ++i) {
      short rule = driver_basis[i].collocation_rule();
      if (rule == CLENSHAW_CURTIS || rule == FEJER2       ||
	  rule == GAUSS_PATTERSON || rule == GENZ_KEISTER ||
	  driver_basis[i].basis_type() == HERMITE_INTERP)
	{ alt_moment_grid = true; break; }
    }
  }
  */

  // initialize expMomentIntDriver if needed
  if (alt_moment_grid) { // modify driver_basis to employ Gaussian integration
    // This approach reduces data reqmts (u_types, adp) and increases reuse,
    // relative to building from scratch using SharedOPAData::construct_basis().
    std::vector<BasisPolynomial> alt_driver_basis(numVars);
    for (size_t i=0; i<numVars; ++i) {
      const BasisPolynomial& basis_i =     driver_basis[i];
      BasisPolynomial&   alt_basis_i = alt_driver_basis[i];
      switch (basis_i.collocation_rule()) {
      case GENZ_KEISTER:
	alt_basis_i = BasisPolynomial(HERMITE_ORTHOG,  GAUSS_HERMITE);  break;
      case CLENSHAW_CURTIS: case FEJER2: case GAUSS_PATTERSON:
	alt_basis_i = BasisPolynomial(LEGENDRE_ORTHOG, GAUSS_LEGENDRE);	break;
      default: // no change: use shallow copy (includes all parametric cases!)
	alt_basis_i = (basis_i.basis_type() == HERMITE_INTERP) ?
	  BasisPolynomial(LEGENDRE_ORTHOG, GAUSS_LEGENDRE) : basis_i;
	break;
      }
      // TO DO: any quad order modifications need to be handled downstream
    }
    /* This approach requires u_types and adp, which aren't readily available
    BasisConfigOptions bc_options(false, false, false, false);
    SharedOrthogPolyApproxData::
      construct_basis(u_types, adp, bc_options, alt_driver_basis);
    */

    // suppress hierarchical sparse grid driver for this purpose
    short driver_type = (expConfigOptions.expCoeffsSolnApproach == QUADRATURE)
                      ? QUADRATURE : COMBINED_SPARSE_GRID;
    expMomentIntDriver = IntegrationDriver(driver_type);
    expMomentIntDriver.mode(INTEGRATION_MODE);
    if (driver_type == COMBINED_SPARSE_GRID) {
      CombinedSparseGridDriver* driver = (CombinedSparseGridDriver*)driverRep;
      CombinedSparseGridDriver* alt_driver
	= (CombinedSparseGridDriver*)expMomentIntDriver.driver_rep();
      alt_driver->growth_rate(driver->growth_rate());
      //alt_driver->refinement_control(driver->refinement_control());
      alt_driver->track_collocation_details(false);
      alt_driver->track_unique_product_weights(true); // non-default
    }
    expMomentIntDriver.initialize_grid(alt_driver_basis);
    // Note: rule/growth pointers initialized in CombinedSparseGridDriver::
    // initialize_grid()
  }
}


void SharedNodalInterpPolyApproxData::allocate_component_sobol()
{
  if (expConfigOptions.vbdFlag) {
    if (expConfigOptions.vbdOrderLimit == 1) // main effects only
      allocate_main_sobol();
    else { // main + interaction effects

      // One approach is to leverage PCE equivalence.  The exact order of the
      // interpolation polynomial (e.g., for nested rules, local or
      // gradient-enhanced interpolants) is not critical for defining
      // interactions; the issue is more the presence of constant dimensions.
      // While the collocation key has a very different meaning from the PCE
      // multi-index, the presence of non-zero's still indicates multi-point
      // interpolation and the presence of dimension effects.
      sobolIndexMap.clear();
      switch (expConfigOptions.expCoeffsSolnApproach) {
      case QUADRATURE: {
	TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
	multi_index_to_sobol_index_map(tpq_driver->collocation_key());
	break;
      }
      case COMBINED_SPARSE_GRID: {
	CombinedSparseGridDriver* csg_driver
	  = (CombinedSparseGridDriver*)driverRep;
	const IntArray&      sm_coeffs  = csg_driver->smolyak_coefficients();
	const UShort3DArray& colloc_key = csg_driver->collocation_key();
	size_t i, num_smolyak_indices = sm_coeffs.size();
	for (i=0; i<num_smolyak_indices; ++i)
	  if (sm_coeffs[i])
	    multi_index_to_sobol_index_map(colloc_key[i]);
	break;
      }
      }
      assign_sobol_index_map_values();

      // another approach: interrogate polynomialBasis[].interpolation_size()
      // or the quadrature/sparse level indices, again focusing on the presence
      // of constant dimensions (size = 1, level = 0).  But given the need to
      // regenerate the effect combinations from this reduced order data, the
      // collocation key idea seems preferable since it's already available.
      //polynomial_basis_to_sobol_indices();
    }
  }
}


void SharedNodalInterpPolyApproxData::increment_component_sobol()
{
  if (expConfigOptions.vbdFlag && expConfigOptions.vbdOrderLimit != 1) {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    if (csg_driver->smolyak_coefficients().back()) {
      reset_sobol_index_map_values();
      multi_index_to_sobol_index_map(csg_driver->collocation_key().back());
      assign_sobol_index_map_values();
    }
  }
}


void SharedNodalInterpPolyApproxData::
accumulate_barycentric(RealVector& t1_accumulator, const UShortArray& lev_index,
		       const UShortArray& key_p)
{
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  unsigned short li_v; size_t v, ei_v; unsigned short key_pv;
  // accumulate sums over variables with max key value
  for (v=1; v<numVars; ++v) {
    li_v = lev_index[v]; key_pv = key_p[v];
    BasisPolynomial& poly_v = polynomialBasis[li_v][v];
    if (li_v == 0) // single integration/interpolation weight = 1
      t1_accumulator[v]  = t1_accumulator[v-1];
    else if (randomVarsKey[v]) // integration
      t1_accumulator[v] += t1_accumulator[v-1] * t1_wts_1d[li_v][v][key_pv];
    else { 
      ei_v = poly_v.exact_index();
      if (ei_v == _NPOS)        // interpolation without exact match
	t1_accumulator[v] += t1_accumulator[v-1] *
	  poly_v.barycentric_value_factor(key_pv);
      else if (ei_v == key_pv)  // interpolation with exact match
	t1_accumulator[v]  = t1_accumulator[v-1];
    }
    t1_accumulator[v-1] = 0.;
    if (key_pv + 1 != poly_v.interpolation_size())
      break;
  }
}


void SharedNodalInterpPolyApproxData::
accumulate_horners(RealVector& t1_accumulator, const UShortArray& lev_index,
		   const UShortArray& key_p, const RealVector& x)
{
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  unsigned short li_v, key_pv;
  // accumulate sums over variables with max key value
  for (size_t v=1; v<numVars; ++v) {
    li_v = lev_index[v]; key_pv = key_p[v];
    BasisPolynomial& poly_v = polynomialBasis[li_v][v];
    if (li_v == 0)   // single integration/interpolation weight = 1
      t1_accumulator[v]  = t1_accumulator[v-1];
    else
      t1_accumulator[v] += (randomVarsKey[v]) ?
	t1_accumulator[v-1] * t1_wts_1d[li_v][v][key_pv] : // integration
	t1_accumulator[v-1] * poly_v.type1_value(x[v], key_pv); // interp
    t1_accumulator[v-1] = 0.;
    if (key_pv + 1 != poly_v.interpolation_size())
      break;
  }
}


void SharedNodalInterpPolyApproxData::
accumulate_horners(RealVector& t1_accumulator, RealMatrix& t2_accumulator,
		   const UShortArray& lev_index, const UShortArray& key_p,
		   const RealVector& x)
{
  Real *t2_accum_v, *t2_accum_vm1, t1_val, t1_wt;
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  const Real3DArray& t2_wts_1d = driverRep->type2_collocation_weights_1d();
  unsigned short li_v; size_t v, vm1, k; unsigned short key_pv;
  // accumulate sums over variables with max key value
  for (v=1; v<numVars; ++v) {
    li_v = lev_index[v]; key_pv = key_p[v]; vm1 = v-1;
    t2_accum_v = t2_accumulator[v]; t2_accum_vm1 = t2_accumulator[vm1];
    BasisPolynomial& poly_v = polynomialBasis[li_v][v];
    if (randomVarsKey[v]) { // integration
      if (li_v == 0) {
	t1_accumulator[v]  = t1_accumulator[vm1];       // t1 weight is 1
	//t2_accum_v[v]    = 0.;                        // t2 weight is 0
	for (k=0; k<numVars; ++k)
	  if (k != v)
	    t2_accum_v[k]  = t2_accum_vm1[k];           // t1 weight is 1
      }
      else {
	t1_wt = t1_wts_1d[li_v][v][key_pv];
	t1_accumulator[v] += t1_accumulator[vm1] * t1_wt;
	t2_accum_v[v]     += t2_accum_vm1[v] * t2_wts_1d[li_v][v][key_pv];
	for (k=0; k<numVars; ++k)
	  if (k != v)
	    t2_accum_v[k] += t2_accum_vm1[k] * t1_wt;
      }
    }
    else {                  // interpolation
      if (li_v == 0) {
	t1_accumulator[v] = t1_accumulator[vm1];         // t1 value is 1
	t2_accum_v[v] = t2_accum_vm1[v] * poly_v.type2_value(x[v],key_pv);
	for (k=0; k<numVars; ++k)
	  if (k != v)
	    t2_accum_v[k] = t2_accum_vm1[k];             // t1 value is 1
      }
      else {
	t1_val = poly_v.type1_value(x[v], key_pv);
	t1_accumulator[v] += t1_accumulator[vm1] * t1_val;
	t2_accum_v[v]     += t2_accum_vm1[v] *
	  poly_v.type2_value(x[v], key_pv);
	for (k=0; k<numVars; ++k)
	  if (k != v)
	    t2_accum_v[k] += t2_accum_vm1[k] * t1_val;
      }
    }
    t1_accumulator[vm1] = 0.;
    for (k=0; k<numVars; ++k)
      t2_accum_vm1[k] = 0.;
    if (key_pv + 1 != poly_v.interpolation_size())
      break;
  }
}


void SharedNodalInterpPolyApproxData::
accumulate_barycentric_gradient(RealMatrix& t1_accumulator,
				const UShortArray& lev_index,
				const UShortArray& key_p, const SizetArray& dvv)
{
  // accumulate sums over variables with max key value
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  Real *accum_v, *accum_vm1, t1_wt_v, bc_vf_v, bc_gf_v;
  size_t v, d, ei_v, num_deriv_vars = dvv.size();
  unsigned short li_v, key_pv;
  for (v=1; v<numVars; ++v) {
    li_v    = lev_index[v];      key_pv    = key_p[v];
    accum_v = t1_accumulator[v]; accum_vm1 = t1_accumulator[v-1];
    BasisPolynomial& poly_v = polynomialBasis[li_v][v];
    if (randomVarsKey[v]) { // cases 1,2
      if (li_v) {
	t1_wt_v = t1_wts_1d[li_v][v][key_pv];
	for (d=0; d<num_deriv_vars; ++d)
	  { accum_v[d] += accum_vm1[d] * t1_wt_v; accum_vm1[d] = 0.; }
      }
      else // t1 weight is 1
	for (d=0; d<num_deriv_vars; ++d)
	  { accum_v[d]  = accum_vm1[d];           accum_vm1[d] = 0.; }
    }
    else if (li_v) {
      ei_v = poly_v.exact_index();
      bc_gf_v = poly_v.barycentric_gradient_factor(key_pv);
      if (ei_v == _NPOS) { // bc_vf_v has general value
	bc_vf_v = poly_v.barycentric_value_factor(key_pv);
	for (d=0; d<num_deriv_vars; ++d) {
	  accum_v[d] += (v == dvv[d] - 1) ? accum_vm1[d] * bc_gf_v
	    : accum_vm1[d] * bc_vf_v;
	  accum_vm1[d] = 0.;
	}
      }
      else if (ei_v == key_pv) // bc_vf_v is 1
	for (d=0; d<num_deriv_vars; ++d) {
	  accum_v[d] += (v == dvv[d] - 1) ? accum_vm1[d] * bc_gf_v
	    : accum_vm1[d];
	  accum_vm1[d] = 0.;
	}
      else // bc_vf_v is 0
	for (d=0; d<num_deriv_vars; ++d) {
	  if (v == dvv[d] - 1) accum_v[d] += accum_vm1[d] * bc_gf_v;
	  accum_vm1[d] = 0.;
	}
    }
    else { // grad factor is zero, value factor is omitted
      for (d=0; d<num_deriv_vars; ++d) {
	if (v != dvv[d] - 1) accum_v[d] += accum_vm1[d];
	accum_vm1[d] = 0.;
      }
    }
    if (key_pv + 1 != poly_v.interpolation_size())
      break;
  }
}


void SharedNodalInterpPolyApproxData::
accumulate_horners_gradient(RealMatrix& t1_accumulator,
			    const UShortArray& lev_index,
			    const UShortArray& key_p, const SizetArray& dvv,
			    const RealVector& x)
{
  // accumulate sums over variables with max key value
  unsigned short key_pv, li_v;
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  Real *accum_v, *accum_vm1, t1_wt_v, t1_val, t1_grad, xv;
  size_t v, v2, d, num_deriv_vars = dvv.size();
  for (v=1; v<numVars; ++v) {
    li_v    = lev_index[v];      key_pv    = key_p[v];
    accum_v = t1_accumulator[v]; accum_vm1 = t1_accumulator[v-1];
    BasisPolynomial& poly_v = polynomialBasis[li_v][v];
    if (randomVarsKey[v]) { // cases 1,2
      if (li_v) {
	t1_wt_v = t1_wts_1d[li_v][v][key_pv];
	for (d=0; d<num_deriv_vars; ++d)
	  { accum_v[d] += accum_vm1[d] * t1_wt_v; accum_vm1[d] = 0.; }
      }
      else // t1 weight is 1
	for (d=0; d<num_deriv_vars; ++d)
	  { accum_v[d]  = accum_vm1[d];           accum_vm1[d] = 0.; }
    }
    else if (li_v) {
      t1_val = poly_v.type1_value(x[v], key_pv);
      for (d=0; d<num_deriv_vars; ++d) {
	accum_v[d] += (v == dvv[d] - 1) ?
	  accum_vm1[d] * poly_v.type1_gradient(x[v], key_pv) :
	  accum_vm1[d] * t1_val;
	accum_vm1[d] = 0.;
      }
    }
    else { // t1_grad is zero, t1_val is 1
      for (d=0; d<num_deriv_vars; ++d) {
	if (v != dvv[d] - 1) accum_v[d] += accum_vm1[d];
	accum_vm1[d] = 0.;
      }
    }
    if (key_pv + 1 != poly_v.interpolation_size())
      break;
  }
}


void SharedNodalInterpPolyApproxData::
accumulate_horners_gradient(RealMatrix& t1_accumulator,
			    RealMatrixArray& t2_accumulators,
			    const UShortArray& lev_index,
			    const UShortArray& key_p, const SizetArray& dvv,
			    const RealVector& x)
{
  // accumulate sums over variables with max key value
  unsigned short key_pv, li_v;
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  const Real3DArray& t2_wts_1d = driverRep->type2_collocation_weights_1d();
  Real *t1_accum_v, *t1_accum_vm1, *t2_accum_v, *t2_accum_vm1,
    t1_wt_v, t2_wt_v, t1_val, t2_val, t1_grad, xv;
  size_t v, v2, vm1, d, num_deriv_vars = dvv.size();
  for (v=1; v<numVars; ++v) {
    li_v = lev_index[v]; key_pv = key_p[v]; xv = x[v]; vm1 = v - 1;
    t1_accum_v = t1_accumulator[v]; t1_accum_vm1 = t1_accumulator[vm1];
    BasisPolynomial& poly_v = polynomialBasis[li_v][v];
    if (randomVarsKey[v]) { // case 2
      if (li_v) {
	t1_wt_v = t1_wts_1d[li_v][v][key_pv];
	t2_wt_v = t2_wts_1d[li_v][v][key_pv];
	for (d=0; d<num_deriv_vars; ++d) {
	  t1_accum_v[d]  += t1_accum_vm1[d] * t1_wt_v;
	  t1_accum_vm1[d] = 0.;
	  t2_accum_v      = t2_accumulators[d][v];
	  t2_accum_vm1    = t2_accumulators[d][vm1];
	  t2_accum_v[v]  += t2_accum_vm1[v] * t2_wt_v;
	  for (v2=0; v2<numVars; ++v2) {
	    if (v2 != v)
	      t2_accum_v[v2] += t2_accum_vm1[v2] * t1_wt_v;
	    t2_accum_vm1[v2]  = 0.;
	  }
	}
      }
      else { // case 2 with t1 weight = 1, t2 weight = 0
	for (d=0; d<num_deriv_vars; ++d) {
	  t1_accum_v[d] = t1_accum_vm1[d]; t1_accum_vm1[d] = 0.;
	  t2_accum_v    = t2_accumulators[d][v];
	  t2_accum_vm1  = t2_accumulators[d][vm1];
	  t2_accum_v[v] = 0.;
	  for (v2=0; v2<numVars; ++v2) {
	    if (v2 != v)
	      t2_accum_v[v2] = t2_accum_vm1[v2];
	    t2_accum_vm1[v2] = 0.;
	  }
	}
      }
    }
    else if (li_v) { // case 4
      t1_val = poly_v.type1_value(xv, key_pv);
      t2_val = poly_v.type2_value(xv, key_pv);
      for (d=0; d<num_deriv_vars; ++d) {
	t2_accum_v   = t2_accumulators[d][v];
	t2_accum_vm1 = t2_accumulators[d][vm1];
	if (dvv[d] - 1 == v) {
	  t1_grad = poly_v.type1_gradient(xv, key_pv);
	  t1_accum_v[d]  += t1_accum_vm1[d] * t1_grad;
	  t1_accum_vm1[d] = 0.;
	  t2_accum_v[v]  += t2_accum_vm1[v] *
	    poly_v.type2_gradient(xv, key_pv);
	  for (v2=0; v2<numVars; ++v2) {
	    if (v2 != v)
	      t2_accum_v[v2] += t2_accum_vm1[v2] * t1_grad;
	    t2_accum_vm1[v2]  = 0.;
	  }
	}
	else {
	  t1_accum_v[d]  += t1_accum_vm1[d] * t1_val;
	  t1_accum_vm1[d] = 0.;
	  t2_accum_v[v]  += t2_accum_vm1[v] * t2_val;
	  for (v2=0; v2<numVars; ++v2) {
	    if (v2 != v)
	      t2_accum_v[v2] += t2_accum_vm1[v2] * t1_val;
	    t2_accum_vm1[v2]  = 0.;
	  }
	}
      }
    }
    else { // case 4 with t1_val = 1., t1_grad = 0., t2_grad = 1
      t2_val = poly_v.type2_value(xv, key_pv);
      for (d=0; d<num_deriv_vars; ++d) {
	t2_accum_v   = t2_accumulators[d][v];
	t2_accum_vm1 = t2_accumulators[d][vm1];
	if (dvv[d] - 1 == v) {
	  t1_accum_v[d] = t1_accum_vm1[d] = 0.;
	  t2_accum_v[v] = t2_accum_vm1[v];
	  for (v2=0; v2<numVars; ++v2) {
	    if (v2 != v) t2_accum_v[v2] = 0.;
	    t2_accum_vm1[v2] = 0.;
	  }
	}
	else {
	  t1_accum_v[d] = t1_accum_vm1[d]; t1_accum_vm1[d] = 0.;
	  t2_accum_v[v] = t2_accum_vm1[v] * t2_val;
	  for (v2=0; v2<numVars; ++v2) {
	    if (v2 != v) t2_accum_v[v2] = t2_accum_vm1[v2];
	    t2_accum_vm1[v2] = 0.;
	  }
	}
      }
    }
    if (key_pv + 1 != poly_v.interpolation_size())
      break;
  }
}


/** Computes the specifics of a higher order grid for reinterpolating
    covariance over dimensions that will not be integrated. */
void SharedNodalInterpPolyApproxData::
reinterpolated_level(const UShortArray& lev_index)
{
  driverRep->reinterpolated_tensor_grid(lev_index, nonRandomIndices);
  update_tensor_interpolation_basis(driverRep->reinterpolated_level_index(),
				    nonRandomIndices);
}


/** Used for pre-computing the expected values of basis polynomial
    products within the PRODUCT_OF_INTERPOLANTS_FULL option for
    all_variables covariance. */
void SharedNodalInterpPolyApproxData::
update_nonzero_basis_products(const UShort2DArray& sm_multi_index)
{
  UShortArray max_lev_index(numVars, 0);
  size_t i, j, num_sm_mi = sm_multi_index.size();
  for (i=0; i<num_sm_mi; ++i) {
    const UShortArray& sm_mi_i = sm_multi_index[i];
    for (j=0; j<numVars; ++j)
      if (sm_mi_i[j] > max_lev_index[j])
	max_lev_index[j] = sm_mi_i[j];
  }

  // TO DO: ensure point mode is Gaussian quadrature, not nested rules
  const Real3DArray&    pts_1d = driverRep->collocation_points_1d();
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();

  // set up bookkeeping for random variables with unique interpolants
  bool empty_nz = nonZerosMapIndices.empty();
  SizetList::iterator it;
  if (empty_nz) {
    size_t num_v = randomIndices.size(), num_nz = 0, v1, v2, v1_cntr, v2_cntr;
    unsigned short min_max_lev;
    nonZerosMapIndices.resize(num_v); nonZerosMapIndices[0] = num_nz++;
    it=randomIndices.begin(); ++it; SizetList::iterator it2; bool found;
    for (v1_cntr=1; it!=randomIndices.end(); ++it, ++v1_cntr) {
      v1 = *it;
      for (it2=randomIndices.begin(), v2_cntr=0; it2!=it; ++it2, ++v2_cntr) {
	v2 = *it2; min_max_lev = std::min(max_lev_index[v1], max_lev_index[v2]);
	if (same_basis(min_max_lev, v1, v2)) {
	  nonZerosMapIndices[v1_cntr] = nonZerosMapIndices[v2_cntr];
	  found = true; break;
	}
      }
      if (!found)
	nonZerosMapIndices[v1_cntr] = num_nz++;
    }
    //nonZerosMapMaxLevels.resize(num_nz);
    nonZerosMapMaxLevels.assign(num_nz, 0);
    nonZerosMapArray.resize(num_nz);
  }
  InterpolationPolynomial *poly_rep_i, *poly_rep_j;
  UShort2DMultiSet map_key; UShortMultiSet mk1, mk2;
  size_t k, l, v, v_cntr = 0, num_pts_i, num_pts_j, nz_index;
  unsigned short j_start, max_lev_v, start; Real basis_prod_v;
  UShortMultiSet::iterator kit, lit; UShort2DMultiSet::iterator mit;
  for (it=randomIndices.begin(); it!=randomIndices.end(); ++it, ++v_cntr) {
    v = *it;
    max_lev_v = max_lev_index[v];
    nz_index  = nonZerosMapIndices[v_cntr];
    start     = nonZerosMapMaxLevels[nz_index] + 1;
    if (start <= max_lev_v) {
      const RealArray&    pts_1d_v =    pts_1d[max_lev_v][v];
      const RealArray& t1_wts_1d_v = t1_wts_1d[max_lev_v][v];
      UShort2DMultiSetRealMap& non_zeros_map = nonZerosMapArray[nz_index];
      for (i=1; i<=max_lev_v; ++i) {
	// InterpolationPolynomial is polynomialBasis[level][var]
	poly_rep_i = (InterpolationPolynomial*)
	  polynomialBasis[i][v].polynomial_rep();
	num_pts_i = poly_rep_i->interpolation_size();
	j_start = (i<start) ? start : 1;
	for (j=j_start; j<i; ++j) {
	  poly_rep_j = (InterpolationPolynomial*)
	    polynomialBasis[j][v].polynomial_rep();
	  num_pts_j = poly_rep_j->interpolation_size();
	  mk1.clear(); mk1.insert(num_pts_i);
	  mk2.clear(); mk2.insert(num_pts_j);
	  for (k=0; k<num_pts_i; ++k) {
	    kit = mk1.insert(k);
	    map_key.clear(); map_key.insert(mk1);
	    for (l=0; l<num_pts_j; ++l) {
	      if (basis_product_1d(poly_rep_i, poly_rep_j, k, l, pts_1d_v,
				   t1_wts_1d_v, basis_prod_v)) {
		lit = mk2.insert(l); mit = map_key.insert(mk2);
		non_zeros_map[map_key] = basis_prod_v;
		map_key.erase(mit); mk2.erase(lit);
#ifdef DEBUG
		PCout << "basis product var " << v << " lev_1 = " << i
		      << " lev_2 = " << j << " key_1 = " << k
		      << " key_2 = " << l << " basis_prod_v = "
		      << basis_prod_v << std::endl;
#endif // DEBUG
	      }
	    }
	    mk1.erase(kit);
	  }
	}
      }
      // Update max level tracking
      nonZerosMapMaxLevels[nz_index] = max_lev_v;
    }
  }
}


bool SharedNodalInterpPolyApproxData::
basis_product(const UShortArray& lev_index_1, const UShortArray& key_1,
	      const UShortArray& lev_index_2, const UShortArray& key_2,
	      Real& prod)
{
  SizetList::iterator it; size_t v, v_cntr = 0;
  unsigned short l1v, l2v, k1v, k2v;
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  prod = 1.;
  for (it=randomIndices.begin(); it!=randomIndices.end(); ++it, ++v_cntr) {
    v = *it;
    l1v = lev_index_1[v]; l2v = lev_index_2[v]; k1v = key_1[v]; k2v = key_2[v];
    // Embed knowledge of special cases:
    // > if one level index is 0 (single colloc pt for which L=1=constant), then
    //   expected val is normal weight from integration of single Lagrange poly
    // > if basis order is same, then E[L_i L_j] = delta_ij * weight_i by 1/0
    //   properties at Gauss points (this fact also used in
    //   tensor_product_covariance() function over the same tensor grid)
    // > if level indices are nonzero and different (for either nested or
    //   non-nested points), then numerical evaluation is needed.
    if (l1v == 0)
      prod *= t1_wts_1d[l2v][v][k2v];
    else if (l2v == 0)
      prod *= t1_wts_1d[l1v][v][k1v];
    else if (l1v == l2v) { // same level: bases eval as delta_ij using interpPts
      if (k1v == k2v) // 1D integral is 1D collocation weight
	prod *= t1_wts_1d[l1v][v][k1v];
      else            // 1D integral is zero
	return false;
    }
    else { // nonzero differing levels: lookup precomputed product
      UShort2DMultiSetRealMap& non_zeros_map
	= nonZerosMapArray[nonZerosMapIndices[v_cntr]];
      // lookup 1D terms
      UShortMultiSet mk1, mk2; mk1.insert(k1v); mk2.insert(k2v);
      mk1.insert(polynomialBasis[l1v][v].interpolation_size());
      mk2.insert(polynomialBasis[l2v][v].interpolation_size());
      UShort2DMultiSet map_key; map_key.insert(mk1); map_key.insert(mk2);
      UShort2DMultiSetRealMap::iterator it = non_zeros_map.find(map_key);
      if (it == non_zeros_map.end())
	return false; // zeros are not stored
      else
	prod *= it->second;
    }
  }
  return true;
}

} // namespace Pecos
