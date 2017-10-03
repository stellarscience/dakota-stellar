/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedInterpPolyApproxData
//- Description:  Implementation code for SharedInterpPolyApproxData class
//-               
//- Owner:        Mike Eldred

#include "SharedInterpPolyApproxData.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "TensorProductDriver.hpp"
#include "CombinedSparseGridDriver.hpp"
#include "HierarchSparseGridDriver.hpp"

namespace Pecos {


/** This version provides the polynomial types needed to retrieve
    collocation points and weights by an integration driver.  These
    may involve orthogonal polynomials which will differ from the
    interpolation polynomial types used in the basis. */
bool SharedInterpPolyApproxData::
initialize_driver_types_rules(const ShortArray& u_types,
			      const BasisConfigOptions& bc_options,
			      ShortArray& basis_types, ShortArray& colloc_rules)
{
  bool extra_dist_params = false;
  size_t i, num_vars = u_types.size();
  basis_types.resize(num_vars); colloc_rules.resize(num_vars);

  // Initialize basis_types, colloc_rules, and extra_dist_params from u_types.
  // interpolation has different requirements from integration/projection; i.e,
  // a good choice for interpolation nodes minimize the Lebesgue constant.  This
  // goal is reflected in the STD_UNIFORM case below.
  for (i=0; i<num_vars; ++i)
    switch (u_types[i]) {
    case STD_UNIFORM: // specialized for interpolation (min Lebesgue constant)
      if (bc_options.piecewiseBasis) {
	basis_types[i] = (bc_options.useDerivs) ?
	  PIECEWISE_CUBIC_INTERP : PIECEWISE_LINEAR_INTERP;
	if (bc_options.openRuleOverride) // closed nested rules required
	  PCerr << "Warning: open rules not currently supported for piecewise "
		<< "polynomial interpolants.  Ignoring override." << std::endl;
	colloc_rules[i] = (bc_options.equidistantRules) ?
	  NEWTON_COTES : CLENSHAW_CURTIS;
      }
      else if (bc_options.gaussRuleOverride) {
        basis_types[i]  = (bc_options.useDerivs) ?
          HERMITE_INTERP : LEGENDRE_ORTHOG;
        colloc_rules[i] = (bc_options.nestedRules) ?
          GAUSS_PATTERSON : GAUSS_LEGENDRE;
      }
      else {
	// LEGENDRE_ORTHOG provides more rule options than CHEBYSHEV_ORTHOG
	// if driver mode is changed (removing need to update basis_types)
	basis_types[i]  = (bc_options.useDerivs) ?
	  HERMITE_INTERP : LEGENDRE_ORTHOG;//CHEBYSHEV_ORTHOG;
	colloc_rules[i] = (bc_options.openRuleOverride) ?
	  FEJER2 : CLENSHAW_CURTIS;
      }
      break;
    default: // all other cases currently rely on Gaussian quadrature rules
      if (initialize_orthogonal_basis_type_rule(u_types[i], bc_options,
						basis_types[i],colloc_rules[i]))
	extra_dist_params = true;
      break;
    }

  return extra_dist_params;
}


void SharedInterpPolyApproxData::
initialize_polynomial_basis_type(short& poly_type_1d, short& rule)
{
  switch (basisType) {
  case PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    poly_type_1d = (basisConfigOptions.useDerivs) ?
      PIECEWISE_CUBIC_INTERP : PIECEWISE_LINEAR_INTERP;
    rule = NEWTON_COTES;                    break;
  case GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL:
  case GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    poly_type_1d = (basisConfigOptions.useDerivs) ?
      HERMITE_INTERP : LAGRANGE_INTERP;
    rule = NO_RULE;                         break;
  default:
    poly_type_1d = NO_POLY; rule = NO_RULE; break;
  }
}


void SharedInterpPolyApproxData::allocate_data()
{
  // use barycentric formulation for global Lagrange basis polynomials.
  // Note: flag needed below in update_{tensor,sparse}_interpolation_basis().
  barycentricFlag = ( !basisConfigOptions.useDerivs &&
    ( basisType == GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL ||
      basisType == GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL ) );

  bool param_update = false;
  const std::vector<BasisPolynomial>& num_int_poly_basis
    = driverRep->polynomial_basis();
  size_t i;
  for (i=0; i<numVars; ++i)
    if (num_int_poly_basis[i].parametric_update())
      { param_update = true; break; }

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray&   quad_order = tpq_driver->quadrature_order();

    // can't use quad_order > quadOrderPrev logic since only 1 pt set is stored
    bool change_order = (quad_order != quadOrderPrev);
    if (change_order || param_update)
      update_tensor_interpolation_basis(tpq_driver->level_index());
    if (change_order) {
      allocate_component_sobol();
      quadOrderPrev = quad_order;
    }
    break;
  }
  case COMBINED_SPARSE_GRID: case HIERARCHICAL_SPARSE_GRID: {
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    unsigned short    ssg_level  = ssg_driver->level();
    const RealVector& aniso_wts  = ssg_driver->anisotropic_weights();

    // Ignore weights since they only reduce the interpolation depth from the
    // level and the basis update uses a coarse increment based on level.  This
    // matches isotropic sparse grids, but forces fewer and larger updates in
    // the case of anisotropic or generalized grids.
    bool change_lev = (ssg_level != ssgLevelPrev), 
       increase_lev = (ssgLevelPrev == USHRT_MAX || ssg_level > ssgLevelPrev);
    if (increase_lev || param_update)
      update_sparse_interpolation_basis(ssg_level);
    if (change_lev) {
      allocate_component_sobol();
      ssgLevelPrev = ssg_level; //ssgAnisoWtsPrev = aniso_wts;
    }
    break;
  }
  }
}


void SharedInterpPolyApproxData::increment_data()
{
  unsigned short max_set_index = 0;
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: {
    // As for allocate_arrays(), increments are performed in coarser steps
    // than may be strictly necessary: all increments are filled in for all
    // vars for a step in level (ignoring anisotropy or generalized indices).
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    const UShortArray& trial_set = csg_driver->trial_set();
    for (size_t i=0; i<numVars; ++i)
      if (trial_set[i] > max_set_index)
	max_set_index = trial_set[i];
    break;
  }
  case HIERARCHICAL_SPARSE_GRID: {
    HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
    switch (expConfigOptions.refinementControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: { // generalized sparse grids
      const UShortArray& trial_set = hsg_driver->trial_set();
      for (size_t i=0; i<numVars; ++i)
	if (trial_set[i] > max_set_index)
	  max_set_index = trial_set[i];
      break;
    }
    default: { // isotropic/anisotropic refinement
      const UShort3DArray&   sm_mi = hsg_driver->smolyak_multi_index();
      const UShortArray& incr_sets = hsg_driver->increment_sets();
      size_t lev, num_lev = sm_mi.size(), set, start_set, num_sets, v;
      for (lev=0; lev<num_lev; ++lev) {
	start_set = incr_sets[lev]; num_sets = sm_mi[lev].size();
	for (set=start_set; set<num_sets; ++set) {
	  const UShortArray& sm_set = sm_mi[lev][set];
	  for (v=0; v<numVars; ++v)
	    if (sm_set[v] > max_set_index)
	      max_set_index = sm_set[v];
	}
      }
      break;
    }
    }
    break;
  }
  default:
    PCerr << "Error: unsupported grid definition in InterpPolyApproximation::"
	  << "increment_coefficients()" << std::endl;
    abort_handler(-1);
    break;
  }

  update_sparse_interpolation_basis(max_set_index);
  increment_component_sobol();
}


void SharedInterpPolyApproxData::decrement_data()
{
  // leave polynomialBasis as is

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: case HIERARCHICAL_SPARSE_GRID: {
    // move previous expansion data to current expansion
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    poppedLevMultiIndex.push_back(ssg_driver->trial_set());
    break;
  }
  }
}


void SharedInterpPolyApproxData::post_push_data()
{
  // leave polynomialBasis as is (a previous increment is being restored)

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: case HIERARCHICAL_SPARSE_GRID: {
    // move previous expansion data to current expansion
    SparseGridDriver* ssg_driver = (SparseGridDriver*)driverRep;
    std::deque<UShortArray>::iterator sit
      = std::find(poppedLevMultiIndex.begin(), poppedLevMultiIndex.end(),
		  ssg_driver->trial_set());
    if (sit != poppedLevMultiIndex.end())
      poppedLevMultiIndex.erase(sit);
    break;
  }
  }
}


void SharedInterpPolyApproxData::post_finalize_data()
{
  // leave polynomialBasis as is (all previous increments are being restored)

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: case HIERARCHICAL_SPARSE_GRID:
    // move previous expansion data to current expansion
    poppedLevMultiIndex.clear();
    break;
  }
}


size_t SharedInterpPolyApproxData::pre_combine_data(short combine_type)
{
  // Sufficient for two grids: if not currently the maximal grid, then swap
  // with the stored grid (only one is stored)
  //bool swap = !driverRep->maximal_grid();
  //if (swap) { driverRep->swap_grid(); allocate_component_sobol(); }

  // For open-ended number of stored grids: retrieve the most refined from the
  // existing grids (from sequence specification + any subsequent refinement).
  // Note: if we assume that multiIndex subsets are enforced across a hierarchy,
  // then the maximal grid is sufficient to allow reinterpolation of all data.
  size_t max_index = driverRep->maximal_grid();
  if (max_index != _NPOS)
    { driverRep->swap_grid(max_index); allocate_component_sobol(); }

  // Most general: overlay all grid refinement levels to create a new superset:
  //size_t new_index = driverRep->overlay_maximal_grid();
  //if (current_grid_index() != new_index) driverRep->swap_grid(new_index);

  return max_index;
}


void SharedInterpPolyApproxData::post_combine_data(short combine_type)
{ driverRep->clear_stored(); }


void SharedInterpPolyApproxData::
update_tensor_interpolation_basis(const UShortArray& lev_index)
{
  // resize if needed (leaving previous levels unmodified)
  resize_polynomial_basis(lev_index);

  // fill any required gaps in polynomialBasis
  for (size_t i=0; i<numVars; ++i)
    update_interpolation_basis(lev_index[i], i);
}


void SharedInterpPolyApproxData::
update_tensor_interpolation_basis(const UShortArray& lev_index,
				  const SizetList& subset_indices)
{
  // resize if needed (leaving previous levels unmodified)
  resize_polynomial_basis(lev_index);

  // fill any required gaps in polynomialBasis
  SizetList::const_iterator cit; size_t i;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit)
    { i = *cit; update_interpolation_basis(lev_index[i], i); }
}


void SharedInterpPolyApproxData::
update_sparse_interpolation_basis(unsigned short max_level)
{
  // resize if needed (leaving previous levels unmodified)
  // j range is 0:w inclusive; i range is 1:w+1 inclusive
  size_t l, v;//, orig_size = polynomialBasis.size();
  resize_polynomial_basis(max_level);

  const std::vector<BasisPolynomial>& num_int_poly_basis
    = driverRep->polynomial_basis();
  // We must currently process all levels, even if not parameterized, since
  // update_interpolation_basis() can only update the polynomial basis if the
  // collocation points are available and the collocation points are only 
  // available if a particular index level is visited for that variable within 
  // IntegrationDriver::compute_tensor_grid() (which calls IntegrationDriver::
  // update_1d_collocation_points_weights()).  That is, there may be gaps that
  // need to be filled in levels that were previously allocated.
  for (v=0; v<numVars; ++v) {
    //if (num_int_poly_basis[v].parameterized()) // check all levels for updates
      for (l=0; l<=max_level; ++l)
	update_interpolation_basis(l, v);
    //else                                       // update only the new levels
    //  for (l=orig_size; l<=max_level; ++l)
    //    update_interpolation_basis(l, v);
  }
}


void SharedInterpPolyApproxData::
update_interpolation_basis(unsigned short lev_index, size_t var_index)
{
  // fill gaps that may exist within any level
  const RealArray& colloc_pts_1d_lv
    = driverRep->collocation_points_1d()[lev_index][var_index];
  if (!colloc_pts_1d_lv.empty()) {
    const BasisPolynomial& num_int_poly_basis_v
      = driverRep->polynomial_basis()[var_index];
    std::vector<BasisPolynomial>& poly_basis_l  = polynomialBasis[lev_index];
    BasisPolynomial&              poly_basis_lv = poly_basis_l[var_index];
    short poly_type_1d, rule;
    // don't share reps in case of parameterized basis or barycentric interp,
    // due to need for individual updates to parameters or interpolated x value.
    if (barycentricFlag || num_int_poly_basis_v.parameterized()) {
      if (poly_basis_lv.is_null()) {
	initialize_polynomial_basis_type(poly_type_1d, rule);
	poly_basis_lv = BasisPolynomial(poly_type_1d, rule);
	poly_basis_lv.interpolation_points(colloc_pts_1d_lv);
      }
      else if (num_int_poly_basis_v.parametric_update())
	poly_basis_lv.interpolation_points(colloc_pts_1d_lv);
    }
    else if (poly_basis_lv.is_null()) { // can share reps for efficiency
      size_t var_index2;
      if (find_basis(lev_index, var_index, var_index2))
	poly_basis_lv = poly_basis_l[var_index2]; // reuse prev via shared rep
      else { // instantiate and initialize a new unique instance
	initialize_polynomial_basis_type(poly_type_1d, rule);
	poly_basis_lv = BasisPolynomial(poly_type_1d, rule);
	poly_basis_lv.interpolation_points(colloc_pts_1d_lv);
      }
    }
  }
}


bool SharedInterpPolyApproxData::
find_basis(unsigned short level, size_t v1, size_t& v2)
{
  std::vector<BasisPolynomial>& poly_basis_l = polynomialBasis[level];
  for (v2=0; v2<numVars; ++v2)
    if (v2 != v1 && !poly_basis_l[v2].is_null() && same_basis(level, v1, v2))
      return true;
  return false; 
}


bool SharedInterpPolyApproxData::
same_basis(unsigned short level, size_t v1, size_t v2)
{
  const ShortArray& rules = driverRep->collocation_rules();
  short rule1 = rules[v1];
  if (rules[v2] == rule1)
    switch (rule1) {
    case GAUSS_JACOBI: case GEN_GAUSS_LAGUERRE: case GOLUB_WELSCH: {
      // rule type insufficient in these cases, check collocation points
      const Real2DArray& colloc_pts_1d
	= driverRep->collocation_points_1d()[level];
      return (colloc_pts_1d[v1] == colloc_pts_1d[v2]); break;
    }
    default:
      return true;                                     break;
    }
  else
    return false;
}


/** Barycentric approach is only valid for value-based global Lagrange
    interpolation, either nodal or hierarchical.  General approach is
    valid for value-based or gradient-enhanced, local or global, and
    nodal or hierarchical. */
Real SharedInterpPolyApproxData::
tensor_product_value(const RealVector& x, const RealVector& exp_t1_coeffs,
		     const RealMatrix& exp_t2_coeffs,
		     const UShortArray& basis_index, const UShort2DArray& key,
		     const SizetArray& colloc_index)
{
  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 1); // value factors needed

    size_t j, num_act_v = barycentric_active_variables(basis_index);
    if (num_act_v == 0) { // convert 1-D exact indices into n-D colloc index
      size_t pt_index = barycentric_exact_index(basis_index);
      if (pt_index == _NPOS) // if exactIndex but not exactDeltaIndex, then
	return 0.;           // at least one dimension has value factor = 0.
      else
	return (colloc_index.empty()) ?
	  exp_t1_coeffs[pt_index] : exp_t1_coeffs[colloc_index[pt_index]];
    }
    else if (num_act_v == numVars) { // interpolation over all variables
      RealVector accumulator(numVars); // init to 0.
      precompute_max_keys(basis_index); // if needed for efficiency
      unsigned short key_i0, key_ij, bi_j, bi_0 = basis_index[0],
	max0 = tensor_product_max_key(0, bi_0);
      const RealVector& bc_vf_0
	= polynomialBasis[bi_0][0].barycentric_value_factors();
      size_t i, num_colloc_pts = key.size();//,ei_0 = poly_0.exact_index(),ei_j;
      for (i=0; i<num_colloc_pts; ++i) {
	const UShortArray& key_i = key[i]; key_i0 = key_i[0];
	//if (ei_0 == _NPOS)
	accumulator[0] += (colloc_index.empty()) ?
	  exp_t1_coeffs[i]               * bc_vf_0[key_i0] :
	  exp_t1_coeffs[colloc_index[i]] * bc_vf_0[key_i0];
	//else if (ei_0 == key_i0)                  // only 1 pt should match
	//  accumulator[0]  = (colloc_index.empty()) ? // (= instead of +=)
	//    exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
	if (key_i0 == max0) {
	  // accumulate sums over variables with max key value
	  for (j=1; j<numVars; ++j) {
	    key_ij = key_i[j]; bi_j = basis_index[j];
	    //ei_j = poly_j.exact_index();
	    //if (ei_j == _NPOS)
	    accumulator[j] += accumulator[j-1] *
	      polynomialBasis[bi_j][j].barycentric_value_factor(key_ij);
	    //else if (ei_j == key_ij)              // only 1 pt should match
	    //  accumulator[j] = accumulator[j-1]; // (= instead of +=)
	    accumulator[j-1] = 0.;
	    if (key_ij != tensor_product_max_key(j, bi_j))
	      break;
	  }
	}
      }
      return accumulator[numVars-1] / barycentric_value_scaling(basis_index);
    }
    else { // partial interpolation over active variables
      SizetList pt_factors, act_v_set; //UShortList num_keys;
      size_t num_act_pts, pt_index;
      barycentric_partial_indexing(basis_index, /*num_keys,*/ pt_factors,
				   act_v_set, num_act_pts, pt_index);
      if (pt_index == _NPOS) return 0.;
      RealVector accumulator(num_act_v); // init to 0.
      accumulate_barycentric_partial(exp_t1_coeffs, basis_index, key,
				     colloc_index, /*num_keys,*/ pt_factors,
				     act_v_set, num_act_pts, pt_index,
				     accumulator);
      return accumulator[num_act_v-1] / barycentric_value_scaling(basis_index);
    }

    /*
    // Simple bc option: delegate loops; cleaner, but some efficiency lost (and
    // precision as well due to subtractive cancellation among large products)
    size_t i, num_colloc_pts = key.size(); Real tp_val = 0.;
    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[i] *
	  barycentric_value_factor(key[i], basis_index);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[colloc_index[i]] *
	  barycentric_value_factor(key[i], basis_index);
    // apply barycentric denominator
    return tp_val / barycentric_value_scaling(basis_index);
    */
  }
  else if (exp_t2_coeffs.empty()) {
    // Horner's rule approach:
    RealVector accumulator(numVars); // init to 0.
    precompute_max_keys(basis_index); // if needed for efficiency
    unsigned short bi_0 = basis_index[0], bi_j, key_i0, key_ij,
      max0 = tensor_product_max_key(0, bi_0);
    BasisPolynomial& poly_0 = polynomialBasis[bi_0][0];
    size_t i, j, num_colloc_pts = key.size(); Real x0 = x[0];
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      if (bi_0)
	accumulator[0] += (colloc_index.empty()) ?
	  exp_t1_coeffs[i]               * poly_0.type1_value(x0, key_i0) :
	  exp_t1_coeffs[colloc_index[i]] * poly_0.type1_value(x0, key_i0);
      else
	accumulator[0] = (colloc_index.empty()) ?
	  exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      if (key_i0 == max0) {
	// accumulate sums over variables with max key value
	for (j=1; j<numVars; ++j) {
	  key_ij = key_i[j]; bi_j = basis_index[j];
	  if (bi_j)
	    accumulator[j] += accumulator[j-1] *
	      polynomialBasis[bi_j][j].type1_value(x[j], key_ij);
	  else
	    accumulator[j]  = accumulator[j-1];
	  accumulator[j-1] = 0.;
	  if (key_ij != tensor_product_max_key(j, bi_j))
	    break;
	}
      }
    }
    return accumulator[numVars-1];

    /*
    // Simpler but more expensive approach:
    size_t i, num_colloc_pts = key.size(); Real tp_val = 0.;
    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[i] * 
	          type1_interpolant_value(x, key[i], basis_index);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += exp_t1_coeffs[colloc_index[i]] *
	          type1_interpolant_value(x, key[i], basis_index);
    return tp_val;
    */
  }
  else {
    // Horner's rule approach:
    RealVector t1_accumulator(numVars);          // init to 0.
    RealMatrix t2_accumulator(numVars, numVars); // init to 0.
    Real *t2_accum_0 = t2_accumulator[0], *t2_accum_j, *t2_accum_jm1;
    precompute_max_keys(basis_index); // if needed for efficiency
    unsigned short     bi_0 = basis_index[0], bi_j, key_i0, key_ij,
      max0 = tensor_product_max_key(0, bi_0);
    BasisPolynomial& poly_0 = polynomialBasis[bi_0][0];
    size_t i, j, jm1, k, c_index, num_colloc_pts = key.size();
    Real t1_val, x0 = x[0];
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      const Real* t2_coeffs_i = exp_t2_coeffs[c_index];
      if (bi_0 == 0) {
	t1_accumulator[0]  = exp_t1_coeffs[c_index];  // t1 value is 1
	t2_accum_0[0]      = t2_coeffs_i[0] * poly_0.type2_value(x0, key_i0);
	for (k=1; k<numVars; ++k)
	  t2_accum_0[k]    = t2_coeffs_i[k];                 // t1 value is 1
      }
      else {
	t1_val = poly_0.type1_value(x0, key_i0);
	t1_accumulator[0] += exp_t1_coeffs[c_index] * t1_val;
	t2_accum_0[0]     += t2_coeffs_i[0] * poly_0.type2_value(x0, key_i0);
	for (k=1; k<numVars; ++k)
	  t2_accum_0[k]   += t2_coeffs_i[k] * t1_val;
      }
      if (key_i0 == max0) {
	// accumulate sums over variables with max key value
	for (j=1; j<numVars; ++j) {
	  bi_j = basis_index[j]; key_ij = key_i[j]; jm1 = j-1;
	  t2_accum_j = t2_accumulator[j]; t2_accum_jm1 = t2_accumulator[jm1];
	  if (bi_j == 0) {
	    t1_accumulator[j] = t1_accumulator[jm1];         // t1 value is 1
	    t2_accum_j[j] = t2_accum_jm1[j] *
	      polynomialBasis[bi_j][j].type2_value(x[j], key_ij);
	    for (k=0; k<numVars; ++k)
	      if (k != j)
		t2_accum_j[k] = t2_accum_jm1[k];             // t1 value is 1
	  }
	  else {
	    BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
	    t1_val = poly_j.type1_value(x[j], key_ij);
	    t1_accumulator[j] += t1_accumulator[jm1] * t1_val;
	    t2_accum_j[j] += t2_accum_jm1[j] * poly_j.type2_value(x[j], key_ij);
	    for (k=0; k<numVars; ++k)
	      if (k != j)
		t2_accum_j[k] += t2_accum_jm1[k] * t1_val;
	  }
	  t1_accumulator[jm1] = 0.;
	  for (k=0; k<numVars; ++k)
	    t2_accum_jm1[k] = 0.;
	  if (key_ij != tensor_product_max_key(j, bi_j))
	    break;
	}
      }
    }
    Real  tp_val   = t1_accumulator[numVars-1];
    Real* t2_accum = t2_accumulator[numVars-1];
    for (j=0; j<numVars; ++j)
      tp_val += t2_accum[j];
    return tp_val;

    /*
    // Simpler but more expensive approach:
    size_t i, num_colloc_pts = key.size(), j, c_index; Real tp_val = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      tp_val += exp_t1_coeffs[c_index] *
	        type1_interpolant_value(x, key_i, basis_index);
      const Real* exp_t2_coeff_i = exp_t2_coeffs[c_index];
      for (j=0; j<numVars; ++j)
	tp_val += exp_t2_coeff_i[j] *
	          type2_interpolant_value(x, j, key_i, basis_index);
    }
    return tp_val;
    */
  }
}


/** All variables version. */
Real SharedInterpPolyApproxData::
tensor_product_value(const RealVector& x, const RealVector& subset_t1_coeffs,
		     const RealMatrix& subset_t2_coeffs,
		     const UShortArray& basis_index,
		     const UShort2DArray& subset_key,
		     const SizetArray& subset_colloc_index,
		     const SizetList& subset_indices)
{
  // Note: subset_* are consistent with the reduced variable subset,
  // but x and basis_index are full space.

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, subset_indices, 1); // value factors needed

    size_t num_subset_v = subset_indices.size(),
      num_act_v = barycentric_active_variables(basis_index, subset_indices);
    if (num_act_v == 0) { // convert 1-D exact indices into n-D colloc index
      size_t pt_index = barycentric_exact_index(basis_index, subset_indices);
      // Note: pt_index calculation utilizes only the subset variables
      if (pt_index == _NPOS) // if exactIndex but not exactDeltaIndex, then
	return 0.;           // at least one dimension has value factor = 0.
      else
	return (subset_colloc_index.empty()) ? subset_t1_coeffs[pt_index] :
	  subset_t1_coeffs[subset_colloc_index[pt_index]];
    }
    else if (num_act_v == num_subset_v) { // interpolation over all of subset
      SizetList::const_iterator cit = subset_indices.begin();
      size_t i, j, num_colloc_pts = subset_key.size(), v0 = *cit, vj;
      RealVector accumulator(num_subset_v); // init to 0.
      precompute_max_keys(basis_index, subset_indices); // if needed
      unsigned short bi_0 = basis_index[v0], bi_j, key_i0, key_ij,
	max0 = tensor_product_max_key(v0, bi_0);
      const RealVector& bc_vf_0
	= polynomialBasis[bi_0][v0].barycentric_value_factors();
      for (i=0; i<num_colloc_pts; ++i) {
	// Note: key is reduced in first dimension (pts), but not second (vars)
	const UShortArray& key_i = subset_key[i]; key_i0 = key_i[v0];
	accumulator[0] += (subset_colloc_index.empty()) ?
	  subset_t1_coeffs[i]                      * bc_vf_0[key_i0] :
	  subset_t1_coeffs[subset_colloc_index[i]] * bc_vf_0[key_i0];
	if (key_i0 == max0) {
	  // accumulate sums over variables with max key value
	  for (j=1, cit=++subset_indices.begin(); j<num_subset_v; ++j, ++cit) {
	    vj = *cit; key_ij = key_i[vj]; bi_j = basis_index[vj];
	    accumulator[j] += accumulator[j-1] *
	      polynomialBasis[bi_j][vj].barycentric_value_factor(key_ij);
	    accumulator[j-1] = 0.;
	    if (key_ij != tensor_product_max_key(vj, bi_j))
	      break;
	  }
	}
      }
      return accumulator[num_subset_v-1]
	/ barycentric_value_scaling(basis_index, subset_indices);
    }
    else { // partial interpolation over active variables
      SizetList pt_factors, act_v_set; //UShortList num_keys;
      size_t num_act_pts, pt_index;
      barycentric_partial_indexing(basis_index, subset_indices, //num_keys,
				   pt_factors, act_v_set, num_act_pts,pt_index);
      if (pt_index == _NPOS) return 0.;
      RealVector accumulator(num_act_v); // init to 0.
      accumulate_barycentric_partial(subset_t1_coeffs, basis_index, subset_key,
				     subset_colloc_index, //num_keys,
				     pt_factors, act_v_set, num_act_pts,
				     pt_index, accumulator);
      return accumulator[num_act_v-1]
	/ barycentric_value_scaling(basis_index, subset_indices);
    }

    /*
    // Simple bc option: delegate loops; cleaner, but less efficient/precise
    size_t i, num_colloc_pts = subset_key.size(); Real tp_val = 0.;
    if (subset_colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += subset_t1_coeffs[i] * 
	  barycentric_value_factor(subset_key[i], basis_index, subset_indices);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += subset_t1_coeffs[subset_colloc_index[i]] *
	  barycentric_value_factor(subset_key[i], basis_index, subset_indices);
    // apply barycentric denominator
    return tp_val / barycentric_value_scaling(basis_index, subset_indices);
    */
  }
  else if (subset_t2_coeffs.empty()) {
    size_t i, num_colloc_pts = subset_key.size(); Real tp_val = 0.;
    if (subset_colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += subset_t1_coeffs[i] * 
	  type1_interpolant_value(x, subset_key[i], basis_index,subset_indices);
    else
      for (i=0; i<num_colloc_pts; ++i)
	tp_val += subset_t1_coeffs[subset_colloc_index[i]] *
	  type1_interpolant_value(x, subset_key[i], basis_index,subset_indices);
    return tp_val;
  }
  else {
    size_t i, j, num_colloc_pts = subset_key.size(), c_index; Real tp_val = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = subset_key[i];
      c_index = (subset_colloc_index.empty()) ? i : subset_colloc_index[i];
      tp_val += subset_t1_coeffs[c_index] *
	type1_interpolant_value(x, key_i, basis_index, subset_indices);
      const Real* subset_t2_coeff_i = subset_t2_coeffs[c_index];
      for (j=0; j<numVars; ++j)
	tp_val += subset_t2_coeff_i[j] *
	  type2_interpolant_value(x, j, key_i, basis_index, subset_indices);
    }
    return tp_val;
  }
}


const RealVector& SharedInterpPolyApproxData::
tensor_product_gradient_basis_variables(const RealVector& x,
					const RealVector& exp_t1_coeffs,
					const RealMatrix& exp_t2_coeffs,
					const UShortArray& basis_index,
					const UShort2DArray& key,
					const SizetArray& colloc_index)
{
  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  size_t i, j, num_colloc_pts = key.size();

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 3); // value+gradient factors needed

    // Note: for basis gradients, one cannot eliminate exact index dimensions
    // as is managed in the barycentric case of tensor_product_value()

    precompute_max_keys(basis_index); // if needed for efficiency
    unsigned short key_i0, key_ij, bi_0 = basis_index[0], bi_j,
      max0 = tensor_product_max_key(0, bi_0);
    BasisPolynomial&   poly_0 = polynomialBasis[bi_0][0];
    const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
    const RealVector& bc_gf_0 = poly_0.barycentric_gradient_factors();
    size_t ei_0 = poly_0.exact_index();
    RealMatrix accumulator(numVars, numVars); // init to 0.
    Real *accum_0 = accumulator[0], t1_coeff;
    for (i=0; i<num_colloc_pts; ++i) {
      t1_coeff = (colloc_index.empty()) ? exp_t1_coeffs[i] :
	exp_t1_coeffs[colloc_index[i]];
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      accumulate_barycentric_gradient(bi_0, key_i0, ei_0, accum_0, t1_coeff,
				      bc_vf_0, bc_gf_0);
      if (key_i0 == max0) {
	// accumulate sums over variables with max key value
	for (j=1; j<numVars; ++j) {
	  bi_j = basis_index[j]; key_ij = key_i[j];
	  accumulate_barycentric_gradient(j, bi_j, key_ij,
	    polynomialBasis[bi_j][j], accumulator);
	  if (key_ij != tensor_product_max_key(j, bi_j))
	    break;
	}
      }
    }
    Real bcg_scale = barycentric_gradient_scaling(basis_index);
    Real* accum = accumulator[numVars-1];
    for (j=0; j<numVars; ++j)
      tpGradient[j] = accum[j] * bcg_scale;

    /*
    size_t i, j, num_colloc_pts = key.size();
    tpGradient = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += t1_coeff_i *
	  barycentric_gradient_factor(j, key_i, basis_index);
    }
    tpGradient.scale(barycentric_gradient_scaling(basis_index));
    */
  }
  else if (exp_t2_coeffs.empty()) {
    tpGradient = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index);
    }
  }
  else {
    size_t k;
    tpGradient = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      const Real* exp_t2_coeff_i = (colloc_index.empty()) ?
	exp_t2_coeffs[i] : exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<numVars; ++j) { // ith contribution to jth grad component
	tpGradient[j] += t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index);
	for (k=0; k<numVars; ++k) // type2 interpolant for kth grad comp
	  tpGradient[j] += exp_t2_coeff_i[k] *
	    type2_interpolant_gradient(x, j, k, key_i, basis_index);
      }
    }
  }
  return tpGradient;
}


const RealVector& SharedInterpPolyApproxData::
tensor_product_gradient_basis_variables(const RealVector& x,
					const RealVector& subset_t1_coeffs,
					const RealMatrix& subset_t2_coeffs,
					const UShortArray& basis_index,
					const UShort2DArray& subset_key,
					const SizetArray& subset_colloc_index,
					const SizetList&  subset_indices)
{
  if (tpGradient.length() != numVars)
    tpGradient.sizeUninitialized(numVars);
  size_t i, j, num_colloc_pts = subset_key.size();

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, subset_indices, 3);//value+grad factors needed

    // Note: for basis gradients, one cannot eliminate exact index dimensions
    // as is managed in the barycentric case of tensor_product_value()

    SizetList::const_iterator cit = subset_indices.begin();
    size_t num_subset_v = subset_indices.size(), v0 = *cit, vj;
    unsigned short key_i0, key_ij, bi_0 = basis_index[v0], bi_j;
    BasisPolynomial&   poly_v0 = polynomialBasis[bi_0][v0];
    const RealVector& bc_vf_v0 = poly_v0.barycentric_value_factors();
    const RealVector& bc_gf_v0 = poly_v0.barycentric_gradient_factors();
    size_t ei_0 = poly_v0.exact_index();
    precompute_max_keys(basis_index, subset_indices);// if needed for efficiency
    unsigned short max0 = tensor_product_max_key(v0, bi_0);
    RealMatrix accumulator(numVars, num_subset_v); // init to 0.
    Real *accum_0 = accumulator[0], t1_coeff;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = subset_key[i]; key_i0 = key_i[v0];
      t1_coeff = (subset_colloc_index.empty()) ? subset_t1_coeffs[i] :
	subset_t1_coeffs[subset_colloc_index[i]];
      accumulate_barycentric_gradient(bi_0, key_i0, ei_0, accum_0, t1_coeff,
				      bc_vf_v0, bc_gf_v0);
      if (key_i0 == max0) {
	// accumulate sums over variables with max key value
	for (j=1, cit=++subset_indices.begin(); j<num_subset_v; ++j, ++cit) {
	  vj = *cit; bi_j = basis_index[vj]; key_ij = key_i[vj];
	  accumulate_barycentric_gradient(j, bi_j, key_ij,
					  polynomialBasis[bi_j][vj],
					  accumulator);
	  if (key_ij != tensor_product_max_key(vj, bi_j))
	    break;
	}
      }
    }
    Real bcg_scale = barycentric_gradient_scaling(basis_index, subset_indices);
    Real* accum = accumulator[num_subset_v-1];
    for (j=0; j<numVars; ++j)
      tpGradient[j] = accum[j] * bcg_scale;

    /*
    tpGradient = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = subset_key[i];
      Real t1_coeff_i = (subset_colloc_index.empty()) ?
	subset_t1_coeffs[i] : subset_t1_coeffs[subset_colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += t1_coeff_i *
	  barycentric_gradient_factor(j, key_i, basis_index, subset_indices);
    }
    tpGradient.scale(barycentric_gradient_scaling(basis_index, subset_indices));
    */
  }
  else if (subset_t2_coeffs.empty()) {
    tpGradient = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = subset_key[i];
      Real t1_coeff_i = (subset_colloc_index.empty()) ?
	subset_t1_coeffs[i] : subset_t1_coeffs[subset_colloc_index[i]];
      for (j=0; j<numVars; ++j)
	tpGradient[j] += t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index, subset_indices);
    }
  }
  else {
    size_t k;
    tpGradient = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = subset_key[i];
      Real t1_coeff_i = (subset_colloc_index.empty()) ?
	subset_t1_coeffs[i] : subset_t1_coeffs[subset_colloc_index[i]];
      const Real* t2_coeff_i = (subset_colloc_index.empty()) ?
	subset_t2_coeffs[i] : subset_t2_coeffs[subset_colloc_index[i]];
      for (j=0; j<numVars; ++j) { // ith contribution to jth grad component
	tpGradient[j] += t1_coeff_i *
	  type1_interpolant_gradient(x, j, key_i, basis_index, subset_indices);
	for (k=0; k<numVars; ++k) // type2 interpolant for kth grad comp
	  tpGradient[j] += t2_coeff_i[k] *
	    type2_interpolant_gradient(x, j, k, key_i, basis_index,
				       subset_indices);
      }
    }
  }
  return tpGradient;
}


const RealVector& SharedInterpPolyApproxData::
tensor_product_gradient_basis_variables(const RealVector& x,
					const RealVector& exp_t1_coeffs,
					const RealMatrix& exp_t2_coeffs,
					const UShortArray& basis_index,
					const UShort2DArray& key,
					const SizetArray& colloc_index,
					const SizetArray& dvv)
{
  size_t i, j, deriv_index, num_colloc_pts = key.size(),
    num_deriv_vars = dvv.size();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);
  if (!num_deriv_vars) return tpGradient;

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 3); // value+gradient factors needed

    // Note: for basis gradients, one cannot eliminate exact index dimensions
    // as is managed in the barycentric case of tensor_product_value()

    precompute_max_keys(basis_index); // if needed for efficiency
    unsigned short key_i0, key_ij, bi_0 = basis_index[0], bi_j,
      max0 = tensor_product_max_key(0, bi_0);
    BasisPolynomial&   poly_0 = polynomialBasis[bi_0][0];
    const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
    const RealVector& bc_gf_0 = poly_0.barycentric_gradient_factors();
    size_t k, ei_0 = poly_0.exact_index(), ei_j, d0_index = dvv[0] - 1, start;
    RealMatrix accumulator(num_deriv_vars, numVars); // init to 0.
    Real *accum_0 = accumulator[0], *accum_j, *accum_jm1,
      t1_coeff_bc_vf_0, bc_vf_j, bc_gf_j, t1_coeff;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      t1_coeff = (colloc_index.empty()) ? exp_t1_coeffs[i] :
	exp_t1_coeffs[colloc_index[i]];
      // mirror version w/o DVV, only logging a subset of accum components
      start = 0;
      if (bi_0) {
	if (d0_index == 0)
	  { start = 1; accum_0[0] += t1_coeff * bc_gf_0[key_i0]; }
        if (ei_0 == _NPOS) {
	  t1_coeff_bc_vf_0 = t1_coeff * bc_vf_0[key_i0];
	  for (j=start; j<num_deriv_vars; ++j)
	    accum_0[j] += t1_coeff_bc_vf_0;
	}
	else if (ei_0 == key_i0) // value factor is 1, else value factor is 0
	  for (j=start; j<num_deriv_vars; ++j)
	    accum_0[j] += t1_coeff;
      }
      else { // grad factor is 0., value factor is omitted
	if (d0_index == 0) start = 1;
	for (j=start; j<num_deriv_vars; ++j)
	  accum_0[j] += t1_coeff;
      }
      if (key_i0 == max0) {
	// accumulate sums over variables with max key value
	for (j=1; j<numVars; ++j) {
	  bi_j = basis_index[j]; key_ij = key_i[j];
	  accum_j = accumulator[j]; accum_jm1 = accumulator[j-1];
	  BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
	  if (bi_j) {
	    ei_j = poly_j.exact_index();
	    bc_gf_j = poly_j.barycentric_gradient_factor(key_ij);
	    if (ei_j == _NPOS) { // bc_vf_j has general value
	      bc_vf_j = poly_j.barycentric_value_factor(key_ij);
	      for (k=0; k<num_deriv_vars; ++k) {
		accum_j[k] += (j == dvv[k] - 1) ? accum_jm1[k] * bc_gf_j
		                                : accum_jm1[k] * bc_vf_j;
		accum_jm1[k] = 0.;
	      }
	    }
	    else if (ei_j == key_ij) // bc_vf_j is 1
	      for (k=0; k<num_deriv_vars; ++k) {
		accum_j[k] += (j == dvv[k] - 1) ? accum_jm1[k] * bc_gf_j
		                                : accum_jm1[k];
		accum_jm1[k] = 0.;
	      }
	    else // bc_vf_j is 0
	      for (k=0; k<num_deriv_vars; ++k) {
		if (j == dvv[k] - 1) accum_j[k] += accum_jm1[k] * bc_gf_j;
		accum_jm1[k] = 0.;
	      }
	  }
	  else { // grad factor is zero, value factor is omitted
	    for (k=0; k<num_deriv_vars; ++k) {
	      if (j != dvv[k] - 1) accum_j[k] += accum_jm1[k];
	      accum_jm1[k] = 0.;
	    }
	  }
	  if (key_ij != tensor_product_max_key(j, bi_j))
	    break;
	}
      }
    }
    Real bcg_scale = barycentric_gradient_scaling(basis_index);
    Real* accum = accumulator[numVars-1];
    for (j=0; j<num_deriv_vars; ++j)
      tpGradient[j] = accum[j] * bcg_scale;

    /*
    tpGradient = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += t1_coeff_i *
	  barycentric_gradient_factor(deriv_index, key_i, basis_index);
      }
    }
    tpGradient.scale(barycentric_gradient_scaling(basis_index));
    */
  }
  else if (exp_t2_coeffs.empty()) {
    tpGradient = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += t1_coeff_i *
	  type1_interpolant_gradient(x, deriv_index, key_i, basis_index);
      }
    }
  }
  else {
    size_t k;
    tpGradient = 0.;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray&   key_i = key[i];
      Real t1_coeff_i = (colloc_index.empty()) ?
	exp_t1_coeffs[i] : exp_t1_coeffs[colloc_index[i]];
      const Real* exp_t2_coeff_i = (colloc_index.empty()) ?
	exp_t2_coeffs[i] : exp_t2_coeffs[colloc_index[i]];
      for (j=0; j<num_deriv_vars; ++j) {
	deriv_index = dvv[j] - 1; // requires an "All" view
	tpGradient[j] += t1_coeff_i *
	  type1_interpolant_gradient(x, deriv_index, key_i, basis_index);
	for (k=0; k<numVars; ++k)
	  tpGradient[j] += exp_t2_coeff_i[k] *
	    type2_interpolant_gradient(x, deriv_index, k, key_i, basis_index);
      }
    }
  }
  return tpGradient;
}


const RealVector& SharedInterpPolyApproxData::
tensor_product_gradient_nonbasis_variables(const RealVector& x,
					   const RealMatrix& exp_t1_coeff_grads,
					   const UShortArray& basis_index,
					   const UShort2DArray& key,
					   const SizetArray& colloc_index)
{
  size_t num_deriv_vars = exp_t1_coeff_grads.numRows();
  if (tpGradient.length() != num_deriv_vars)
    tpGradient.sizeUninitialized(num_deriv_vars);

  if (barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    set_new_point(x, basis_index, 1); // value factors needed for coeff grads

    size_t num_act_v = barycentric_active_variables(basis_index);
    if (num_act_v == 0) { // convert 1-D exact indices into n-D colloc index
      size_t pt_index = barycentric_exact_index(basis_index);
      if (pt_index == _NPOS) // if exactIndex but not exactDeltaIndex, then
	tpGradient = 0.;     // at least one dimension has value factor = 0.
      else if (colloc_index.empty())
	copy_data(exp_t1_coeff_grads[pt_index],(int)num_deriv_vars, tpGradient);
      else
	copy_data(exp_t1_coeff_grads[colloc_index[pt_index]],
		  (int)num_deriv_vars, tpGradient);
    }
    else if (num_act_v == numVars) { // interpolation over all variables
      RealMatrix accumulator(numVars, num_deriv_vars); // init to 0.
      precompute_max_keys(basis_index); // if needed for efficiency
      unsigned short key_i0, key_ij, bi_0 = basis_index[0], bi_j,
	max0 = tensor_product_max_key(0, bi_0);
      const RealVector& bc_vf_0
	= polynomialBasis[bi_0][0].barycentric_value_factors();
      size_t i, j, k, num_colloc_pts = key.size();
      Real *accum_0 = accumulator[0], *accum_j, *accum_jm1, bc_vf_00, bc_vf_jj;
      for (i=0; i<num_colloc_pts; ++i) {
	const UShortArray& key_i = key[i]; key_i0 = key_i[0];
	const Real* grad = (colloc_index.empty()) ? exp_t1_coeff_grads[i] :
	  exp_t1_coeff_grads[colloc_index[i]];
	bc_vf_00 = bc_vf_0[key_i0];
	for (j=0; j<num_deriv_vars; ++j)
	  accum_0[j] += grad[j] * bc_vf_00;
	if (key_i0 == max0) {
	  // accumulate sums over variables with max key value
	  for (j=1; j<numVars; ++j) {
	    key_ij = key_i[j]; bi_j = basis_index[j];
	    accum_j = accumulator[j]; accum_jm1 = accumulator[j-1];
	    bc_vf_jj =
	      polynomialBasis[bi_j][j].barycentric_value_factor(key_ij);
	    for (k=0; k<num_deriv_vars; ++k) {
	      accum_j[k]  += accum_jm1[k] * bc_vf_jj;
	      accum_jm1[k] = 0.;
	    }
	    if (key_ij != tensor_product_max_key(j, bi_j))
	      break;
	  }
	}
      }
      Real scale = 1. / barycentric_value_scaling(basis_index);
      Real* accum = accumulator[numVars-1];
      for (j=0; j<num_deriv_vars; ++j)
	tpGradient[j] = accum[j] * scale;
    }
    else { // partial interpolation over active variables
      SizetList pt_factors, act_v_set; //UShortList num_keys;
      size_t num_act_pts, pt_index;
      barycentric_partial_indexing(basis_index, /*num_keys,*/ pt_factors,
				   act_v_set, num_act_pts, pt_index);
      if (pt_index == _NPOS) { tpGradient = 0.; return tpGradient; }
      RealMatrix accumulator(num_deriv_vars, num_act_v); // init to 0.
      size_t i, j, k, v0 = act_v_set.front(), vj;
      unsigned short key_i0, key_ij, bi_0 = basis_index[v0], bi_j,
	max0 = tensor_product_max_key(v0, bi_0);
      const RealVector& bc_vf_v0
	= polynomialBasis[bi_0][v0].barycentric_value_factors();
      size_t pf0 = pt_factors.front(), pfj, prev_pt_set,
	pts_v0_pf0 = pf0 * tensor_product_num_key(v0, bi_0);//*num_keys.front();
      Real *accum_0, *accum_j, *accum_jm1, bc_vf_00, bc_vf_jj;
      SizetList::iterator pf_it, av_it; //UShortList::iterator nk_it;
      // loop over active pts, summing contributions from active variables
      for (i=0; i<num_act_pts; ++i) {
	const UShortArray& key_i = key[pt_index]; key_i0 = key_i[v0];
	const Real* grad = (colloc_index.empty()) ? exp_t1_coeff_grads[pt_index]
	  : exp_t1_coeff_grads[colloc_index[pt_index]];
	accum_0 = accumulator[0]; bc_vf_00 = bc_vf_v0[key_i0];
	for (j=0; j<num_deriv_vars; ++j)
	  accum_0[j] += grad[j] * bc_vf_00;
	pt_index += pf0;
	if (key_i0 == max0) {
	  // accumulate sums over variables with max key value
	  for (j=1, prev_pt_set=pts_v0_pf0, //nk_it=++num_keys.begin(),
	       pf_it=++pt_factors.begin(), av_it=++act_v_set.begin();
	       j<num_act_v; ++j, /*++nk_it,*/ ++pf_it, ++av_it) {
	    vj = *av_it; pfj = *pf_it;
	    key_ij = key_i[vj]; bi_j = basis_index[vj];
	    // update accumulators: push [j-1] entry up to [j] level
	    accum_j = accumulator[j]; accum_jm1 = accumulator[j-1];
	    bc_vf_jj =
	      polynomialBasis[bi_j][vj].barycentric_value_factor(key_ij);
	    for (k=0; k<num_deriv_vars; ++k) {
	      accum_j[k]  += accum_jm1[k] * bc_vf_jj;
	      accum_jm1[k] = 0.;
	    }
	    // update pt_index: prev index rolls back to 0 and curr index +1.
	    // index increment is zero unless active vars are nonconsecutive.
	    if (pfj != prev_pt_set)
	      pt_index += pfj - prev_pt_set;
	    if (key_ij == tensor_product_max_key(vj, bi_j))
	      prev_pt_set = tensor_product_num_key(vj, bi_j) * pfj;//*nk_it*pfj;
	    else
	      break;
	  }
	}
      }
      Real scale = 1. / barycentric_value_scaling(basis_index);
      Real* accum = accumulator[num_act_v-1];
      for (j=0; j<num_deriv_vars; ++j)
	tpGradient[j] = accum[j] * scale;
    }

    /*
    tpGradient = 0.;
    size_t i, num_colloc_pts = key.size();
    for (i=0; i<num_colloc_pts; ++i) {
      const Real* t1_coeff_grad_i = (colloc_index.empty()) ?
	exp_t1_coeff_grads[i] : exp_t1_coeff_grads[colloc_index[i]];
      Real bc_fact = barycentric_value_factor(key[i], basis_index);
      for (j=0; j<num_deriv_vars; ++j)
	tpGradient[j] += t1_coeff_grad_i[j] * bc_fact;
    }
    // apply barycentric denominator
    tpGradient.scale(1. / barycentric_value_scaling(basis_index));
    */
  }
  else {
    tpGradient = 0.;
    size_t i, j, num_colloc_pts = key.size(); Real t1_val;
    for (i=0; i<num_colloc_pts; ++i) {
      const Real* t1_coeff_grad_i = (colloc_index.empty()) ?
	exp_t1_coeff_grads[i] : exp_t1_coeff_grads[colloc_index[i]];
      t1_val = type1_interpolant_value(x, key[i], basis_index);
      for (j=0; j<num_deriv_vars; ++j)
	tpGradient[j] += t1_coeff_grad_i[j] * t1_val;
    }
  }
  return tpGradient;
}


void SharedInterpPolyApproxData::
accumulate_barycentric_partial(const RealVector& t1_coeffs,
			       const UShortArray& basis_index,
			       const UShort2DArray& key,
			       const SizetArray& colloc_index,
			       //const UShortList& num_keys,
			       const SizetList& pt_factors,
			       const SizetList& act_v_set,
			       size_t num_act_pts, size_t pt_index,
			       RealVector& accumulator)
{
  size_t num_act_v = accumulator.length();
  size_t i, j, v0 = act_v_set.front(), vj;
  unsigned short key_i0, key_ij, bi_0 = basis_index[v0], bi_j,
    max0 = tensor_product_max_key(v0, bi_0);
  const RealVector& bc_vf_v0
    = polynomialBasis[bi_0][v0].barycentric_value_factors();
  size_t pf0 = pt_factors.front(), pfj, prev_pt_set,
    pts_v0_pf0 = pf0 * tensor_product_num_key(v0, bi_0);// * num_keys.front();
  SizetList::const_iterator pf_it, av_it; //UShortList::const_iterator nk_it;
  // loop over active pts, summing contributions from active variables
  for (i=0; i<num_act_pts; ++i) {
    const UShortArray& key_i = key[pt_index]; key_i0 = key_i[v0];
    accumulator[0] += (colloc_index.empty()) ? 
      t1_coeffs[pt_index]               * bc_vf_v0[key_i0] :
      t1_coeffs[colloc_index[pt_index]] * bc_vf_v0[key_i0];
    pt_index += pf0;
    if (key_i0 == max0) {
      // accumulate sums over variables with max key value
      for (j=1, prev_pt_set=pts_v0_pf0, //nk_it=++num_keys.begin(),
	   pf_it=++pt_factors.begin(), av_it=++act_v_set.begin();
	   j<num_act_v; ++j, /*++nk_it,*/ ++pf_it, ++av_it) {
	vj = *av_it; pfj = *pf_it; key_ij = key_i[vj]; bi_j = basis_index[vj];
	// update accumulators: push [j-1] entry up to [j] level
	accumulator[j]  += accumulator[j-1] *
	  polynomialBasis[bi_j][vj].barycentric_value_factor(key_ij);
	accumulator[j-1] = 0.;
	// update pt_index: prev index rolls back to 0 and curr index +1.
	// index increment is zero unless active vars are nonconsecutive.
	if (pfj != prev_pt_set)
	  pt_index += pfj - prev_pt_set;
	if (key_ij == tensor_product_max_key(vj, bi_j))
	  prev_pt_set = tensor_product_num_key(vj, bi_j) * pfj;//*nk_it * pfj;
	else
	  break;
      }
    }
  }
}


void SharedInterpPolyApproxData::
accumulate_barycentric_gradient(unsigned short bi_0, unsigned short key_i0,
				size_t ei_0, Real* accum_0, Real t1_coeff,
				const RealVector& bc_vf_0,
				const RealVector& bc_gf_0)
{
  if (bi_0) {
    // accumulate grad factor for comp 0 and value factor for other comps
    accum_0[0] += t1_coeff * bc_gf_0[key_i0];
    if (ei_0 == _NPOS) {
      Real t1_coeff_bc_vf_00 = t1_coeff * bc_vf_0[key_i0];
      for (size_t j=1; j<numVars; ++j)
	accum_0[j] += t1_coeff_bc_vf_00;
    }
    else if (ei_0 == key_i0) // value factor is 1
      for (size_t j=1; j<numVars; ++j)
	accum_0[j] += t1_coeff;
    //else value factor is 0
  }
  else // grad factor is zero, value factor is omitted
    for (size_t j=1; j<numVars; ++j)
      accum_0[j] += t1_coeff;
}


void SharedInterpPolyApproxData::
accumulate_barycentric_gradient(size_t j, unsigned short bi_j,
				unsigned short key_ij, BasisPolynomial& poly_j,
				RealMatrix& accumulator)
{
  Real *accum_j = accumulator[j], *accum_jm1 = accumulator[j-1];
  if (bi_j) {
    size_t ei_j = poly_j.exact_index();
    accum_j[j] += accum_jm1[j] * poly_j.barycentric_gradient_factor(key_ij);
    if (ei_j == _NPOS) { // bc_vf_jj has general value
      Real bc_vf_jj = poly_j.barycentric_value_factor(key_ij);
      for (size_t k=0; k<numVars; ++k) {
	if (k != j) accum_j[k] += accum_jm1[k] * bc_vf_jj;
	accum_jm1[k] = 0.;
      }
    }
    else if (ei_j == key_ij) { // bc_vf_jj is 1
      for (size_t k=0; k<numVars; ++k) {
	if (k != j) accum_j[k] += accum_jm1[k];
	accum_jm1[k] = 0.;
      }
    }
    else // bc_vf_jj is 0
      for (size_t k=0; k<numVars; ++k)
	accum_jm1[k] = 0.;
  }
  else { // grad factor is zero, value factor is omitted
    for (size_t k=0; k<numVars; ++k) {
      if (k != j) accum_j[k] += accum_jm1[k];
      accum_jm1[k] = 0.;
    }
  }
}

} // namespace Pecos
