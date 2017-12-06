/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedPolyApproxData
//- Description:  Implementation code for SharedPolyApproxData class
//-               
//- Owner:        Mike Eldred

#include "SharedPolyApproxData.hpp"
#include "BasisPolynomial.hpp"
#include "DistributionParams.hpp"
#include "NumericGenOrthogPolynomial.hpp"
#include "HahnOrthogPolynomial.hpp"
#include "pecos_stat_util.hpp"

//#define DEBUG


namespace Pecos {

/** This version supports only orthogonal polynomial types.  In this case,
    the polynomial types needed for an orthogonal basis and for computing
    collocation points and weights in an integration driver are the same. */
bool SharedPolyApproxData::
initialize_orthogonal_basis_types_rules(const ShortArray& u_types,
					const BasisConfigOptions& bc_options,
					ShortArray& basis_types,
					ShortArray& colloc_rules)
{
  bool extra_dist_params = false;
  size_t i, num_vars = u_types.size();
  basis_types.resize(num_vars); colloc_rules.resize(num_vars);

  // Initialize basis_types, colloc_rules, and extra_dist_params from u_types:
  for (i=0; i<num_vars; ++i)
    if (initialize_orthogonal_basis_type_rule(u_types[i], bc_options,
					      basis_types[i], colloc_rules[i]))
      extra_dist_params = true;

  return extra_dist_params;
}


bool SharedPolyApproxData::
initialize_orthogonal_basis_type_rule(short u_type,
				      const BasisConfigOptions& bc_options,
				      short& basis_type, short& colloc_rule)
{
  // > open Gauss rules are used for all global cases
  // > nested closed rules are used for piecewise approximations
  bool extra_dist_params;
  switch (u_type) {
  case STD_NORMAL:
    basis_type  = HERMITE_ORTHOG;
    colloc_rule = (bc_options.nestedRules) ? GENZ_KEISTER : GAUSS_HERMITE;
    extra_dist_params = false; break;
  case STD_UNIFORM:
    basis_type = LEGENDRE_ORTHOG;
    // no piecewise orthogonal polynomial bases at this time
    //if (bc_options.piecewiseBasis) // closed nested rules required
    //  colloc_rule = (bc_options.equidistantRules) ?
    //    NEWTON_COTES : CLENSHAW_CURTIS;
    //else
    colloc_rule = (bc_options.nestedRules) ? GAUSS_PATTERSON : GAUSS_LEGENDRE;
    extra_dist_params = false; break;
  case STD_EXPONENTIAL:
    basis_type = LAGUERRE_ORTHOG; colloc_rule = GAUSS_LAGUERRE;
    extra_dist_params = false; break;
  case STD_BETA:
    basis_type = JACOBI_ORTHOG;   colloc_rule = GAUSS_JACOBI;
    extra_dist_params = true; break;
    // a special case of Jacobi is Chebyshev polynomials for alpha=beta=-1/2,
    // which are orthogonal to the weight fn 1/sqrt(1-x^2).
    //case BETA_CHEBY:
    //  basis_type = CHEBYSHEV_ORTHOG; colloc_rule = CLENSHAW_CURTIS; break;
    // Note 1: CLENSHAW_CURTIS/FEJER2 are not currently used for orthogonal
    // projection in Pecos, although they *ARE* used for interpolation (see
    // SharedInterpPolyApproxData::initialize_interpolation_rules()).
    // Note 2: CLENSHAW_CURTIS/FEJER2 are not Gauss rules for Chebyshev
    // polynomials, as these points correspond to extrema rather than roots
    // (point sets from Fejer's first rule are roots, but are not nested).
  case STD_GAMMA:
    basis_type = GEN_LAGUERRE_ORTHOG; colloc_rule = GEN_GAUSS_LAGUERRE;
    extra_dist_params = true; break;
  case POISSON:
    basis_type = CHARLIER_DISCRETE; colloc_rule = GAUSS_CHARLIER;
    extra_dist_params = true; break;
  case BINOMIAL:
    basis_type = KRAWTCHOUK_DISCRETE; colloc_rule = GAUSS_KRAWTCHOUK;
    extra_dist_params = true; break;
  case NEGATIVE_BINOMIAL:
  case GEOMETRIC: // special case of NEGATIVE_BINOMIAL
    basis_type = MEIXNER_DISCRETE; colloc_rule = GAUSS_MEIXNER;
    extra_dist_params = true; break;
  case HYPERGEOMETRIC:
    basis_type = HAHN_DISCRETE; colloc_rule = GAUSS_HAHN;
    extra_dist_params = true; break;
  default:
    basis_type = NUM_GEN_ORTHOG; colloc_rule = GOLUB_WELSCH;
    extra_dist_params = true; break;
  }

  return extra_dist_params;
}


void SharedPolyApproxData::
initialize_polynomial_basis(const ShortArray& basis_types,
			    const ShortArray& colloc_rules,
			    std::vector<BasisPolynomial>& poly_basis)
{
  size_t i, num_vars = basis_types.size(), num_rules = colloc_rules.size();

  // instantiate poly_basis using basis_types and colloc_rules
  if (poly_basis.size() != num_vars) {
    poly_basis.resize(num_vars);
    if (num_rules == num_vars)
      for (i=0; i<num_vars; ++i)
	poly_basis[i] = BasisPolynomial(basis_types[i], colloc_rules[i]);
    else if (num_rules == 0)   // use default rules for each basis type
      for (i=0; i<num_vars; ++i)
	poly_basis[i] = BasisPolynomial(basis_types[i]);
    else if (num_rules == 1) { // cubature utilizes a single rule
      short colloc_rule = colloc_rules[0];
      for (i=0; i<num_vars; ++i)
	poly_basis[i] = BasisPolynomial(basis_types[i], colloc_rule);
    }

    /*
    // Could reuse objects as in InterpPolyApproximation, but this would require
    // caching dist params and moving code from distribution_parameters() to
    // here, which isn't yet well motivated.
    size_t i, j;
    for (i=0; i<numVars; ++i) {
      // reuse prev instance via shared rep or instantiate new unique instance
      short basis_type_i = basis_types[i];
      bool found = false;
      for (j=0; j<i; ++j)
	if ( basis_type_i == basis_types[j] && ( basis_type_i <= LAGUERRE ||
	     ( basis_type_i == JACOBI && jacobiAlphas[] == jacobiAlphas[] &&
	       jacobiBetas[] == jacobiBetas[] ) ||
	     ( basis_type_i == GENERALIZED_LAGUERRE &&
	       genLagAlphas[] == genLagAlphas[] ) ) )
	  { found = true; break; }
      polynomialBasis[i]
	= (found) ? polynomialBasis[j] : BasisPolynomial(basis_type_i);
    }
    */
  }
}


void SharedPolyApproxData::
update_basis_distribution_parameters(const ShortArray& u_types,
				     const AleatoryDistParams& adp,
				     std::vector<BasisPolynomial>& poly_basis)
{
  size_t i, num_vars = u_types.size(), nuv_cntr = 0, lnuv_cntr = 0,
    luuv_cntr = 0, tuv_cntr = 0, beuv_cntr = 0, gauv_cntr = 0, guuv_cntr = 0,
    fuv_cntr = 0, wuv_cntr = 0, hbuv_cntr = 0, biuv_cntr=0, nbuv_cntr=0,
    puv_cntr=0, geuv_cntr=0, hguv_cntr=0, 
    hpiuv_cntr=0, hpsuv_cntr=0, hpruv_cntr=0;

  // update poly_basis using distribution data from dp
  for (i=0; i<num_vars; ++i)
    switch (u_types[i]) {
    case STD_NORMAL:
      ++nuv_cntr; break;
    case STD_UNIFORM: case STD_EXPONENTIAL:
      break;
    case STD_BETA:
      poly_basis[i].alpha_stat(adp.beta_alpha(beuv_cntr));
      poly_basis[i].beta_stat(adp.beta_beta(beuv_cntr));
      ++beuv_cntr; break;
    case STD_GAMMA:
      poly_basis[i].alpha_stat(adp.gamma_alpha(gauv_cntr));
      ++gauv_cntr; break;
    case BOUNDED_NORMAL:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	bounded_normal_distribution(adp.normal_mean(nuv_cntr),
				    adp.normal_std_deviation(nuv_cntr),
				    adp.normal_lower_bound(nuv_cntr),
				    adp.normal_upper_bound(nuv_cntr));
      ++nuv_cntr; break;
    case LOGNORMAL: {
      Real lambda, zeta;
      params_from_lognormal_spec(adp.lognormal_means(),
				 adp.lognormal_std_deviations(),
				 adp.lognormal_lambdas(), adp.lognormal_zetas(),
				 adp.lognormal_error_factors(), lnuv_cntr,
				 lambda, zeta);
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	lognormal_distribution(lambda, zeta);
      ++lnuv_cntr; break;
    }
    case BOUNDED_LOGNORMAL: {
      Real lambda, zeta;
      params_from_lognormal_spec(adp.lognormal_means(),
				 adp.lognormal_std_deviations(),
				 adp.lognormal_lambdas(), adp.lognormal_zetas(),
				 adp.lognormal_error_factors(), lnuv_cntr,
				 lambda, zeta);
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	bounded_lognormal_distribution(lambda, zeta,
				       adp.lognormal_lower_bound(lnuv_cntr),
				       adp.lognormal_upper_bound(lnuv_cntr));
      ++lnuv_cntr; break;
    }
    case LOGUNIFORM:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	loguniform_distribution(adp.loguniform_lower_bound(luuv_cntr),
				adp.loguniform_upper_bound(luuv_cntr));
      ++luuv_cntr; break;
    case TRIANGULAR:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	triangular_distribution(adp.triangular_lower_bound(tuv_cntr),
				adp.triangular_mode(tuv_cntr),
				adp.triangular_upper_bound(tuv_cntr));
      ++tuv_cntr; break;
    case GUMBEL:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	gumbel_distribution(adp.gumbel_alpha(guuv_cntr),
			    adp.gumbel_beta(guuv_cntr));
      ++guuv_cntr; break;
    case FRECHET:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	frechet_distribution(adp.frechet_alpha(fuv_cntr),
			     adp.frechet_beta(fuv_cntr));
      ++fuv_cntr; break;
    case WEIBULL:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	weibull_distribution(adp.weibull_alpha(wuv_cntr),
			     adp.weibull_beta(wuv_cntr));
      ++wuv_cntr; break;
    case HISTOGRAM_BIN:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
	histogram_bin_distribution(adp.histogram_bin_pairs(hbuv_cntr));
      ++hbuv_cntr; break;
    case BINOMIAL:
      poly_basis[i].alpha_stat(adp.binomial_probability_per_trial(biuv_cntr));
      poly_basis[i].beta_stat(adp.binomial_num_trials(biuv_cntr));
      ++biuv_cntr; break;
    case NEGATIVE_BINOMIAL:
      poly_basis[i].alpha_stat(
	   adp.negative_binomial_probability_per_trial(nbuv_cntr));
      poly_basis[i].beta_stat(adp.negative_binomial_num_trials(nbuv_cntr));
      ++nbuv_cntr; break;
    case POISSON:
      poly_basis[i].alpha_stat(adp.poisson_lambda(puv_cntr));
      ++puv_cntr; break;
    case GEOMETRIC:
      // Use Meixner polynomial for Geometric variables setting num_trials to 1
      poly_basis[i].alpha_stat(
	   adp.geometric_probability_per_trial(geuv_cntr));
      poly_basis[i].beta_stat(1);
      ++geuv_cntr; break;
    case HYPERGEOMETRIC:
      poly_basis[i].alpha_stat(adp.hypergeometric_total_population(hguv_cntr));
      poly_basis[i].beta_stat(adp.hypergeometric_selected_population(hguv_cntr));
      ((HahnOrthogPolynomial*)poly_basis[i].polynomial_rep())->gamma_stat(adp.hypergeometric_num_drawn(hguv_cntr));
      ++hguv_cntr; break;
    case HISTOGRAM_PT_INT:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
        histogram_pt_distribution(adp.histogram_point_int_pairs(hpiuv_cntr));
      ++hpiuv_cntr; break;
    case HISTOGRAM_PT_STRING:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
        histogram_pt_distribution(adp.histogram_point_string_pairs(hpsuv_cntr));
      ++hpsuv_cntr; break;
    case HISTOGRAM_PT_REAL:
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())->
        histogram_pt_distribution(adp.histogram_point_real_pairs(hpruv_cntr));
      ++hpruv_cntr; break;
    default:
      PCerr << "Error: unsupported u-space type in SharedPolyApproxData::"
	    << "distribution_parameters()" << std::endl;
      abort_handler(-1);
      break;
    }
}


void SharedPolyApproxData::increment_data()
{
  PCerr << "Error: derived class does not redefine increment_data()."
	<< std::endl;
  abort_handler(-1);
}


void SharedPolyApproxData::decrement_data()
{
  PCerr << "Error: derived class does not redefine decrement_data()."
	<< std::endl;
  abort_handler(-1);
}


void SharedPolyApproxData::pre_push_data()
{ } // default implementation is no op


void SharedPolyApproxData::post_push_data()
{ } // default implementation is no op


void SharedPolyApproxData::pre_finalize_data()
{ } // default implementation is no op


void SharedPolyApproxData::post_finalize_data()
{ } // default implementation is no op


size_t SharedPolyApproxData::pre_combine_data(short combine_type)
{ return _NPOS; } // default implementation is no op, no swap


void SharedPolyApproxData::post_combine_data(short combine_type)
{ } // default implementation is no op


/*
void SharedPolyApproxData::allocate_component_sobol()
{
  // default implementation is reasonable for tensor expansions, but is
  // wasteful (and should be overridden) for total-order and sparse 
  // sum-of-tensor (which emulate total-order) expansions.
  if (expConfigOptions.vbdFlag && expConfigOptions.expansionCoeffFlag &&
      sobolIndices.empty()) {
    switch (expConfigOptions.vbdOrderLimit) {
    case 1: // main effects only
      allocate_main_sobol();                    break;
    default: // main + interactions
      allocate_main_interaction_sobol(numVars); break;
    }
  }
}

void SharedPolyApproxData::
allocate_main_interaction_sobol(unsigned short max_order)
{
  // include m-way interactions for m <= max_order.  Note: sobol index 0 is
  // unused (corresponds to constant exp term with no variable dependence).
  unsigned long sobol_len = 0; size_t v;
  BitArray set(numVars, 0);
  if (max_order >= 1) {
    sobol_len += numVars;
    for (v=0; v<numVars; ++v)
      { set.set(v); sobolIndexMap[set] = v; set.reset(v); }
  }
  unsigned long cntr = sobol_len;
  for (unsigned short ord=2; ord<=max_order; ++ord) {
    // compute number of terms from n_choose_k() for each exp order
    sobol_len += BasisPolynomial::n_choose_k(numVars, ord);
    // use increment_terms() and exclude term equalities to get interactions
    UShortArray terms(ord, 1); // # of terms = level
    bool order_complete = false;
    while (!order_complete) {
      size_t last_index = ord - 1, prev_index = ord - 2;
      for (terms[last_index]=1; terms[last_index]<terms[prev_index];
	   ++terms[last_index]) {
	// convert orders (within terms) to variable indices (within set)
	set.reset();
	for (size_t i=0; i<ord; ++i)
	  set.set(terms[i]-1);
	sobolIndexMap[set] = cntr; ++cntr;
      }
      increment_terms(terms, last_index, prev_index, numVars,
		      order_complete, true);
    }
  }

  sobolIndices.sizeUninitialized(sobol_len);
}
*/


void SharedPolyApproxData::allocate_main_sobol()
{
  if (sobolIndexMap.empty()) {
    // define binary sets corresponding to main effects
    BitArray set(numVars, 0);
    // prepend the 0-way interaction for indexing consistency
    sobolIndexMap[set] = 0;
    for (size_t v=0; v<numVars; ++v) // activate bit for variable v
      { set.set(v); sobolIndexMap[set] = v+1; set.reset(v); }
  }
}


void SharedPolyApproxData::
multi_index_to_sobol_index_map(const UShort2DArray& mi)
{
  // component Sobol' indices are limited by two factors:
  // (1) the terms appearing in the multi-index
  // (2) an additional user specification (interactionOrderLimit)

  BitArray set(numVars);
  size_t i, j, interactions, num_mi = mi.size();
  for (i=0; i<num_mi; ++i) {
    // determine the bit set corresponding to this expansion term
    interactions = 0;
    for (j=0; j<numVars; ++j)
      if (mi[i][j]) { set.set(j); ++interactions; } // activate bit j
      else          set.reset(j);                 // deactivate bit j

    // If set does not already exist, insert it.
    // For value in key-value pair, initially use the interaction order;
    // will be updated below in assign_sobol_index_map_values().  The
    // 0-way interaction is included to support child lookup requirements
    // in InterpPolyApproximation::compute_partial_variance().
    if ( ( !expConfigOptions.vbdOrderLimit ||                // no limit
	   interactions <= expConfigOptions.vbdOrderLimit )  // within limit
	 && sobolIndexMap.find(set) == sobolIndexMap.end() ) // new set
      sobolIndexMap[set] = interactions; // order of interaction
  }
}


void SharedPolyApproxData::reset_sobol_index_map_values()
{
  // return sobolIndexMap values to the order of interaction
  // (the starting point for assign_sobol_index_map_values())
  for (BAULMIter it=sobolIndexMap.begin(); it!=sobolIndexMap.end(); ++it)
    it->second = it->first.count(); // order of interaction
}


void SharedPolyApproxData::assign_sobol_index_map_values()
{
  // Compute total counts for interactions within each interaction group.
  // The 0-way interaction is included to support child lookup requirements
  // in InterpPolyApproximation::compute_partial_variance().
  ULongArray counters(numVars+1, 0);
  for (BAULMIter it=sobolIndexMap.begin(); it!=sobolIndexMap.end(); ++it)
    ++counters[it->second];

  // aggregate the totals into starting indices for each interaction group
  ULongArray indices(numVars+1); indices[0] = 0;
  for (size_t i=1; i<=numVars; ++i)
    indices[i] = indices[i-1] + counters[i-1];
  //unsigned long sobol_len = indices[numVars] + counters[numVars];

  // Now that the sobolIndexMap sets are fully defined and sorted, define the
  // mapping of these sets to sobolIndices based on counters.  We choose to
  // group the sobolIndices by the number of interactions to simplify output.
  unsigned long interaction_order;
  for (BAULMIter it=sobolIndexMap.begin(); it!=sobolIndexMap.end(); ++it) {
    interaction_order = it->second;
    it->second = indices[interaction_order]++;
  }
}


/** Return the number of terms in a tensor-product expansion.  For
    isotropic and anisotropic expansion orders, calculation of the
    number of expansion terms is straightforward: Prod(p_i + 1). */
size_t SharedPolyApproxData::
tensor_product_terms(const UShortArray& order, bool include_upper_bound)
{
  size_t i, n = order.size(), num_terms = 1;
  if (include_upper_bound)
    for (i=0; i<n; ++i)
      num_terms *= order[i] + 1; // multi-index from expansion order p (default)
  else
    for (i=0; i<n; ++i)
      num_terms *= order[i];     // collocation index from quadrature order m
  return num_terms;
}


void SharedPolyApproxData::
tensor_product_multi_index(const UShortArray& order, UShort2DArray& multi_index,
			   bool include_upper_bound)
{
  // This function may be used for defining collocation indices for quadrature
  // points (indices = 0:m-1) or polynomial order multi-indices for tensor
  // expansions (indices = 0:p).  The inclusion or exclusion of the upper bound
  // (m or p) is managed with include_upper_bound (false for m, true for p).

  // rather than inserting into multi_index, go ahead and estimate its length
  // (since its inexpensive) and use term-by-term assignment.
  size_t i, n = order.size(),
    mi_len = tensor_product_terms(order, include_upper_bound);
  if (mi_len != multi_index.size())
    multi_index.resize(mi_len);
  UShortArray indices(n, 0); multi_index[0] = indices;
  for (i=1; i<mi_len; ++i) {
    // increment the n-dimensional index set
    increment_indices(indices, order, include_upper_bound);
    multi_index[i] = indices;
  }
}


void SharedPolyApproxData::
hierarchical_tensor_product_multi_index(const UShort2DArray& delta_quad,
					UShort2DArray& multi_index)
{
  // delta_quad are non-sequential indices of hierarchical collocation point
  // increments
  size_t i, j, n = delta_quad.size(), mi_len = 1;
  UShortArray index_bound(n);
  for (i=0; i<n; ++i)
    mi_len *= index_bound[i] = delta_quad[i].size();
  if (mi_len != multi_index.size())
    multi_index.resize(mi_len);
  UShortArray indices(n, 0);
  for (i=0; i<mi_len; ++i) {
    multi_index[i].resize(n);
    for (j=0; j<n; ++j)
      multi_index[i][j] = delta_quad[j][indices[j]];
    if (i < mi_len-1)
      increment_indices(indices, index_bound, false); // 0 <= index < bound
  }
}


/** Return the number of terms in a total-order expansion.  For anisotropic
    expansion order, no simple expression is currently available and the
    number of expansion terms is computed using the multiIndex recursion. */
size_t SharedPolyApproxData::
total_order_terms(const UShortArray& upper_bound, short lower_bound_offset)
{
  size_t i, n = upper_bound.size();
  if (!n) {
    PCerr << "Error: empty upper_bound in SharedPolyApproxData::"
	  << "total_order_terms()." << std::endl;
    abort_handler(-1);
  }

  bool isotropic = true;
  unsigned short order = upper_bound[0];
  for (i=1; i<n; ++i)
    if (upper_bound[i] != order)
      { isotropic = false; break; }

  size_t num_terms;
  if (isotropic) {
    num_terms = (size_t)std::floor(
      BasisPolynomial::n_choose_k(order+n, order)+.5); // round non-integral
    if (lower_bound_offset >= 0) { // default is -1
      int omit_order = order - lower_bound_offset - 1;
      if (omit_order >= 0)
	num_terms -= (size_t)std::floor(
	  BasisPolynomial::n_choose_k(omit_order+n, omit_order)+.5); // round
    }
  }
  else { // anisotropic: use multiIndex recursion to compute
    bool mi_lower_bound = (lower_bound_offset >= 0); // default is -1
    if (mi_lower_bound) {
      // Smolyak combinatorial form is invalid for anisotropic levels
      PCerr << "Error: anisotropic orders not currently supported with "
	    << "multi-index lower bound\n       in SharedPolyApproxData::"
	    << "total_order_terms()." << std::endl;
      abort_handler(-1);
    }
    unsigned short max_order = order, order_nd;
    for (i=1; i<n; ++i)
      max_order = std::max(max_order, upper_bound[i]);
    num_terms = 1; // order 0
    if (max_order >= 1)   // order 1
      for (i=0; i<n; ++i)
	if (upper_bound[i] >= 1) // only upper bound may be anisotropic
	  ++num_terms;
    for (order_nd=2; order_nd<=max_order; ++order_nd) { // order 2 through max
      UShortArray terms(order_nd, 1); // # of terms = current order
      bool order_complete = false;
      while (!order_complete) {
	size_t last_index = order_nd - 1, prev_index = order_nd - 2;
	for (terms[last_index]=1; terms[last_index]<=terms[prev_index]; 
	     ++terms[last_index]) {
	  bool include = true;
	  for (i=0; i<n; ++i) {
	    if (std::count(terms.begin(), terms.end(), i+1) > upper_bound[i]) {
	      // only upper bound may be anisotropic
	      include = false;
	      break;
	    }
	  }
	  if (include)
	    ++num_terms;
	}
	// increment term hierarchy
	increment_terms(terms, last_index, prev_index, n, order_complete);
      }
    }
  }
  return num_terms;
}


/** Overloaded version for defining the multi-indices for a single
    scalar level.  Anisotropy is not supported, so this version is not
    usable as a kernel within other overloaded versions. */
void SharedPolyApproxData::
total_order_multi_index(unsigned short level, size_t num_vars, 
			UShort2DArray& multi_index)
{
  UShortArray mi(num_vars, 0);
  multi_index.clear();
  // special logic required for level < 2 due to prev_index defn below
  switch (level) {
  case 0:
    multi_index.push_back(mi); break;
  case 1:
    for (size_t i=0; i<num_vars; ++i)
      { mi[i] = 1; multi_index.push_back(mi); mi[i] = 0; }
    break;
  default: {
    UShortArray terms(level, 1); // # of terms = level
    bool order_complete = false;
    while (!order_complete) {
      // this is the inner-most loop w/i the nested looping managed by terms
      size_t last_index = level - 1, prev_index = level - 2;
      for (terms[last_index]=1;
	   terms[last_index]<=terms[prev_index]; ++terms[last_index]) {
	for (size_t i=0; i<num_vars; ++i)
	  mi[i] = std::count(terms.begin(), terms.end(), i+1);
	multi_index.push_back(mi);
      }
      increment_terms(terms, last_index, prev_index, num_vars, order_complete);
    }
    break;
  }
  }
}


void SharedPolyApproxData::
total_order_multi_index(const UShortArray& upper_bound,
			UShort2DArray& multi_index, short lower_bound_offset, 
			size_t max_terms)
{
  // populate multi_index: implementation follows ordering of Eq. 4.1 in
  // [Xiu and Karniadakis, 2002].
  // To handle anisotropy, we currently perform a complete total-order
  // recursion based on max_order and reject multi-indices with components
  // greater than those in upper_bound.
  size_t i, cntr = 0, n = upper_bound.size();
  if (!n) {
    PCerr << "Error: empty upper_bound in SharedPolyApproxData::"
	  << "total_order_multi_index()." << std::endl;
    abort_handler(-1);
  }

  // Since total_order_terms() is expensive, use insertions into multi_index
  // rather than sizing with entry-by-entry assignment
  bool isotropic = true;
  unsigned short order = upper_bound[0];
  for (i=1; i<n; ++i)
    if (upper_bound[i] != order)
      { isotropic = false; break; }

  bool mi_lower_bound = (lower_bound_offset >= 0); // default is -1
  unsigned short max_order = order, min_order = 0, order_nd;
  if (isotropic) {
    if (mi_lower_bound)
      min_order = (lower_bound_offset >= order)
	? 0 : order - lower_bound_offset; // e.g., Smolyak l.b. = w-N+1
  }
  else {
    if (mi_lower_bound) {
      // Smolyak combinatorial form is invalid for anisotropic levels
      PCerr << "Error: anisotropic orders not currently supported with "
	    << "multi-index lower bound\n       in SharedPolyApproxData::"
	    << "total_order_multi_index()." << std::endl;
      abort_handler(-1);
    }
    for (i=1; i<n; ++i)
      max_order = std::max(max_order, upper_bound[i]);
  }

  // special logic required for order_nd < 2 due to prev_index defn below
  UShortArray mi(n,0);
  multi_index.clear();
  if (min_order == 0) // && max_order >= 0
    multi_index.push_back(mi); // order 0
  if (min_order <= 1 && max_order >= 1) { // order 1
    for (i=0; i<n && cntr<max_terms; ++i) {
      if (upper_bound[i] >= 1) { // only upper bound may be anisotropic
	mi[i] = 1; // ith entry is nonzero
	multi_index.push_back(mi);
	mi[i] = 0; // reset
      }
    }
  }
  for (order_nd=std::max(min_order,(unsigned short)2);
       order_nd<=max_order; ++order_nd) {
    UShortArray terms(order_nd, 1);//# of terms = current order
    bool order_complete = false;
    while (!order_complete) {
      // this is the inner-most loop w/i the nested looping managed by terms
      size_t last_index = order_nd - 1, prev_index = order_nd - 2;
      for (terms[last_index]=1;
	   terms[last_index]<=terms[prev_index] && cntr<max_terms;
	   ++terms[last_index]) {
	// store the orders of the univariate polynomials to be used for
	// constructing the current multivariate basis function
	bool include = true;
	for (i=0; i<n; ++i) {
	  mi[i] = std::count(terms.begin(), terms.end(), i+1);
	  if (mi[i] > upper_bound[i]) { // only upper bound may be anisotropic
	    include = false;
	    break;
	  }
	}
	if (include)
	  multi_index.push_back(mi);
#ifdef DEBUG
	PCout << "terms:\n" << terms << std::endl;
#endif // DEBUG
      }
      if (cntr == max_terms)
	order_complete = true;
      else // increment term hierarchy
	increment_terms(terms, last_index, prev_index, n, order_complete);
    }
  }

#ifdef DEBUG
  size_t mi_len = multi_index.size();
  for (i=0; i<mi_len; ++i)
    PCout << "multiIndex[" << i << "]:\n" << multi_index[i] << std::endl;
#endif // DEBUG
}

} // namespace Pecos
