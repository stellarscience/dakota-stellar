/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedOrthogPolyApproxData
//- Description:  Implementation code for SharedOrthogPolyApproxData class
//-               
//- Owner:        Mike Eldred

#include "SharedOrthogPolyApproxData.hpp"
#include "pecos_global_defs.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG
//#define DECAY_DEBUG

namespace Pecos {


void SharedOrthogPolyApproxData::allocate_data()
{
  // detect changes since previous construction
  bool update_exp_form = (approxOrder != approxOrderPrev);
  //bool restore_exp_form = (multiIndex.size() != t*_*_terms(approxOrder));

  if (update_exp_form) { //|| restore_exp_form) {
    inflate_scalar(approxOrder, numVars); // promote scalar->vector, if needed
    switch (expConfigOptions.expBasisType) {
    case DEFAULT_BASIS: // should not occur (reassigned in NonDPCE ctor)
    case TOTAL_ORDER_BASIS:
      total_order_multi_index(approxOrder, multiIndex);    break;
    case TENSOR_PRODUCT_BASIS:
      tensor_product_multi_index(approxOrder, multiIndex); break;
    }
    precompute_maximal_rules(approxOrder);
    allocate_component_sobol(multiIndex);
    // Note: defer this if update_exp_form is needed downstream
    approxOrderPrev = approxOrder;
  }

  // output (candidate) expansion form
  PCout << "Orthogonal polynomial approximation order = { ";
  for (size_t i=0; i<numVars; ++i)
    PCout << approxOrder[i] << ' ';
  switch (expConfigOptions.expBasisType) {
  case DEFAULT_BASIS: // should not occur (reassigned in NonDPCE ctor)
  case TOTAL_ORDER_BASIS:
    PCout << "} using total-order expansion of ";         break;
  case TENSOR_PRODUCT_BASIS:
    PCout << "} using tensor-product expansion of ";      break;
  }
  PCout << multiIndex.size() << " terms\n";
}


void SharedOrthogPolyApproxData::allocate_data(const UShort2DArray& mi)
{
  multiIndex = mi;
  allocate_component_sobol(mi);

  // output form of imported expansion
  PCout << "Orthogonal polynomial approximation using imported expansion of "
	<< multiIndex.size() << " terms\n";
}


void SharedOrthogPolyApproxData::
allocate_component_sobol(const UShort2DArray& multi_index)
{
  if (expConfigOptions.vbdFlag) {
    if (expConfigOptions.vbdOrderLimit == 1) // main effects only
      allocate_main_sobol();
    else { // main + interaction effects
      sobolIndexMap.clear();
      multi_index_to_sobol_index_map(multi_index);
      assign_sobol_index_map_values();
      /*
      unsigned short max_order = approxOrder[0];
      size_t v, num_v = sharedDataRep->numVars;
      for (v=1; v<num_v; ++v)
	if (approxOrder[v] > max_order)
	  max_order = approxOrder[v];
      if (max_order >= num_v)	{
	if (sobolIndices.empty())
	  allocate_main_interaction_sobol(num_v); // all n-way interactions
      }
      else {
	bool update_exp_form = (approxOrder != approxOrderPrev);
	if (update_exp_form)
	  allocate_main_interaction_sobol(max_order);
      }
      */
    }
  }
}


void SharedOrthogPolyApproxData::
precompute_maximal_rules(const UShort2DArray& multi_index)
{
  // This version scans the multiIndex for maximal order per variable.
  // Due to the overhead of processing multiIndex, protect this processing
  // according to polynomial basis type (would be better to encapsulate
  // within the BasisPolynomial hierarchy).

  size_t i, num_mi_terms = multi_index.size();
  for (i=0; i<numVars; ++i)
    switch (polynomialBasis[i].basis_type()) {
    case NUM_GEN_ORTHOG: {
      unsigned short max_order = multi_index[0][i];
      for (size_t j=1; j<num_mi_terms; ++j)
	if (multi_index[j][i] > max_order)
	  max_order = multi_index[j][i];
      polynomialBasis[i].precompute_rules(max_order);
      break;
    }
    // default is no-op
    }
}


void SharedOrthogPolyApproxData::
precompute_maximal_rules(const UShortArray& approx_order)
{
  // This version employs an incoming approx order, eliminating the overhead
  // of scanning the multiIndex contents.  Since the over head is reduced,
  // can call for each basis polynomial and rely on virtual precompute_rules()
  // to target polynomials that support precomputation optimizations.

  for (size_t i=0; i<numVars; ++i)
    //switch (polynomialBasis[i].basis_type()) {
    //case NUM_GEN_ORTHOG:
    polynomialBasis[i].precompute_rules(approx_order[i]); //break;
    // default is no-op
    //}
}


void SharedOrthogPolyApproxData::
increment_trial_set(CombinedSparseGridDriver* csg_driver,
		    UShort2DArray& aggregated_mi)
{
  size_t last_index = tpMultiIndex.size();
  // increment tpMultiIndex{,Map,MapRef} arrays
  UShort2DArray new_us2a; SizetArray new_sa;
  tpMultiIndex.push_back(new_us2a);
  tpMultiIndexMap.push_back(new_sa); tpMultiIndexMapRef.push_back(0);
  // update tpMultiIndex
  UShortArray exp_order(numVars);
  sparse_grid_level_to_expansion_order(csg_driver, csg_driver->trial_set(),
				       exp_order);
  tensor_product_multi_index(exp_order, tpMultiIndex[last_index]);
  // update multiIndex and append bookkeeping
  append_multi_index(tpMultiIndex[last_index], aggregated_mi,
		     tpMultiIndexMap[last_index],
		     tpMultiIndexMapRef[last_index]);
}


void SharedOrthogPolyApproxData::
decrement_trial_set(const UShortArray& trial_set,
		    UShort2DArray& aggregated_mi, bool save_map)
{
  // reset the aggregated multi-index
  size_t num_exp_terms = tpMultiIndexMapRef.back();
  aggregated_mi.resize(num_exp_terms); // truncate previous increment

  // reset tensor-product bookkeeping and save restorable data
  poppedLevMultiIndex.push_back(trial_set);
  poppedTPMultiIndex.push_back(tpMultiIndex.back());
  if (save_map) { // always needed if we want to mix and match
    poppedTPMultiIndexMap.push_back(tpMultiIndexMap.back());
    poppedTPMultiIndexMapRef.push_back(num_exp_terms);
  }

  tpMultiIndex.pop_back();
  tpMultiIndexMap.pop_back();
  tpMultiIndexMapRef.pop_back();
}


void SharedOrthogPolyApproxData::
pre_push_trial_set(const UShortArray& trial_set,
		      UShort2DArray& aggregated_mi, bool monotonic)
{
  pushIndex = find_index(poppedLevMultiIndex, trial_set);
  size_t last_index = tpMultiIndex.size();

  std::deque<UShort2DArray>::iterator iit = poppedTPMultiIndex.begin();
  std::advance(iit, pushIndex); tpMultiIndex.push_back(*iit);

  // update multiIndex
  if (monotonic) { // reuse previous Map,MapRef bookkeeping if possible
    std::deque<SizetArray>::iterator mit = poppedTPMultiIndexMap.begin();
    std::deque<size_t>::iterator     rit = poppedTPMultiIndexMapRef.begin();
    std::advance(mit, pushIndex);    tpMultiIndexMap.push_back(*mit);
    std::advance(rit, pushIndex);    tpMultiIndexMapRef.push_back(*rit);
    append_multi_index(tpMultiIndex[last_index], tpMultiIndexMap[last_index],
		       tpMultiIndexMapRef[last_index], aggregated_mi);
  }
  else { // replace previous Map,MapRef bookkeeping with new
    SizetArray sa; tpMultiIndexMap.push_back(sa);
    tpMultiIndexMapRef.push_back(0);
    append_multi_index(tpMultiIndex[last_index], aggregated_mi,
		       tpMultiIndexMap[last_index],
		       tpMultiIndexMapRef[last_index]);
  }
}


void SharedOrthogPolyApproxData::
post_push_trial_set(const UShortArray& trial_set,
		       UShort2DArray& aggregated_mi, bool save_map)
{
  std::deque<UShortArray>::iterator   sit = poppedLevMultiIndex.begin();
  std::deque<UShort2DArray>::iterator iit = poppedTPMultiIndex.begin();
  std::advance(sit, pushIndex);       poppedLevMultiIndex.erase(sit);
  std::advance(iit, pushIndex);       poppedTPMultiIndex.erase(iit);
  if (save_map) { // always needed if we want to mix and match
    std::deque<SizetArray>::iterator  mit = poppedTPMultiIndexMap.begin();
    std::deque<size_t>::iterator      rit = poppedTPMultiIndexMapRef.begin();
    std::advance(mit, pushIndex);     poppedTPMultiIndexMap.erase(mit);
    std::advance(rit, pushIndex);     poppedTPMultiIndexMapRef.erase(rit);
  }
}


void SharedOrthogPolyApproxData::
update_component_sobol(const UShort2DArray& multi_index)
{
  if (expConfigOptions.vbdFlag && expConfigOptions.vbdOrderLimit != 1) {
    reset_sobol_index_map_values();
    multi_index_to_sobol_index_map(multi_index);
    assign_sobol_index_map_values();
  }
}


/** Default storage, specialized in derived classes. */
void SharedOrthogPolyApproxData::store_data(size_t index)
{
  // Storing used for multifidelity; popping used for generalized sparse grids

  bool push = (index == _NPOS || index == storedMultiIndex.size());
  if (push) storedMultiIndex.push_back(multiIndex);
  else      storedMultiIndex[index] = multiIndex;
  
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: // both approx order and driver
    if (push) storedApproxOrder.push_back(approxOrder);
    else      storedApproxOrder[index] = approxOrder;
    driverRep->store_grid(index); break;
  case COMBINED_SPARSE_GRID: // driver only
    driverRep->store_grid(index); break;
  default: // approx order only
    if (push) storedApproxOrder.push_back(approxOrder);
    else      storedApproxOrder[index] = approxOrder;
    break;
  }
}


void SharedOrthogPolyApproxData::restore_data(size_t index)
{
  multiIndex = (index == _NPOS)
    ? storedMultiIndex.back() : storedMultiIndex[index];
  
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: // both approx order and driver
    approxOrder = (index == _NPOS)
      ? storedApproxOrder.back() : storedApproxOrder[index];
    driverRep->restore_grid(index); break;
  case COMBINED_SPARSE_GRID: // driver only
    driverRep->restore_grid(index); break;
  default: // approx order only
    approxOrder = (index == _NPOS)
      ? storedApproxOrder.back() : storedApproxOrder[index];
    break;
  }
}


void SharedOrthogPolyApproxData::remove_stored_data(size_t index)
{
  bool pop = (index == _NPOS || index == storedMultiIndex.size() - 1);
  if (pop) storedMultiIndex.pop_back();
  else {
    UShort3DArray::iterator it = storedMultiIndex.begin();
    std::advance(it, index); storedMultiIndex.erase(it);
  }
  
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: // both approx order and driver
    if (pop) storedApproxOrder.pop_back();
    else {
      UShort2DArray::iterator it = storedApproxOrder.begin();
      std::advance(it, index); storedApproxOrder.erase(it);
    }
    driverRep->remove_stored_grid(index); break;
  case COMBINED_SPARSE_GRID: // driver only
    driverRep->remove_stored_grid(index); break;
  default: // approx order only
    if (pop) storedApproxOrder.pop_back();
    else {
      UShort2DArray::iterator it = storedApproxOrder.begin();
      std::advance(it, index); storedApproxOrder.erase(it);
    }
    break;
  }
}


size_t SharedOrthogPolyApproxData::maximal_expansion()
{
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case COMBINED_SPARSE_GRID:
    return driverRep->maximal_grid(); break;
  //case : Not supported: different expansionSamples with same exp order
  default: {
    if (storedApproxOrder.empty()) return _NPOS; // active is maximal
    size_t i, j, max_index = _NPOS, num_stored = storedApproxOrder.size(),
      len = approxOrder.size();
    UShortArray max_ao = approxOrder;
    for (i=0; i<num_stored; ++i) {
      // first test for strict =, < , or >.  If inconclusive, resort to
      // computing total_order_terms().
      bool strict_eq = true, strict_less_eq = true, strict_great_eq = true;
      UShortArray& stored_ao = storedApproxOrder[i];
      for (j=0; j<len; ++j) {
	if      (stored_ao[j] < max_ao[j]) strict_eq = strict_great_eq = false;
	else if (stored_ao[j] > max_ao[j]) strict_eq = strict_less_eq  = false;
      }
      if (strict_eq || strict_less_eq) { }
      else if (strict_great_eq ||
	       total_order_terms(stored_ao) > total_order_terms(max_ao))
	{ max_index = i; max_ao = stored_ao; }
    }
    return max_index; break;
  }
  }
}


void SharedOrthogPolyApproxData::swap_data(size_t index)
{
  std::swap(storedMultiIndex[index], multiIndex);
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    std::swap(storedApproxOrder[index], approxOrder);
    driverRep->swap_grid(index);                             break;
  case COMBINED_SPARSE_GRID: driverRep->swap_grid(index);    break;
  default: std::swap(storedApproxOrder[index], approxOrder); break;
  }
}


size_t SharedOrthogPolyApproxData::pre_combine_data(short combine_type)
{
  // based on incoming combine_type, combine the data stored previously
  // by store_coefficients()

  // Sufficient for two grids: if not currently the maximal grid, then swap
  // with the stored grid (only one is stored)
  //bool swap = !maximal_expansion();
  //if (swap) swap_data();
  
  // For open-ended number of stored grids: retrieve the most refined from the
  // existing grids (from sequence specification + any subsequent refinement)
  size_t max_index = maximal_expansion();
  if (max_index != _NPOS) swap_data(max_index);

  // Most general: overlay all grid refinement levels to create a new superset
  //size_t new_index = overlay_maximal_grid();
  //if (current_grid_index() != new_index) swap_data(new_index);

  switch (combine_type) {
  case ADD_COMBINE: {
    // update multiIndex with any storedMultiIndex terms not yet included.
    // An update in place is sufficient.
    size_t i, stored_mi_map_ref, num_stored = storedMultiIndex.size();
    storedMultiIndexMap.resize(num_stored);
    for (i=0; i<num_stored; ++i)
      //append_multi_index(multiIndex, storedMultiIndex[i], combinedMultiIndex,
      //                   storedMultiIndexMap[i], stored_mi_map_ref);
      append_multi_index(storedMultiIndex[i], multiIndex,
			 storedMultiIndexMap[i], stored_mi_map_ref);
    // reset sobolIndexMap from aggregated multiIndex
    allocate_component_sobol(multiIndex);
    break;
  }
  case MULT_COMBINE: {
    // default implementation: product of total-order expansions
    // (specialized in SharedProjectOrthogPolyApproxData::pre_combine_data())

    // update approxOrder and define combinedMultiIndex
    size_t i, j, num_stored = storedApproxOrder.size();
    for (i=0; i<num_stored; ++i)
      for (j=0; j<numVars; ++j)
	approxOrder[j] += storedApproxOrder[i][j];
    total_order_multi_index(approxOrder, combinedMultiIndex);
    // define sobolIndexMap from combinedMultiIndex
    allocate_component_sobol(combinedMultiIndex);
    break;
  }
  case ADD_MULT_COMBINE:
    PCerr << "Error : additive+multiplicative combination not yet implemented "
	  << "in SharedOrthogPolyApproxData::combine_data()" << std::endl;
    abort_handler(-1);
    break;
  }

  return max_index;
}


void SharedOrthogPolyApproxData::post_combine_data(short combine_type)
{
  // storedMultiIndex and storedApproxOrder used downstream in
  // ProjectOrthogPolyApproximation::integrate_response_moments(),
  // which calls ProjectOrthogPolyApproximation::stored_value()

  //storedMultiIndex.clear(); // needed in {OPA,POPA,ROPA}::stored_value()
  storedMultiIndexMap.clear();
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID:
    driverRep->clear_stored(); break;
  case QUADRATURE:
    //storedApproxOrder.clear(); // needed in ProjectOPA::stored_value()
    driverRep->clear_stored(); break;
  default: // total-order expansions
    storedApproxOrder.clear(); break;
  }

  switch (combine_type) {
  case MULT_COMBINE:
    std::swap(multiIndex, combinedMultiIndex); // pointer swap for efficiency
    combinedMultiIndex.clear();
    break;
  }
}


/** The optional growth_rate supports the option of forcing the computed
    integrand order to be conservative in the presence of exponential growth
    due to nested quadrature rules.  This avoids aggressive formulation of PCE
    expansion orders when an exponential rule takes a large jump that is not
    balanced by the other index set component mappings.  Note that restricting
    the expansion growth directly from the level (*_RESTRICTED_GROWTH cases
    below used by SPARSE_INT_RESTR_TENSOR_SUM_EXP) is similar but not identical
    to restricting the quadrature order growth from the level and then computing
    the integrand and expansion orders from the restricted quadrature order
    (default UNRESTRICTED_GROWTH case below used by SPARSE_INT_TENSOR_SUM_EXP
    and TENSOR_INT_TENSOR_SUM_EXP, where quadrature rule restriction happens
    elsewhere).  In particular, these approaches differ in granularity of
    control, since the former approach grows linearly and the latter approach
    selects the minimal quadrature order (from nonlinear growth or lookup) that
    meets a linear target. */
void SharedOrthogPolyApproxData::
sparse_grid_level_to_expansion_order(CombinedSparseGridDriver* csg_driver,
				     const UShortArray& level,
				     UShortArray& exp_order)
                                     //,short growth_rate)
{
  size_t n = level.size();
  UShortArray int_order(n);
  //switch (growth_rate) {
  //case UNRESTRICTED_GROWTH: { // used for {SPARSE,TENSOR}_INT_TENSOR_SUM_EXP
    // Best option for TENSOR_INT_TENSOR_SUM_EXP, but SPARSE_INT_TENSOR_SUM_EXP
    // is generally too aggressive for nested rules and exponential growth
    // (SPARSE_INT_RESTR_TENSOR_SUM_EXP is preferred).
    UShortArray quad_order(n);
    csg_driver->level_to_order(level, quad_order);
    quadrature_order_to_integrand_order(csg_driver, quad_order, int_order);
    //break;
  //}
  /*
  case SLOW_RESTRICTED_GROWTH: // not currently used
    for (size_t i=0; i<n; ++i) // synch with slow linear growth: i = 2l + 1
      int_order[i] =  2*level[i] + 1;
    break;
  case MODERATE_RESTRICTED_GROWTH: // used for SPARSE_INT_RESTR_TENSOR_SUM_EXP
    // mitigate uneven integrand coverage due to exponential rule growth by
    // enforcing moderate linear expansion growth.
    for (size_t i=0; i<n; ++i) // synch with moderate linear growth: i = 4l + 1
      int_order[i] =  4*level[i] + 1;
    break;
  }
  */
  integrand_order_to_expansion_order(int_order, exp_order);
}


void SharedOrthogPolyApproxData::
quadrature_order_to_integrand_order(IntegrationDriver* int_driver,
				    const UShortArray& quad_order,
				    UShortArray& int_order)
{
  // Need to know exact polynomial resolution for each mixed tensor grid:
  //   Gaussian integrand resolution:        2m-1
  //   Gauss-Patterson integrand resolution: 2m-1 - previous constraints + 1
  //   Clenshaw-Curtis integrand resolution: m (odd m), m-1 (even m)

  // Burkardt monomial test logic:
  //   sparse_grid_monomial_test: resolve monomials of total degree 2*level + 1
  //   for all rules --> doesn't make sense for exponential growth rules where
  //   order grows faster for Gauss than CC (level_to_order exponential is
  //   2^{w+1}-1 for Gauss and 2^w+1 for CC) --> estimate appears valid for CC
  //   (although it does not define a crisp boundary, since some monomials above
  //   the boundary are resolved) but overly conservative for Gauss (whole
  //   orders above the boundary estimate are resolved).

  size_t i, n = quad_order.size();
  if (int_order.size() != n)
    int_order.resize(n);
  const ShortArray& colloc_rules = int_driver->collocation_rules();
  if (colloc_rules.empty()) // use orthogPolyTypes with default modes
    for (i=0; i<n; ++i)
      switch (orthogPolyTypes[i]) {
      case CHEBYSHEV_ORTHOG: // default mode is Clenshaw-Curtis
	int_order[i] = (quad_order[i] % 2) ? quad_order[i] : quad_order[i] - 1;
	break;
      default: // default mode is standard non-nested Gauss rules
	int_order[i] =  2*quad_order[i] - 1; // i = 2m - 1
	break;
      }
  else {
    const UShortArray& gk_order = int_driver->genz_keister_order();
    const UShortArray& gk_prec  = int_driver->genz_keister_precision();
    for (i=0; i<n; ++i)
      switch (colloc_rules[i]) {
      case CLENSHAW_CURTIS: case FEJER2:
	// i = m (odd m), m-1 (even m).  Note that growth rule enforces odd.
	// TO DO: verify FEJER2 same as CC
	int_order[i] = (quad_order[i] % 2) ? quad_order[i] : quad_order[i] - 1;
	break;
      case GAUSS_PATTERSON: {
	// for o(l)=2^{l+1}-1, o(l-1) = (o(l)-1)/2
	unsigned short prev_o = std::max(1,(quad_order[i] - 1)/2);
	int_order[i] = 2*quad_order[i] - prev_o;
	break;
      }
      case GENZ_KEISTER: {
	// same relationship as Gauss-Patterson, except prev_o does not follow
	// simple pattern and requires lookup
	unsigned short lev = 0, max_lev = 5;
	for (lev=0; lev<=max_lev; ++lev)
	  if (gk_order[lev] == quad_order[i])
	    { int_order[i] = gk_prec[lev]; break; }
	/*
	int lev, o, prev_o = 1, max_lev = 4, i_rule = GENZ_KEISTER,
	  g_rule = FULL_EXPONENTIAL; // map l->o directly without restriction
	for (lev=0; lev<=max_lev; ++lev) {
	  webbur::level_growth_to_order(1, &lev, &i_rule, &g_rule, &o);
	  if (o == quad_order[i])
	    { int_order[i] = 2*quad_order[i] - prev_o; break; }
	  else
	    prev_o = o;
	}
	*/
	if (lev > max_lev) {
	  PCerr << "Error: maximum GENZ_KEISTER level exceeded in ProjectOrthog"
		<< "PolyApproximation::quadrature_order_to_integrand_order()."
		<< std::endl;
	  abort_handler(-1);
	}
	break;
      }
      default: // standard non-nested Gauss rules
	int_order[i] =  2*quad_order[i] - 1; break; // i = 2m - 1
      }
  }
}


void SharedOrthogPolyApproxData::
integrand_order_to_expansion_order(const UShortArray& int_order,
				   UShortArray& exp_order)
{
  // reserve half of the integrand order for the expansion and half for the
  // response function (integrand = 2p)
  size_t i, n = int_order.size();
  if (exp_order.size() != n)
    exp_order.resize(n);
  for (i=0; i<n; ++i)
    exp_order[i] = int_order[i] / 2; // remainder truncated
}


/** Append to combined_mi based on append_mi. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& append_mi, UShort2DArray& combined_mi)
{
  if (combined_mi.empty())
    combined_mi = append_mi;
  else {
    size_t i, num_app_mi = append_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      if (std::find(combined_mi.begin(), combined_mi.end(), search_mi) ==
	  combined_mi.end())
	combined_mi.push_back(search_mi);
    }
  }
}


/** Append to combined_mi based on append_mi. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShortArraySet& append_mi, UShort2DArray& combined_mi)
{
  UShortArraySet::const_iterator cit;
  for (cit=append_mi.begin(); cit!=append_mi.end(); ++cit) {
    const UShortArray& search_mi = *cit;
    if (std::find(combined_mi.begin(), combined_mi.end(), search_mi) ==
	combined_mi.end())
      combined_mi.push_back(search_mi);
  }
}


/** Append append_mi to combined_mi, and update append_mi_map
    (SizetArray) and append_mi_map_ref to facilitate related
    aggregations without repeated searching. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& append_mi, UShort2DArray& combined_mi,
		   SizetArray& append_mi_map, size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.resize(num_app_mi);
  if (combined_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map[i] = i;
  }
  else {
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      size_t index = find_index(combined_mi, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	append_mi_map[i] = combined_mi.size();
	combined_mi.push_back(search_mi);
      }
      else
	append_mi_map[i] = index;
    }
  }
}


/** Append append_mi to combined_mi, and update append_mi_map
    (SizetSet) and append_mi_map_ref to facilitate related
    aggregations without repeated searching. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& append_mi, UShort2DArray& combined_mi,
		   SizetSet& append_mi_map, size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.clear();
  if (combined_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map.insert(i);
  }
  else {
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      size_t index = find_index(combined_mi, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	append_mi_map.insert(combined_mi.size());
	combined_mi.push_back(search_mi);
      }
      else
	append_mi_map.insert(index);
    }
  }
}


/** Append append_mi to combined_mi, and update append_mi_map (SizetSet)
    and append_mi_map_ref to facilitate related aggregations without
    repeated searching.  This case is used when append_mi and
    combined_mi follow a consistent order without gaps. */
void SharedOrthogPolyApproxData::
append_leading_multi_index(const UShort2DArray& append_mi,
			   UShort2DArray& combined_mi,
			   SizetSet& append_mi_map, size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.clear();
  if (combined_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map.insert(i);
  }
  else {
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      append_mi_map.insert(i);
      if (i < append_mi_map_ref) {
	// verify that append_mi is a leading subset with consistent ordering
	if (append_mi[i] != combined_mi[i]) {
	  PCerr << "Error: leading subset assumption violated in SharedOrthog"
		<< "PolyApproxData::append_leading_multi_index()." << std::endl;
	  abort_handler(-1);
	}
      }
      else
	combined_mi.push_back(append_mi[i]);
    }
  }
}


/** Append to combined_mi based on append_mi and previously defined
    append_mi_map and append_mi_map_ref.  If necessary, update
    append_mi_map and append_mi_map_ref. */
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& append_mi, SizetArray& append_mi_map,
		   size_t& append_mi_map_ref,   UShort2DArray& combined_mi)
{
  if (combined_mi.empty())
    combined_mi = append_mi; // assume append_mi_map{,_ref} are up to date
  else {
    size_t i, num_app_mi = append_mi.size(), num_mi = combined_mi.size();
    if (num_mi == append_mi_map_ref) { // current mi corresponds to ref
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref)
	  combined_mi.push_back(append_mi[i]);
    }
    else if (num_mi > append_mi_map_ref) { // mi has grown since ref taken
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref) { // previously appended
	  const UShortArray& search_mi = append_mi[i];
	  // search from reference pt forward
	  UShort2DArray::iterator it, it_start = combined_mi.begin();
	  std::advance(it_start, append_mi_map_ref);
	  it = std::find(it_start, combined_mi.end(), search_mi);
	  if (it == combined_mi.end()) { // still an append: update map, append
	    append_mi_map[i] = combined_mi.size();
	    combined_mi.push_back(append_mi[i]);
	  }
	  else // no longer an append: only update map
	    append_mi_map[i] = append_mi_map_ref + std::distance(it_start, it);
	}
      append_mi_map_ref = num_mi; // reference point now updated
    }
    else { // combined_mi is not allowed to shrink since ref taken
      PCerr << "Error: combined_mi inconsistent with reference size in "
	    << "OrthogPolyApproximation::append_multi_index()." << std::endl;
      abort_handler(-1);
    }
  }
}


/** Append to combined_mi based on append_mi and previously defined
    append_mi_map and append_mi_map_ref.  If necessary, update
    append_mi_map and append_mi_map_ref. */
void SharedOrthogPolyApproxData::
append_multi_index(SizetSet& sparse_indices, const UShort2DArray& append_mi,
		   UShort2DArray& combined_mi, RealVector& exp_coeffs,
		   RealMatrix& exp_coeff_grads)
{
  if (combined_mi.empty())
    combined_mi = append_mi; // sparse indices & exp coeffs are up to date
  else { // merge multi-indices; update sparse_indices and exp coeffs

    bool sparse_append = !sparse_indices.empty(); // empty if over-determined LS
    bool coeff_flag = !exp_coeffs.empty(), grad_flag = !exp_coeff_grads.empty();
    RealVector old_exp_coeffs; RealMatrix old_exp_coeff_grads;
    if (coeff_flag) old_exp_coeffs      = exp_coeffs;
    if (grad_flag)  old_exp_coeff_grads = exp_coeff_grads;
 
    size_t i, combined_index, coeff_index, num_app_mi = append_mi.size(),
      num_coeff = (sparse_append) ? sparse_indices.size() : num_app_mi;
    SizetArray append_mi_map(num_app_mi);
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      combined_index = find_index(combined_mi, search_mi);
      if (combined_index == _NPOS) { // search_mi does not exist in combined_mi
	combined_index = combined_mi.size();
	combined_mi.push_back(search_mi);
      }
      append_mi_map[i] = combined_index;
      if (!sparse_append)
	sparse_indices.insert(combined_index); // becomes resorted
    }

    SizetSet old_sparse_indices; SizetSet::iterator it;
    if (sparse_append) {
      old_sparse_indices = sparse_indices;
      sparse_indices.clear();
      for (it=old_sparse_indices.begin(); it!=old_sparse_indices.end(); ++it)
	sparse_indices.insert(append_mi_map[*it]); // becomes resorted
      it = old_sparse_indices.begin(); // reset for loop to follow
    }

    // now that resorting is completed, reorder exp_coeff{s,_grads} to match
    for (i=0; i<num_coeff; ++i) {
      if (sparse_append) { combined_index = append_mi_map[*it]; ++it; }
      else                 combined_index = append_mi_map[i];
      coeff_index = std::distance(sparse_indices.begin(),
				  sparse_indices.find(combined_index));
      if (coeff_flag) exp_coeffs[coeff_index] = old_exp_coeffs[i];
      if (grad_flag) {
	Real *exp_coeff_grad     = exp_coeff_grads[coeff_index],
	     *old_exp_coeff_grad = old_exp_coeff_grads[i];
	for (size_t j=0; j<numVars; ++j)
	  exp_coeff_grad[j] = old_exp_coeff_grad[j];
      }
    }
  }
}


// The following variants maintain a separation between ref_mi and combined_mi,
// rather than updating in place.  In current use cases, append_mi_map has
// provided sufficient bookkeeping to allow in-place updating.

/*  Create combined_mi by appending append_mi to ref_mi.
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& ref_mi, const UShort2DArray& append_mi,
		   UShort2DArray& combined_mi)
{
  if (ref_mi.empty())
    combined_mi = append_mi;
  else {
    combined_mi = ref_mi;
    size_t i, num_app_mi = append_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = append_mi[i];
      if (std::find(combined_mi.begin(), combined_mi.end(),
		    search_mi) == combined_mi.end())
	combined_mi.push_back(search_mi);
    }
  }
}
*/

/*  Append append_mi to ref_mi to create combined_mi, and update
    append_mi_map and append_mi_map_ref to facilitate related
    aggregations without repeated searching.
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& ref_mi, const UShort2DArray& append_mi,
		   UShort2DArray& combined_mi,  SizetArray& append_mi_map,
		   size_t& append_mi_map_ref)
{
  size_t i, num_app_mi = append_mi.size();
  append_mi_map.resize(num_app_mi);
  if (ref_mi.empty()) {
    combined_mi = append_mi;
    append_mi_map_ref = 0;
    for (i=0; i<num_app_mi; ++i)
      append_mi_map[i] = i;
  }
  else {
    combined_mi = ref_mi;
    append_mi_map_ref = combined_mi.size();
    for (i=0; i<num_app_mi; ++i) {
      const UShortArray& search_mi = app_mi[i];
      size_t index = find_index(combined_mi, search_mi);
      if (index == _NPOS) { // search_mi does not yet exist in multi_index
	append_mi_map[i] = combined_mi.size();
	combined_mi.push_back(search_mi);
      }
      else
	append_mi_map[i] = index;
    }
  }
}
*/

/*  Append append_mi to ref_mi to create combined_mi using previously
    defined append_mi_map and append_mi_map_ref.  If necessary, update
    append_mi_map and append_mi_map_ref.
void SharedOrthogPolyApproxData::
append_multi_index(const UShort2DArray& ref_mi, const UShort2DArray& append_mi,
		   SizetArray& append_mi_map, size_t& append_mi_map_ref,
		   UShort2DArray& combined_mi)
{
  if (ref_mi.empty())
    combined_mi = append_mi; // assume append_mi_map{,_ref} are up to date
  else {
    combined_mi = ref_mi;
    size_t i, num_app_mi = append_mi.size(), num_mi = combined_mi.size();
    if (num_mi == append_mi_map_ref) { // current mi corresponds to ref
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref)
	  combined_mi.push_back(append_mi[i]);
    }
    else if (num_mi > append_mi_map_ref) { // mi has grown since ref taken
      for (i=0; i<num_app_mi; ++i)
	if (append_mi_map[i] >= append_mi_map_ref) { // previously appended
	  const UShortArray& search_mi = append_mi[i];
	  // search from reference pt forward
	  UShort2DArray::iterator it, it_start = combined_mi.begin();
	  std::advance(it_start, append_mi_map_ref);
	  it = std::find(it_start, combined_mi.end(), search_mi);
	  if (it == combined_mi.end()) { // still an append: update map, append
	    append_mi_map[i] = combined_mi.size();
	    combined_mi.push_back(append_mi[i]);
	  }
	  else // no longer an append: only update map
	    append_mi_map[i] = append_mi_map_ref + std::distance(it_start, it);
	}
      append_mi_map_ref = num_mi; // reference point now updated
    }
    else { // combined_mi is not allowed to shrink since ref taken
      PCerr << "Error: ref_mi inconsistent with reference size in "
	    << "OrthogPolyApproximation::append_multi_index()." << std::endl;
      abort_handler(-1);
    }
  }
}
*/


void SharedOrthogPolyApproxData::
update_pareto_set(const UShortArray& mi_i, UShort2DArray& combined_pareto)
{
  std::list<UShort2DArray::iterator> removes;
  UShort2DArray::iterator jit;
  bool i_dominated = false, j_dominated;
  for (jit=combined_pareto.begin(); jit!=combined_pareto.end(); ++jit) {
    assess_dominance(mi_i, *jit, i_dominated, j_dominated);
    if (i_dominated) break;
    if (j_dominated) removes.push_back(jit);
  }
  // 
  // prune newly dominated in reverse order (vector iterators
  // following a point of insertion or deletion are invalidated)
  while (!removes.empty())
    { combined_pareto.erase(removes.back()); removes.pop_back(); }
  // add nondominated
  if (!i_dominated)
    combined_pareto.push_back(mi_i);
}


void SharedOrthogPolyApproxData::
update_frontier(const UShortArray& mi_i, UShortArraySet& mi_frontier)
{
  std::list<UShortArraySet::iterator> removes;
  UShortArraySet::iterator jit;
  bool i_dominated = false, j_dominated;
  for (jit=mi_frontier.begin(); jit!=mi_frontier.end(); ++jit) {
    assess_strong_dominance(mi_i, *jit, i_dominated, j_dominated);
    if (i_dominated) break;
    if (j_dominated) removes.push_back(jit);
  }
  // prune newly dominated in any order
  std::list<UShortArraySet::iterator>::iterator rm_iter;
  for (rm_iter=removes.begin(); rm_iter!=removes.end(); ++rm_iter)
    mi_frontier.erase(*rm_iter);
  // add nondominated
  if (!i_dominated)
    mi_frontier.insert(mi_i);
}


/*
bool SharedOrthogPolyApproxData::
assess_dominance(const UShort2DArray& pareto,
		 const UShort2DArray& combined_pareto)
{
  bool new_dominated = true, i_dominated, j_dominated;
  size_t i, j, num_p = pareto.size(),
    num_combined_p = combined_pareto.size();
  for (i=0; i<num_p; ++i) {
    const UShortArray& pareto_i = pareto[i];
    i_dominated = false;
    for (j=0; j<num_combined_p; ++j) {
      assess_dominance(pareto_i, combined_pareto[j], i_dominated, j_dominated);
      if (i_dominated) break;
    }
    if (!i_dominated) {
      new_dominated = false;
#ifdef DEBUG
      PCout << "Nondominated new pareto member =\n" << pareto_i;
#else
      break;
#endif // DEBUG
    }
  }
  return new_dominated;
}
*/


/** Weak Pareto dominance: multi_index a weakly dominates multi_index b
    iff a_i >= b_i for all i and a_i > b_i for at least one dimension.
    Here we add the additional distinction of a challenger versus an 
    incumbent: tie goes to the incumbent (the challenger is dominated
    and is not added redundantly to the Pareto set). */
void SharedOrthogPolyApproxData::
assess_dominance(const UShortArray& new_order,
		 const UShortArray& existing_order,
		 bool& new_dominated, bool& existing_dominated)
{
  // can't use std::vector::operator< (used for component-wise sorting)
  size_t i, n = new_order.size();
  bool equal = true, existing_dominated_temp = true;
  new_dominated = true;
  for (i=0; i<n; ++i)
    if (new_order[i] > existing_order[i])
      { equal = false; new_dominated = false; }
    else if (existing_order[i] > new_order[i])
      { equal = false; existing_dominated_temp = false; }
  // asymmetric logic since incumbent wins a tie
  existing_dominated = (!equal && existing_dominated_temp);
}


/** Strong Pareto dominance: multi_index a strongly dominates
    multi_index b iff a_i > b_i for all i.  This case needs no notion
    of challenger versus incumbent. */
void SharedOrthogPolyApproxData::
assess_strong_dominance(const UShortArray& order_a,
			const UShortArray& order_b,
			bool& a_dominated, bool& b_dominated)
{
  // can't use std::vector::operator< (used for component-wise sorting)
  size_t i, n = order_a.size();
  a_dominated = b_dominated = true;
  for (i=0; i<n; ++i)
    if (order_a[i] == order_b[i])
      { a_dominated = b_dominated = false; break; }
    else if (order_a[i] > order_b[i])
      a_dominated = false;
    else // order_b[i] > order_a[i]
      b_dominated = false;
}


/** This test works in combination with DEBUG settings in
    (Legendre,Laguerre,Jacobi,GenLaguerre)OrthogPolynomial::type1_gradient(). */
void SharedOrthogPolyApproxData::gradient_check()
{
  BasisPolynomial hermite_poly(HERMITE_ORTHOG), legendre_poly(LEGENDRE_ORTHOG),
    laguerre_poly(LAGUERRE_ORTHOG), jacobi_poly(JACOBI_ORTHOG),
    gen_laguerre_poly(GEN_LAGUERRE_ORTHOG), chebyshev_poly(CHEBYSHEV_ORTHOG);
  // alpha/beta selections mirror dakota_uq_rosenbrock_pce.in
  jacobi_poly.alpha_stat(1.5);
  jacobi_poly.beta_stat(2.);
  gen_laguerre_poly.alpha_stat(2.5);

  Real x = 0.5; // valid for all support ranges: [-1,1], [0,Inf], [-Inf, Inf]
  PCout << "-------------------------------------------------\n";
  for (size_t n=0; n<=10; n++) {
    PCout << "Gradients at " << x << " for order " << n << '\n';
    hermite_poly.type1_gradient(x, n);
    legendre_poly.type1_gradient(x, n);
    laguerre_poly.type1_gradient(x, n);
    jacobi_poly.type1_gradient(x, n);
    gen_laguerre_poly.type1_gradient(x, n);
    chebyshev_poly.type1_gradient(x, n);
    PCout << "-------------------------------------------------\n";
  }
}


void SharedOrthogPolyApproxData::
get_tag(char* tag, size_t j, unsigned short order) const
{
  switch (orthogPolyTypes[j]) {
  case HERMITE_ORTHOG:
    std::sprintf(tag,  "He%i", order); break;
  case LEGENDRE_ORTHOG:
    std::sprintf(tag,   "P%i", order); break;
  case LAGUERRE_ORTHOG:
    std::sprintf(tag,   "L%i", order); break;
  case JACOBI_ORTHOG:
    std::sprintf(tag, "Pab%i", order); break;
  case GEN_LAGUERRE_ORTHOG:
    std::sprintf(tag,  "La%i", order); break;
  case CHEBYSHEV_ORTHOG:
    std::sprintf(tag,   "T%i", order); break;
  case NUM_GEN_ORTHOG:
    std::sprintf(tag, "Num%i", order); break;
  default:
    PCerr << "Error: bad polynomial type = " << orthogPolyTypes[j]
	  << " in SharedOrthogPolyApproxData::get_tag()." << std::endl;
    abort_handler(-1);
    break; 
  }
}

} // namespace Pecos
