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
#include "MultivariateDistribution.hpp"
#include "pecos_global_defs.hpp"
#include "pecos_math_util.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG
//#define DECAY_DEBUG

namespace Pecos {


void SharedOrthogPolyApproxData::allocate_data()
{
  UShortArray&  approx_order =  approxOrdIter->second;
  UShort2DArray& multi_index = multiIndexIter->second;

  // detect changes since previous construction
  // *** TO DO: replace with a flag updated in set functions, once sobol
  // ***        index bookkeeping is more modular.
  bool update_exp_form
    = (approx_order != prevApproxOrder || activeKey != prevActiveKey);
  //bool restore_exp_form = (multiIndex.size() != t*_*_terms(approxOrder));

  if (/*updateExpForm*/update_exp_form) { //|| restore_exp_form) {
    inflate_scalar(approx_order, numVars); // promote scalar->vector, if needed
    switch (expConfigOptions.expBasisType) {
    case DEFAULT_BASIS: // should not occur (reassigned in NonDPCE ctor)
    case TOTAL_ORDER_BASIS:
      total_order_multi_index(approx_order, multi_index);    break;
    case TENSOR_PRODUCT_BASIS:
      tensor_product_multi_index(approx_order, multi_index); break;
    }
    precompute_maximal_rules(approx_order);
    allocate_component_sobol(multi_index);
    // Note: defer this if update_exp_form is needed downstream
    prevApproxOrder = approx_order;
    prevActiveKey   = activeKey;
  }

  // output (candidate) expansion form
  PCout << "Orthogonal polynomial approximation order = { ";
  for (size_t i=0; i<numVars; ++i)
    PCout << approx_order[i] << ' ';
  switch (expConfigOptions.expBasisType) {
  case DEFAULT_BASIS: // should not occur (reassigned in NonDPCE ctor)
  case TOTAL_ORDER_BASIS:
    PCout << "} using total-order expansion of ";         break;
  case TENSOR_PRODUCT_BASIS:
    PCout << "} using tensor-product expansion of ";      break;
  }
  PCout << multi_index.size() << " terms\n";
}


void SharedOrthogPolyApproxData::active_key(const UShortArray& key)
{
  if (activeKey != key) {
    activeKey = key; // SharedPolyApproxData::active_key(key);
    update_active_iterators();

    switch (expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE: case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID:
      driverRep->active_key(key); break;
    }
  }
}


void SharedOrthogPolyApproxData::clear_keys()
{
  SharedPolyApproxData::clear_keys();

  approxOrder.clear(); approxOrdIter = approxOrder.end();
  multiIndex.clear(); multiIndexIter = multiIndex.end();
  tpMultiIndex.clear(); tpMultiIndexMap.clear(); tpMultiIndexMapRef.clear();
  poppedMultiIndex.clear(); poppedMultiIndexMap.clear();
  poppedMultiIndexMapRef.clear();

  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID:
    driverRep->clear_keys(); break;
  }
}


void SharedOrthogPolyApproxData::allocate_data(const UShort2DArray& multi_index)
{
  multiIndexIter->second = multi_index;
  allocate_component_sobol(multi_index);

  // output form of imported expansion
  PCout << "Orthogonal polynomial approximation using imported expansion of "
	<< multi_index.size() << " terms\n";
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
	bool update_exp_form
	  = (approxOrder != prevApproxOrder || activeKey != prevActiveKey);
	if (update_exp_form)
	  allocate_main_interaction_sobol(max_order);
      }
      */
    }
  }
}


bool SharedOrthogPolyApproxData::push_available()
{
  switch (expConfigOptions.refineControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
    IncrementalSparseGridDriver* isg_driver
      = (IncrementalSparseGridDriver*)driverRep;
    return isg_driver->push_trial_available();
    break;
  }
  //case UNIFORM_CONTROL:  case DIMENSION_ADAPTIVE_CONTROL_SOBOL:
  //case DIMENSION_ADAPTIVE_CONTROL_DECAY:
  default:
    return !poppedMultiIndex[activeKey].empty(); break;
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
  UShort3DArray& tp_mi         = tpMultiIndex[activeKey];
  Sizet2DArray&  tp_mi_map     = tpMultiIndexMap[activeKey];
  SizetArray&    tp_mi_map_ref = tpMultiIndexMapRef[activeKey];
  size_t last_index = tp_mi.size();
  // increment tpMultiIndex{,Map,MapRef} arrays
  UShort2DArray new_us2a; SizetArray new_sa;
  tp_mi.push_back(new_us2a);
  tp_mi_map.push_back(new_sa); tp_mi_map_ref.push_back(0);
  // update tpMultiIndex
  UShortArray exp_order(numVars);
  sparse_grid_level_to_expansion_order(csg_driver, csg_driver->trial_set(),
				       exp_order);
  tensor_product_multi_index(exp_order, tp_mi[last_index]);
  // update multiIndex and append bookkeeping
  append_multi_index(tp_mi[last_index], aggregated_mi,
		     tp_mi_map[last_index], tp_mi_map_ref[last_index]);
}


void SharedOrthogPolyApproxData::
decrement_trial_set(const UShortArray& trial_set,
		    UShort2DArray& aggregated_mi, bool save_map)
{
  UShort3DArray& tp_mi         = tpMultiIndex[activeKey];
  Sizet2DArray&  tp_mi_map     = tpMultiIndexMap[activeKey];
  SizetArray&    tp_mi_map_ref = tpMultiIndexMapRef[activeKey];
  // reset the aggregated multi-index
  size_t num_exp_terms = tp_mi_map_ref.back();
  aggregated_mi.resize(num_exp_terms); // truncate previous increment

  // reset tensor-product bookkeeping and save restorable data
  //poppedLevMultiIndex[activeKey].push_back(trial_set);
  poppedMultiIndex[activeKey].push_back(tp_mi.back());
  if (save_map) { // always needed if we want to mix and match
    poppedMultiIndexMap[activeKey].push_back(tp_mi_map.back());
    poppedMultiIndexMapRef[activeKey].push_back(num_exp_terms);
  }

  tp_mi.pop_back();  tp_mi_map.pop_back();  tp_mi_map_ref.pop_back();
}


void SharedOrthogPolyApproxData::
pre_push_trial_set(const UShortArray& trial_set,
		      UShort2DArray& aggregated_mi, bool monotonic)
{
  UShort3DArray& tp_mi         = tpMultiIndex[activeKey];
  Sizet2DArray&  tp_mi_map     = tpMultiIndexMap[activeKey];
  SizetArray&    tp_mi_map_ref = tpMultiIndexMapRef[activeKey];

  // retrieve restoration index, previously-computed by ISGDriver
  size_t p_index = push_index(), last_index = tp_mi.size();

  UShort2DArrayDeque::iterator iit
    = poppedMultiIndex[activeKey].begin() + p_index;
  tp_mi.push_back(*iit);

  // update multiIndex
  if (monotonic) { // reuse previous Map,MapRef bookkeeping if possible
    SizetArrayDeque::iterator mit
      = poppedMultiIndexMap[activeKey].begin() + p_index;
    SizetDeque::iterator rit
      = poppedMultiIndexMapRef[activeKey].begin() + p_index;
    tp_mi_map.push_back(*mit);  tp_mi_map_ref.push_back(*rit);
    append_multi_index(tp_mi[last_index], tp_mi_map[last_index],
		       tp_mi_map_ref[last_index], aggregated_mi);
  }
  else { // replace previous Map,MapRef bookkeeping with new
    SizetArray sa; tp_mi_map.push_back(sa); tp_mi_map_ref.push_back(0);
    append_multi_index(tp_mi[last_index], aggregated_mi,
		       tp_mi_map[last_index], tp_mi_map_ref[last_index]);
  }
}


void SharedOrthogPolyApproxData::
post_push_trial_set(const UShortArray& trial_set,
		       UShort2DArray& aggregated_mi, bool save_map)
{
  //UShortArrayDeque& popped_lev_mi = poppedLevMultiIndex[activeKey];
  //UShortArrayDeque::iterator sit = popped_lev_mi.begin() + p_index;
  //popped_lev_mi.erase(sit);

  // retrieve restoration index, previously-computed by ISGDriver
  size_t p_index = push_index();

  UShort2DArrayDeque& popped_tp_mi = poppedMultiIndex[activeKey];
  UShort2DArrayDeque::iterator iit = popped_tp_mi.begin() + p_index;
  popped_tp_mi.erase(iit);

  if (save_map) { // always needed if we want to mix and match
    SizetArrayDeque& popped_tp_mi_map = poppedMultiIndexMap[activeKey];
    SizetArrayDeque::iterator mit = popped_tp_mi_map.begin() + p_index;
    popped_tp_mi_map.erase(mit);

    SizetDeque& popped_tp_mi_map_ref = poppedMultiIndexMapRef[activeKey];
    SizetDeque::iterator rit = popped_tp_mi_map_ref.begin() + p_index;
    popped_tp_mi_map_ref.erase(rit);
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


const UShortArray& SharedOrthogPolyApproxData::maximal_expansion()
{
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID:
    return driverRep->maximal_grid(); break;
  //case : Not supported: different expansionSamples with same exp order
  default: {
    std::map<UShortArray, UShortArray>::iterator
      ao_it = approxOrder.begin(), max_it = ao_it;
    size_t j, len = ao_it->second.size();
    ++ao_it;
    for (; ao_it!=approxOrder.end(); ++ao_it) {
      // first test for strict =, < , or >; if inconclusive, resort to
      // computing total_order_terms()
      bool strict_eq = true, strict_less_eq = true, strict_great_eq = true;
      UShortArray& ao = ao_it->second; UShortArray& max_ao = max_it->second;
      for (j=0; j<len; ++j) {
	if      (ao[j] < max_ao[j]) strict_eq = strict_great_eq = false;
	else if (ao[j] > max_ao[j]) strict_eq = strict_less_eq  = false;
      }
      if (strict_eq || strict_less_eq) { }
      else if (strict_great_eq ||
	       total_order_terms(ao) > total_order_terms(max_ao))
	max_it = ao_it;
    }
    return max_it->first; // return form/level key
    break;
  }
  }
}


void SharedOrthogPolyApproxData::pre_combine_data()
{
  // For open-ended number of stored grids: retrieve the most refined from the
  // existing grids (from sequence specification + any subsequent refinement)
  //active_key(maximal_expansion());

  switch (expConfigOptions.combineType) {
  case MULT_COMBINE: {
    // default implementation: product of total-order expansions
    // (specialized in SharedProjectOrthogPolyApproxData::pre_combine_data())

    // roll up approxOrders to define combinedMultiIndex
    size_t cntr, j, num_seq = approxOrder.size() - 2; // bridge from 1st to last
    if (num_seq) combinedMultiIndexSeq.resize(num_seq);
    std::map<UShortArray, UShortArray>::iterator ao_it = approxOrder.begin();
    UShortArray combined_ao = ao_it->second;   ++ao_it; // copy
    for (cntr=0; ao_it!=approxOrder.end(); ++ao_it, ++cntr) {
      const UShortArray& ao = ao_it->second;
      for (j=0; j<numVars; ++j)
	combined_ao[j] += ao[j];
      UShort2DArray& combined_mi = (cntr < num_seq) ?
	combinedMultiIndexSeq[cntr] : combinedMultiIndex;
      total_order_multi_index(combined_ao, combined_mi);
    }
    break;
  }
  case ADD_MULT_COMBINE:
    PCerr << "Error : additive+multiplicative combination not yet implemented "
	  << "in SharedOrthogPolyApproxData::combine_data()" << std::endl;
    abort_handler(-1);
    break;
  default: { //case ADD_COMBINE:
    // combine all multiIndex keys into combinedMultiIndex{,Map}
    size_t i, num_combine = multiIndex.size(), combine_mi_map_ref;
    combinedMultiIndex.clear();  combinedMultiIndexMap.resize(num_combine);
    std::map<UShortArray, UShort2DArray>::iterator mi_it;
    for (mi_it=multiIndex.begin(), i=0; mi_it!=multiIndex.end(); ++mi_it, ++i)
      append_multi_index(mi_it->second, combinedMultiIndex,
			 combinedMultiIndexMap[i], combine_mi_map_ref);
    /*
    // update active (maximal) multiIndex with any non-active multiIndex
    // terms not yet included.  An update in place is sufficient.
    size_t i, num_combine = multiIndex.size() - 1, cntr = 0, combine_mi_map_ref;
    combinedMultiIndexMap.resize(num_combine);
    std::map<UShortArray, UShort2DArray>::iterator mi_it;
    combinedMultiIndex = multiIndexIter->second; // copy
    for (mi_it = multiIndex.begin(); mi_it != multiIndex.end(); ++mi_it)
      if (mi_it->first != activeKey) {
	append_multi_index(mi_it->second, combinedMultiIndex,
			   combinedMultiIndexMap[cntr], combine_mi_map_ref);
	++cntr;
      }
    */
    break;
  }
  }

  // reset sobolIndexMap from aggregated multiIndex
  //allocate_component_sobol(combinedMultiIndex);
}


void SharedOrthogPolyApproxData::combined_to_active(bool clear_combined)
{
  // retrieve the most refined from the existing grids (from sequence
  // specification + any subsequent refinement)
  // *** This might still be a good idea (e.g., avoid temporarily inflating
  // *** memory footprint), though not strictly required...
  //active_key(maximal_expansion());

  // combine level data for Driver into combinedSmolyak{MultiIndex,Coeffs},
  // combined{Var,T1Weight,T2Weight}Sets, et al.
  // Note: unlike Nodal SC, PCE only requires grid combination at the end, so
  //       incremental combined grid updates are not necessary for efficiency.
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: case QUADRATURE:
  //case CUBATURE: // key-based data not implemented for multilevel
    driverRep->combine_grid();
    driverRep->combined_to_active(clear_combined);
    break;
  }

  if (clear_combined) {
    std::swap(multiIndexIter->second, combinedMultiIndex); // pointer swap
    combinedMultiIndex.clear();
    combinedMultiIndexMap.clear(); combinedMultiIndexSeq.clear();
  }
  else
    multiIndexIter->second = combinedMultiIndex; // copy

  allocate_component_sobol();//(multiIndexIter->second);
}


void SharedOrthogPolyApproxData::clear_inactive_data()
{
  bool ao = false, tp = false;
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID:
    tp = true; driverRep->clear_inactive(); break;
  case QUADRATURE:
    ao = true; driverRep->clear_inactive(); break;
  default: // total-order expansions
    ao = true;
    if (expConfigOptions.expBasisType == ADAPTED_BASIS_GENERALIZED)
      tp = true;
    break;
  }

  std::map<UShortArray, UShort2DArray>::iterator  mi_it = multiIndex.begin();
  std::map<UShortArray, UShortArray>::iterator    ao_it = approxOrder.begin();
  std::map<UShortArray, UShort3DArray>::iterator tp1_it = tpMultiIndex.begin();
  std::map<UShortArray, Sizet2DArray>::iterator  tp2_it
    = tpMultiIndexMap.begin();
  std::map<UShortArray, SizetArray>::iterator    tp3_it
    = tpMultiIndexMapRef.begin();
  while (mi_it != multiIndex.end())
    if (mi_it == multiIndexIter) { // preserve active
      ++mi_it;
      if (ao) ++ao_it;
      if (tp) { ++tp1_it; ++tp2_it;  ++tp3_it; }
    }
    else { // clear inactive: postfix increments manage iterator invalidations
      multiIndex.erase(mi_it++);
      if (ao) approxOrder.erase(ao_it++);
      if (tp) {
	tpMultiIndex.erase(tp1_it++);
	tpMultiIndexMap.erase(tp2_it++);
	tpMultiIndexMapRef.erase(tp3_it++);
      }
    }
}


void SharedOrthogPolyApproxData::
product_multi_index(const UShort2DArray& multi_index_a,
		    const UShort2DArray& multi_index_b,
		    UShort2DArray& multi_index_c)
{
  // c = a * b

  // Note: pareto set approach works when there is a conversion from sm_mi
  // to tensor-product multi-index.  Working only with the multi-indices,
  // we need to carry lower order (dominated) terms in order to fill in
  // lower-order product terms.

  size_t i, j, v, num_mi_a = multi_index_a.size(),
    num_mi_b = multi_index_b.size();
  UShortArray prod_mi(numVars); UShortArraySet prod_mi_sets;
  for (i=0; i<num_mi_a; ++i)
    for (j=0; j<num_mi_b; ++j) {
      for (v=0; v<numVars; ++v)
	prod_mi[v] = multi_index_a[i][v] + multi_index_b[j][v];
      prod_mi_sets.insert(prod_mi); // sorted, unique
    }

  size_t num_mi_c = prod_mi_sets.size();
  multi_index_c.resize(num_mi_c);
  UShortArraySet::iterator mi_it;
  for (i=0, mi_it=prod_mi_sets.begin(); i<num_mi_c; ++i, ++mi_it)
    multi_index_c[i] = *mi_it;
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


/** This test works in combination with DEBUG settings in
    (Legendre,Laguerre,Jacobi,GenLaguerre)OrthogPolynomial::type1_gradient(). */
void SharedOrthogPolyApproxData::gradient_check()
{
  BasisPolynomial hermite_poly(HERMITE_ORTHOG), legendre_poly(LEGENDRE_ORTHOG),
    laguerre_poly(LAGUERRE_ORTHOG), jacobi_poly(JACOBI_ORTHOG),
    gen_laguerre_poly(GEN_LAGUERRE_ORTHOG), chebyshev_poly(CHEBYSHEV_ORTHOG);
  // alpha/beta selections mirror dakota_uq_rosenbrock_pce.in
  jacobi_poly.push_parameter(BE_ALPHA, 1.5);
  jacobi_poly.push_parameter(BE_BETA, 2.);
  gen_laguerre_poly.push_parameter(GA_ALPHA, 2.5);

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
