/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:       SharedRegressOrthogPolyApproxData
//- Description: Implementation code for SharedRegressOrthogPolyApproxData class
//-               
//- Owner:       John Jakeman

#include "SharedRegressOrthogPolyApproxData.hpp"
#include "pecos_global_defs.hpp"
#include "pecos_math_util.hpp"
#include "SurrogateData.hpp"

namespace Pecos {


void SharedRegressOrthogPolyApproxData::allocate_data()
{
  UShortArray&   approx_order =  approxOrdIter->second;
  UShort2DArray& multi_index  = multiIndexIter->second;

  if (expConfigOptions.expCoeffsSolnApproach == ORTHOG_LEAST_INTERPOLATION) {
    // clear history from previous expansion; new pts -> new least interpolant
    approx_order.clear();  // for update_approx_order() -> exp combination logic
    multi_index.clear();   // for reuse check in ROPA::least_interpolation()
    sobolIndexMap.clear(); // for update_component_sobol()

    if (expConfigOptions.vbdFlag && expConfigOptions.vbdOrderLimit == 1)
      allocate_main_sobol(); // main effects only

    PCout << "Orthogonal polynomial approximation of least order\n";
    return;
  }
  else if (expConfigOptions.expBasisType == DEFAULT_BASIS     ||
	   expConfigOptions.expBasisType == TOTAL_ORDER_BASIS ||
	   expConfigOptions.expBasisType == TENSOR_PRODUCT_BASIS)
    { SharedOrthogPolyApproxData::allocate_data(); return; }

  // Adapted basis selection:

  // detect changes since previous construction
  bool update_exp_form
    = (approx_order != prevApproxOrder || activeKey != prevActiveKey);
  //bool restore_exp_form = (multi_index.size() != t*_*_terms(approx_order));

  if (update_exp_form) { //|| restore_exp_form) {
    switch (expConfigOptions.expBasisType) {
    // ADAPTED_BASIS_GENERALIZED starts from a reference sparse grid basis
    // (defined from initSGLevel), as for generalized sparse grids.  A
    // key difference is the use of multiIndexGrowthFactor, since the actual
    // quadrature growth rules provide no benefit in this context.
    case ADAPTED_BASIS_GENERALIZED: {
      // We could assume that:
      // (1) exp_order defines the upper bound of the basis (not the starting
      //     point) --> e.g., we are in 100D and would like to recover terms up
      //     to order 5, but can't form a candidate multiIndex that large. So we
      //     start from SGL 0 and select components up to this upper bnd.  In 
      //     this case, we specify either colloc_points or colloc_ratio << 1.
      // (2) exp_order defines the starting point and there is no explicit
      //     upper bound --> colloc_ratio is more (initially) meaningful.
      // (3) neither: exp_order only used to define colloc_pts from colloc_ratio
      //     and we hard-wire the starting point (e.g., level 0) for adaptation
      // For now, we use case 3 as it is the simplest

      // initialize the sparse grid driver (lightweight mode) for generating
      // candidate index sets
      lsgDriver.initialize_grid(numVars, regressConfigOptions.initSGLevel);

      // define reference multi_index and tpMultiIndex{,Map,MapRef} from 
      // initial sparse grid level
      //sparse_grid_multi_index(&lsgDriver, multi_index); // heavyweight mapping
      multi_index.clear();
      tpMultiIndex[activeKey].clear();
      tpMultiIndexMap[activeKey].clear();
      tpMultiIndexMapRef[activeKey].clear();
      const UShort2DArray& sm_mi = lsgDriver.smolyak_multi_index();
      size_t i, num_sm_mi = sm_mi.size();
      for (i=0; i<num_sm_mi; ++i)
	increment_trial_set(sm_mi[i], multi_index); // lightweight mapping
      break;
    }
    case ADAPTED_BASIS_EXPANDING_FRONT:
      inflate_scalar(approx_order, numVars);// promote scalar->vector, if needed
      total_order_multi_index(approx_order, multi_index);
      break;
    }
    allocate_component_sobol(multi_index);
    // Note: defer this if update_exp_form is needed downstream
    prevApproxOrder = approx_order;
    prevActiveKey   = activeKey;
  }

  // output (candidate) expansion form
  PCout << "Orthogonal polynomial approximation order = { ";
  for (size_t i=0; i<numVars; ++i) PCout << approx_order[i] << ' ';
  PCout << "} using adapted expansion initiated from " << multi_index.size()
	<< " terms\n";
}


void SharedRegressOrthogPolyApproxData::
update_approx_order(unsigned short new_order)
{
  UShortArray& approx_order = approxOrdIter->second;
  if (approx_order.empty() || new_order > approx_order[0])
    approx_order.assign(numVars, new_order);
}


void SharedRegressOrthogPolyApproxData::approx_order_to_multi_index()
{
  // approxOrder updated from NonDPCE --> propagate to multiIndex
  UShortArray&  approx_order =  approxOrdIter->second;
  UShort2DArray& multi_index = multiIndexIter->second;

  switch (expConfigOptions.expBasisType) {
  case TENSOR_PRODUCT_BASIS:
    tensor_product_multi_index(approx_order, multi_index); break;
  default:
    total_order_multi_index(approx_order, multi_index);    break;
  }
}


void SharedRegressOrthogPolyApproxData::increment_data()
{
  // TO DO: ADAPTED_BASIS_GENERALIZED case ???

  // To automatically update approxOrder, would need to either infer a
  // collocation ratio (based on initial data size + initial approxOrder)
  // or have it set from NonDPCE.  Inferring is problematic due to rounding 
  // effects with discrete sample counts and access to settings for the 
  // general case (useDerivs, termsOrder).  Then increment_order() would be
  // justified by incremented data size.

  // for decrement
  prevMultiIndex = multiIndexIter->second;

  // Better: manage increments from NonDPCE using ratio_samples_to_order()
  //         (e.g., see NonDQUESOBayesCalibration::update_model())
  approx_order_to_multi_index();
  allocate_component_sobol(multiIndexIter->second);
  //prevApproxOrder = approx_order;
  //prevActiveKey   = activeKey;
}


void SharedRegressOrthogPolyApproxData::decrement_data()
{
  // TO DO: ADAPTED_BASIS_GENERALIZED case ???

  poppedMultiIndex[activeKey].push_back(multiIndexIter->second);

  //approx_order_to_multi_index();
  multiIndexIter->second = prevMultiIndex;

  //allocate_component_sobol(multiIndexIter->second);
  //prevApproxOrder = approx_order;
  //prevActiveKey   = activeKey;
}


void SharedRegressOrthogPolyApproxData::pre_push_data()
{
  // TO DO: ADAPTED_BASIS_GENERALIZED case ???

  // currently returns 0 for all cases other than generalized sparse grids
  size_t p_index = push_index();

  // for decrement
  prevMultiIndex = multiIndexIter->second;

  std::map<UShortArray, UShort2DArrayDeque>::iterator pop_it
    = poppedMultiIndex.find(activeKey);
  UShort2DArrayDeque::iterator u2a_it;
  if (pop_it == poppedMultiIndex.end() || pop_it->second.size() <= p_index) {
    PCerr << "Error: lookup failure in SharedRegressOrthogPolyApproxData::"
	  << "pre_push_data()." << std::endl;
    abort_handler(-1);
  }
  else {
    u2a_it = pop_it->second.begin();   std::advance(u2a_it, p_index);
    multiIndexIter->second = *u2a_it;  pop_it->second.erase(u2a_it);
  }
  //approx_order_to_multi_index();
  allocate_component_sobol(multiIndexIter->second);
  //prevApproxOrder = approx_order;
  //prevActiveKey   = activeKey;
}


/** Append append_mi to combined_mi, and update append_mi_map (SizetSet)
    and append_mi_map_ref to facilitate related aggregations without
    repeated searching.  This case is used when append_mi and
    combined_mi follow a consistent order without gaps. */
void SharedRegressOrthogPolyApproxData::
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
	  PCerr << "Error: leading subset assumption violated in SharedRegress"
		<< "OrthogPolyApproxData::append_leading_multi_index()."
		<< std::endl;
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
void SharedRegressOrthogPolyApproxData::
append_sparse_multi_index(SizetSet& sparse_indices,
			  const UShort2DArray& append_mi,
			  UShort2DArray& combined_mi,
			  RealVector& exp_coeffs, RealMatrix& exp_coeff_grads)
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


void SharedRegressOrthogPolyApproxData::
increment_trial_set(const UShortArray& trial_set, UShort2DArray& aggregated_mi)
{
  UShort3DArray& tp_mi         = tpMultiIndex[activeKey];
  Sizet2DArray&  tp_mi_map     = tpMultiIndexMap[activeKey];
  SizetArray&    tp_mi_map_ref = tpMultiIndexMapRef[activeKey];

  size_t i, last_index = tp_mi.size();
  // increment tpMultiIndex{,Map,MapRef} arrays
  UShort2DArray new_us2a; SizetArray new_sa;
  tp_mi.push_back(new_us2a);
  tp_mi_map.push_back(new_sa); tp_mi_map_ref.push_back(0);
  // update tpMultiIndex
  UShortArray exp_order(numVars);
  // linear growth in Gaussian rules would normally result in a factor of 2:
  //   m = 2l+1 -> i = 2m-1 = 4l+1 -> o = i/2 = 2l
  // This is the default, but finer and coarser grain growth can be used.
  unsigned short growth_fact = regressConfigOptions.multiIndexGrowthFactor;
  for (i=0; i<numVars; ++i)
    exp_order[i] = growth_fact * trial_set[i];
  tensor_product_multi_index(exp_order, tp_mi[last_index]);
  // update multiIndex and append bookkeeping
  append_multi_index(tp_mi[last_index], aggregated_mi,
		     tp_mi_map[last_index], tp_mi_map_ref[last_index]);
}


bool SharedRegressOrthogPolyApproxData::
set_restriction(UShort2DArray& aggregated_mi, SizetSet& sparse_indices,
		SizetSet& save_tp)
{
  if (sparse_indices.empty()) // dense multi-index: no restriction possible
    return false;

  UShort3DArray& tp_mi         = tpMultiIndex[activeKey];
  Sizet2DArray&  tp_mi_map     = tpMultiIndexMap[activeKey];
  SizetArray&    tp_mi_map_ref = tpMultiIndexMapRef[activeKey];

  // determine the TP multi-indices to save
  StSCIter cit;
  size_t i, num_tp_mi = tp_mi_map_ref.size(), last_index = num_tp_mi - 1,
    sparse_ind, save_tp_cntr, new_tp_cntr;
  bool map_ref_find;
  for (cit=sparse_indices.begin(); cit!=sparse_indices.end(); ++cit) {
    sparse_ind = *cit; map_ref_find = false;
    for (i=0; i<last_index; ++i)
      // upper bound for i-th multi-index append defined by MapRef[i+1]
      if (sparse_ind < tp_mi_map_ref[i+1])
	{ map_ref_find = true; break; }
    if (map_ref_find) save_tp.insert(i); // occurs within i-th appended TP
    else              save_tp.insert(last_index); // must be last TP
  }

  // prune unneeded TP multi-index sets and update Map, MapRef, sparse_indices
  size_t num_save = save_tp.size();
  if (num_save == num_tp_mi)
    return false;
  else {
    UShort2DArray old_aggregated_mi(aggregated_mi);   aggregated_mi.clear();
    SizetSet      old_sparse_indices(sparse_indices); sparse_indices.clear();
    for (cit=save_tp.begin(), new_tp_cntr=0; cit!=save_tp.end();
	 ++cit, ++new_tp_cntr) {
      save_tp_cntr = *cit;
      if (save_tp_cntr != new_tp_cntr) // reuse previous tpMultiIndex entries
	tp_mi[new_tp_cntr] = tp_mi[save_tp_cntr];
      append_multi_index(tp_mi[new_tp_cntr], aggregated_mi,
			 tp_mi_map[new_tp_cntr], tp_mi_map_ref[new_tp_cntr]);
    }
    // prune old records off the end
    tp_mi.resize(num_save);
    tp_mi_map.resize(num_save); tp_mi_map_ref.resize(num_save);
    // update sparse_indices using old_aggregated_mi
    // TO DO: review SharedPolyApproxData::append_multi_index(SizetSet&)
    //        for more efficient logic?
    for (cit=old_sparse_indices.begin(); cit!=old_sparse_indices.end(); ++cit)
      sparse_indices.insert(find_index(aggregated_mi, old_aggregated_mi[*cit]));

    return true;
  }
}


void SharedRegressOrthogPolyApproxData::
pack_polynomial_data(const RealVector& c_vars, const UShortArray& mi,
		     bool add_val,  double* pack_val,  size_t& pv_cntr,
		     bool add_grad, double* pack_grad, size_t& pg_cntr)
{
  if (add_val)
    { pack_val[pv_cntr] = multivariate_polynomial(c_vars, mi); ++pv_cntr; }
  if (add_grad) {
    const RealVector& mvp_grad
      = multivariate_polynomial_gradient_vector(c_vars, mi);
    for (size_t j=0; j<numVars; ++j, ++pg_cntr)
      pack_grad[pg_cntr] = mvp_grad[j];
  }
}


void SharedRegressOrthogPolyApproxData::
pack_response_data(const SurrogateDataResp& sdr,
		   bool add_val,  double* pack_val,  size_t& pv_cntr,
		   bool add_grad, double* pack_grad, size_t& pg_cntr)
{
  if (add_val)
    { pack_val[pv_cntr] = sdr.response_function(); ++pv_cntr; }
  if (add_grad) {
    const RealVector& resp_grad = sdr.response_gradient();
    for (size_t j=0; j<numVars; ++j, ++pg_cntr)
      pack_grad[pg_cntr] = resp_grad[j];
  }
}

} // namespace Pecos
