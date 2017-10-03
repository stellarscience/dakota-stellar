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
#include "SurrogateData.hpp"


namespace Pecos {

void SharedRegressOrthogPolyApproxData::allocate_data()
{
  if (expConfigOptions.expCoeffsSolnApproach == ORTHOG_LEAST_INTERPOLATION) {
    // clear history from previous expansion; new pts -> new least interpolant
    approxOrder.clear();   // for update_approx_order() -> exp combination logic
    multiIndex.clear();    // for reuse check in ROPA::least_interpolation()
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
  bool update_exp_form = (approxOrder != approxOrderPrev);
  //bool restore_exp_form = (multiIndex.size() != t*_*_terms(approxOrder));

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

      // define reference multiIndex and tpMultiIndex{,Map,MapRef} from 
      // initial sparse grid level
      //sparse_grid_multi_index(&lsgDriver, multiIndex); // heavyweight mapping
      multiIndex.clear();
      tpMultiIndex.clear(); tpMultiIndexMap.clear(); tpMultiIndexMapRef.clear();
      const UShort2DArray& sm_mi = lsgDriver.smolyak_multi_index();
      size_t i, num_sm_mi = sm_mi.size();
      for (i=0; i<num_sm_mi; ++i)
	increment_trial_set(sm_mi[i], multiIndex); // lightweight mapping
      break;
    }
    case ADAPTED_BASIS_EXPANDING_FRONT:
      inflate_scalar(approxOrder, numVars); // promote scalar->vector, if needed
      total_order_multi_index(approxOrder, multiIndex);
      break;
    }
    allocate_component_sobol(multiIndex);
    // Note: defer this if update_exp_form is needed downstream
    approxOrderPrev = approxOrder;
  }

  // output (candidate) expansion form
  PCout << "Orthogonal polynomial approximation order = { ";
  for (size_t i=0; i<numVars; ++i) PCout << approxOrder[i] << ' ';
  PCout << "} using adapted expansion initiated from " << multiIndex.size()
	<< " terms\n";
}


void SharedRegressOrthogPolyApproxData::
update_approx_order(unsigned short new_order)
{
  if (approxOrder.empty() || new_order > approxOrder[0])
    approxOrder.assign(numVars, new_order);
}


void SharedRegressOrthogPolyApproxData::increment_data()
{
  // To automatically update approxOrder, would need to either infer a
  // collocation ratio (based on initial data size + initial approxOrder)
  // or have it set from NonDPCE.  Inferring is problematic due to rounding 
  // effects with discrete sample counts and access to settings for the 
  // general case (useDerivs, termsOrder).  Then increment_order() would be
  // justified by incremented data size.

  // Better: manage increments from NonDPCE using ratio_samples_to_order()
  //         (e.g., see NonDQUESOBayesCalibration::update_model())

  // approxOrder updated from NonDPolynomialChaos --> propagate to multiIndex
  bool update_exp_form = (approxOrder != approxOrderPrev);
  if (update_exp_form) {
    switch (expConfigOptions.expBasisType) {
    case TENSOR_PRODUCT_BASIS:
      tensor_product_multi_index(approxOrder, multiIndex); break;
    default:
      total_order_multi_index(approxOrder, multiIndex);    break;
    }
    allocate_component_sobol(multiIndex);
    approxOrderPrev = approxOrder;
  }
}


void SharedRegressOrthogPolyApproxData::
increment_trial_set(const UShortArray& trial_set, UShort2DArray& aggregated_mi)
{
  size_t i, last_index = tpMultiIndex.size();
  // increment tpMultiIndex{,Map,MapRef} arrays
  UShort2DArray new_us2a; SizetArray new_sa;
  tpMultiIndex.push_back(new_us2a);
  tpMultiIndexMap.push_back(new_sa); tpMultiIndexMapRef.push_back(0);
  // update tpMultiIndex
  UShortArray exp_order(numVars);
  // linear growth in Gaussian rules would normally result in a factor of 2:
  //   m = 2l+1 -> i = 2m-1 = 4l+1 -> o = i/2 = 2l
  // This is the default, but finer and coarser grain growth can be used.
  unsigned short growth_fact = regressConfigOptions.multiIndexGrowthFactor;
  for (i=0; i<numVars; ++i)
    exp_order[i] = growth_fact * trial_set[i];
  tensor_product_multi_index(exp_order, tpMultiIndex[last_index]);
  // update multiIndex and append bookkeeping
  append_multi_index(tpMultiIndex[last_index], aggregated_mi,
		     tpMultiIndexMap[last_index],
		     tpMultiIndexMapRef[last_index]);
}


bool SharedRegressOrthogPolyApproxData::
set_restriction(UShort2DArray& aggregated_mi, SizetSet& sparse_indices,
		SizetSet& save_tp)
{
  if (sparse_indices.empty()) // dense multi-index: no restriction possible
    return false;

  // determine the TP multi-indices to save
  StSCIter cit;
  size_t i, num_tp_mi = tpMultiIndexMapRef.size(), last_index = num_tp_mi - 1,
    sparse_ind, save_tp_cntr, new_tp_cntr;
  bool map_ref_find;
  for (cit=sparse_indices.begin(); cit!=sparse_indices.end(); ++cit) {
    sparse_ind = *cit; map_ref_find = false;
    for (i=0; i<last_index; ++i)
      // upper bound for i-th multi-index append defined by MapRef[i+1]
      if (sparse_ind < tpMultiIndexMapRef[i+1])
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
	tpMultiIndex[new_tp_cntr] = tpMultiIndex[save_tp_cntr];
      append_multi_index(tpMultiIndex[new_tp_cntr], aggregated_mi,
			 tpMultiIndexMap[new_tp_cntr],
			 tpMultiIndexMapRef[new_tp_cntr]);
    }
    // prune old records off the end
    tpMultiIndex.resize(num_save);
    tpMultiIndexMap.resize(num_save);
    tpMultiIndexMapRef.resize(num_save);
    // update sparse_indices using old_aggregated_mi
    // TO DO: review SharedOrthogPolyApproxData::append_multi_index(SizetSet&)
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
