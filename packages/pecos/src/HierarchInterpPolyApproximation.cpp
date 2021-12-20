/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchInterpPolyApproximation
//- Description:  Implementation code for InterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "HierarchInterpPolyApproximation.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "pecos_stat_util.hpp"

//#define DEBUG
//#define VBD_DEBUG
//#define INTERPOLATION_TEST

namespace Pecos {


void HierarchInterpPolyApproximation::allocate_arrays()
{
  InterpPolyApproximation::allocate_arrays();

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const ExpansionConfigOptions& ec_options = data_rep->expConfigOptions;
  const UShort4DArray& colloc_key = data_rep->hsg_driver()->collocation_key();
  size_t i, j, k, num_levels = colloc_key.size(), num_sets, num_tp_pts,
    num_deriv_v = modSurrData.num_derivative_variables();

  RealVector2DArray& exp_t1_coeffs = expT1CoeffsIter->second;
  RealMatrix2DArray& exp_t2_coeffs = expT2CoeffsIter->second;
  RealMatrix2DArray& exp_t1_coeff_grads = expT1CoeffGradsIter->second;
  
  if (exp_t1_coeffs.size() != num_levels)
    exp_t1_coeffs.resize(num_levels);
  if (exp_t2_coeffs.size() != num_levels)
    exp_t2_coeffs.resize(num_levels);
  if (exp_t1_coeff_grads.size() != num_levels)
    exp_t1_coeff_grads.resize(num_levels);
  for (i=0; i<num_levels; ++i) {
    const UShort3DArray& key_i = colloc_key[i];
    num_sets = key_i.size();
    if (exp_t1_coeffs[i].size() != num_sets)
      exp_t1_coeffs[i].resize(num_sets);
    if (exp_t2_coeffs[i].size() != num_sets)
      exp_t2_coeffs[i].resize(num_sets);
    if (exp_t1_coeff_grads[i].size() != num_sets)
      exp_t1_coeff_grads[i].resize(num_sets);
    for (j=0; j<num_sets; ++j) {
      num_tp_pts = key_i[j].size();
      for (k=0; k<num_tp_pts; ++k) {
	if (expansionCoeffFlag) {
	  exp_t1_coeffs[i][j].sizeUninitialized(num_tp_pts);
	  if (data_rep->basisConfigOptions.useDerivs)
	    exp_t2_coeffs[i][j].shapeUninitialized(num_deriv_v, num_tp_pts);
	}
	if (expansionCoeffGradFlag)
	  exp_t1_coeff_grads[i][j].shapeUninitialized(num_deriv_v, num_tp_pts);
      }
    }
  }

  // checking num_points is insufficient due to anisotropy --> changes in
  // anisotropic weights could move points around without changing the total
  //size_t num_points = modSurrData.points();
  //bool update_exp_form =
  //  ( (expansionCoeffFlag     && exp_t1_coeffs.length() != num_points) ||
  //    (expansionCoeffGradFlag && exp_t1_coeff_grads.numCols() != num_points));

  if (ec_options.refineControl) {
    size_t num_moments = (data_rep->nonRandomIndices.empty()) ? 4 : 2;
    if (refMomentsIter->second.empty())
      refMomentsIter->second.sizeUninitialized(num_moments);
    if (deltaMomentsIter->second.empty())
      deltaMomentsIter->second.sizeUninitialized(num_moments);
  }
}


void HierarchInterpPolyApproximation::compute_coefficients()
{
  PolynomialApproximation::compute_coefficients();
  if (!expansionCoeffFlag && !expansionCoeffGradFlag)
    return;

  allocate_arrays();

  /*
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  if (data_rep->expConfigOptions.discrepancyType == RECURSIVE_DISCREP)
    compute_recursive_coefficients();
  else
    // modSurrData is comprised of discrepancy data (for active keys beyond
    // the first) and we estimate hierarchical coefficients for each key
    // independent of all others
    compute_distinct_coefficients();
  */

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key   = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, v, num_levels = colloc_key.size(), num_sets, num_tp_pts,
    cntr = 0, c_index, num_deriv_vars = modSurrData.num_derivative_variables();
  bool empty_c_index = colloc_index.empty(),
    use_derivs = data_rep->basisConfigOptions.useDerivs;

  RealVector2DArray& exp_t1_coeffs = expT1CoeffsIter->second;
  RealMatrix2DArray& exp_t2_coeffs = expT2CoeffsIter->second;
  RealMatrix2DArray& exp_t1_coeff_grads = expT1CoeffGradsIter->second;

  // always use modified surrogate data to build interpolants for either
  // the base or discrepancy levels (according to the active key)
  const SDVArray& sdv_array = modSurrData.variables_data();
  const SDRArray& sdr_array = modSurrData.response_data();

  // level 0
  c_index = (empty_c_index) ? cntr++ : colloc_index[0][0][0];
  const SurrogateDataResp& sdr_0 = sdr_array[c_index];
  if (expansionCoeffFlag) {
    exp_t1_coeffs[0][0][0] = sdr_0.response_function();
    if (use_derivs)
      Teuchos::setCol(sdr_0.response_gradient(), 0, exp_t2_coeffs[0][0]);
  }
  if (expansionCoeffGradFlag)
    Teuchos::setCol(sdr_0.response_gradient(), 0, exp_t1_coeff_grads[0][0]);

  // levels 1 to num_levels
  for (lev=1; lev<num_levels; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];
    num_sets = key_l.size();
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = key_l[set].size();
      for (pt=0; pt<num_tp_pts; ++pt) {
	c_index = (empty_c_index) ? cntr++ : colloc_index[lev][set][pt];
	const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	const SurrogateDataResp& sdr = sdr_array[c_index];
	// coefficients are hierarchical surpluses
	if (expansionCoeffFlag) {
	  exp_t1_coeffs[lev][set][pt] = sdr.response_function() -
	    value(c_vars, sm_mi, colloc_key, exp_t1_coeffs,exp_t2_coeffs,lev-1);
	  if (use_derivs) {
	    const RealVector& data_grad = sdr.response_gradient();
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, colloc_key, exp_t1_coeffs, exp_t2_coeffs, lev-1);
	    Real* hier_grad = exp_t2_coeffs[lev][set][pt];
	    for (v=0; v<num_deriv_vars; ++v)
	      hier_grad[v] = data_grad[v] - prev_grad[v];
	  }
	}
	if (expansionCoeffGradFlag) {
	  const RealVector& data_grad = sdr.response_gradient();
	  const RealVector& prev_grad = gradient_nonbasis_variables(c_vars,
	    sm_mi, colloc_key, exp_t1_coeff_grads, lev-1);
	  Real* hier_grad = exp_t1_coeff_grads[lev][set][pt];
	  for (v=0; v<num_deriv_vars; ++v)
	    hier_grad[v] = data_grad[v] - prev_grad[v];
	}
      }
    }
  }
#ifdef DEBUG
  PCout << "Hierarchical expansion coeffs (T1 in compute_coefficients()):\n"
	<< exp_t1_coeffs;
#endif
#ifdef INTERPOLATION_TEST
  test_interpolation();
#endif

  // if efficient deltas needed, compute coefficients of product interpolants
  short ref_metric = data_rep->expConfigOptions.refineMetric;
  if (ref_metric == COVARIANCE_METRIC || ref_metric == MIXED_STATS_METRIC)
    initialize_products(); // initialize/update prodType{1,2}Coeffs

  clear_computed_bits();
}


void HierarchInterpPolyApproximation::increment_coefficients()
{
  // TO DO: partial sync for new TP data set, e.g. update_surrogate_data() ?
  synchronize_surrogate_data();

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();

  /* Prior to activeKey tracking for computed bits and cached moments:
  bool updated = update_active_iterators(data_rep->activeKey);
  if (updated) clear_computed_bits();
  else         increment_reference_to_current();
  */
  update_active_iterators(data_rep->activeKey);
  increment_reference_to_current();

  // increment expansionType{1,2}Coeff{s,Grads}
  switch (data_rep->expConfigOptions.refineControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED:
    increment_coefficients(hsg_driver->trial_set());  break;
  default: {
    const UShort3DArray&   sm_mi = hsg_driver->smolyak_multi_index();
    const UShortArray& incr_sets = hsg_driver->increment_sets();
    size_t lev, num_lev = sm_mi.size(), set, start_set, num_sets;
    for (lev=0; lev<num_lev; ++lev) {
      start_set = incr_sets[lev]; num_sets = sm_mi[lev].size();
      for (set=start_set; set<num_sets; ++set)
	increment_coefficients(sm_mi[lev][set]);
    }
    break;
  }
  }
#ifdef DEBUG
  PCout << "Hierarchical expansion coeffs (T1 in increment_coefficients()):\n"
	<< expT1CoeffsIter->second;
#endif

  // size sobolIndices based on shared sobolIndexMap
  allocate_component_sobol();

  // increment productType{1,2}Coeffs, if needed
  if (product_interpolants()) {
    UShort2DArray incr_key;
    if (data_rep->expConfigOptions.refineControl ==
	DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)
      hsg_driver->partition_increment_key(incr_key);
    else
      hsg_driver->increment_sets_to_increment_key(
	hsg_driver->increment_sets(), incr_key);
    increment_products(incr_key);
  }
}


/** Lower level helper function to process a single index set. */
void HierarchInterpPolyApproximation::
increment_coefficients(const UShortArray& index_set)
{
  RealVector2DArray& exp_t1_coeffs = expT1CoeffsIter->second;
  RealMatrix2DArray& exp_t2_coeffs = expT2CoeffsIter->second;
  RealMatrix2DArray& exp_t1_coeff_grads = expT1CoeffGradsIter->second;
  const SDVArray& sdv_array = modSurrData.variables_data();
  const SDRArray& sdr_array = modSurrData.response_data();

  size_t lev, old_levels = exp_t1_coeffs.size(), set, old_sets,
    pt, old_pts = 0;
  for (lev=0; lev<old_levels; ++lev) {
    old_sets = exp_t1_coeffs[lev].size();
    for (set=0; set<old_sets; ++set)
      old_pts += (expansionCoeffFlag) ?	exp_t1_coeffs[lev][set].length() :
	exp_t1_coeff_grads[lev][set].numCols();
  }
  lev = l1_norm(index_set);
  if (lev >= old_levels) {
    exp_t1_coeffs.resize(lev+1);
    exp_t2_coeffs.resize(lev+1);
    exp_t1_coeff_grads.resize(lev+1);
  }
  set = exp_t1_coeffs[lev].size();
  // append empty and update in place
  RealVector fns; RealMatrix grads;
  exp_t1_coeffs[lev].push_back(fns);
  exp_t2_coeffs[lev].push_back(grads);
  exp_t1_coeff_grads[lev].push_back(grads);
  RealVector& t1_coeffs      = exp_t1_coeffs[lev][set];
  RealMatrix& t2_coeffs      = exp_t2_coeffs[lev][set];
  RealMatrix& t1_coeff_grads = exp_t1_coeff_grads[lev][set];

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray& colloc_key = hsg_driver->collocation_key();
  size_t index, num_trial_pts = colloc_key[lev][set].size(),
    v, num_deriv_vars = 0;
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;
  if (expansionCoeffFlag) {
    t1_coeffs.sizeUninitialized(num_trial_pts);
    if (use_derivs) {
      num_deriv_vars = exp_t2_coeffs[0][0].numRows();
      t2_coeffs.shapeUninitialized(num_deriv_vars, num_trial_pts);
    }
  }
  if (expansionCoeffGradFlag) {
    num_deriv_vars = exp_t1_coeff_grads[0][0].numRows();
    t1_coeff_grads.shapeUninitialized(num_deriv_vars, num_trial_pts);
  }
 
  for (pt=0, index=old_pts; pt<num_trial_pts; ++pt, ++index) {
    const RealVector& c_vars = sdv_array[index].continuous_variables();
    const SurrogateDataResp& sdr = sdr_array[index];
    if (expansionCoeffFlag) {
      t1_coeffs[pt] = sdr.response_function() -
	value(c_vars, sm_mi, colloc_key, exp_t1_coeffs, exp_t2_coeffs, lev-1);
      if (use_derivs) {
	const RealVector& data_grad = sdr.response_gradient();
	const RealVector& prev_grad = gradient_basis_variables(c_vars, sm_mi,
	  colloc_key, exp_t1_coeffs, exp_t2_coeffs, lev-1);
	Real* hier_grad = t2_coeffs[pt];
	for (v=0; v<num_deriv_vars; ++v)
	  hier_grad[v] = data_grad[v] - prev_grad[v];
      }
    }
    if (expansionCoeffGradFlag) {
      const RealVector& data_grad = sdr.response_gradient();
      const RealVector& prev_grad = gradient_nonbasis_variables(c_vars, sm_mi,
	colloc_key, exp_t1_coeff_grads, lev-1);
      Real* hier_grad = t1_coeff_grads[pt];
      for (v=0; v<num_deriv_vars; ++v)
	hier_grad[v] = data_grad[v] - prev_grad[v];
    }
  }
}


void HierarchInterpPolyApproximation::pop_coefficients(bool save_data)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;

  /* Prior to activeKey tracking for computed bits and cached moments:
  bool updated = update_active_iterators(key);
  if (updated) clear_computed_bits();
  else         decrement_current_to_reference();
  */
  update_active_iterators(key);
  decrement_current_to_reference();

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  RealVector2DArray& exp_t1c  = expT1CoeffsIter->second;
  RealMatrix2DArray& exp_t2c  = expT2CoeffsIter->second;
  RealMatrix2DArray& exp_t1cg = expT1CoeffGradsIter->second;
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;
  size_t lev, num_lev, set, start_set;
  switch (data_rep->expConfigOptions.refineControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
    // Note: trial_set() is valid so long as expansion decrement precedes
    // grid decrement (reversing order from increment + update/push), but
    // trial_level() is less volatile / safer.
    lev = hsg_driver->trial_level();
    if (expansionCoeffFlag) {
      if (save_data) {
	RealVectorDequeArray& pop_t1c = poppedExpT1Coeffs[key];
	if (pop_t1c.size() <= lev) pop_t1c.resize(lev+1);
	push_back_to_back(exp_t1c[lev], pop_t1c[lev]);
      }
      else exp_t1c[lev].pop_back();
      if (use_derivs) {
	if (save_data) {
	  RealMatrixDequeArray& pop_t2c = poppedExpT2Coeffs[key];
	  if (pop_t2c.size() <= lev) pop_t2c.resize(lev+1);
	  push_back_to_back(exp_t2c[lev], pop_t2c[lev]);
	}
	else exp_t2c[lev].pop_back();
      }
    }
    if (expansionCoeffGradFlag) {
      if (save_data) {
	RealMatrixDequeArray& pop_t1cg = poppedExpT1CoeffGrads[key];
	if (pop_t1cg.size() <= lev) pop_t1cg.resize(lev+1);
	push_back_to_back(exp_t1cg[lev], pop_t1cg[lev]);
      }
      else exp_t1cg[lev].pop_back();
    }
    break;
  }
  default: {
    // See note above about ordering of expansion / grid decrements
    const UShortArray& incr_sets = data_rep->hsg_driver()->increment_sets();
    num_lev = incr_sets.size();
    // save_data is generally on, so don't bother avoiding lookup if off
    RealVectorDequeArray& pop_t1c  = poppedExpT1Coeffs[key];
    RealMatrixDequeArray& pop_t2c  = poppedExpT2Coeffs[key];
    RealMatrixDequeArray& pop_t1cg = poppedExpT1CoeffGrads[key];
    if (save_data && expansionCoeffFlag) {
      if (pop_t1c.size() <= num_lev) pop_t1c.resize(num_lev+1);
      if (data_rep->basisConfigOptions.useDerivs && pop_t2c.size() <= num_lev)
	pop_t2c.resize(num_lev+1);
    }
    if (save_data && expansionCoeffGradFlag && pop_t1cg.size() <= num_lev)
      pop_t1cg.resize(num_lev+1);
    for (lev=0; lev<num_lev; ++lev) {
      start_set = incr_sets[lev];
      if (save_data) {
	if (expansionCoeffFlag) {
	  push_range_to_back(exp_t1c[lev], start_set, pop_t1c[lev]);
	  if (use_derivs)
	    push_range_to_back(exp_t2c[lev], start_set, pop_t2c[lev]);
	}
	if (expansionCoeffGradFlag)
	  push_range_to_back(exp_t1cg[lev], start_set, pop_t1cg[lev]);
      }
      else {
	if (expansionCoeffFlag) {
	  exp_t1c[lev].resize(start_set);
	  if (use_derivs) exp_t2c[lev].resize(start_set);
	}
	if (expansionCoeffGradFlag)
	  exp_t1cg[lev].resize(start_set);
      }
    }
    break;
  }
  }
#ifdef DEBUG
  PCout << "Hierarchical expansion coeffs (T1 in pop_coefficients()):\n"
	<< exp_t1c;
#endif

  // decrement productType{1,2}Coeffs if in use
  // Note: we define all entries in productType{1,2}Coeffs for convenience
  // even if use_derivs is inactive
  if (product_interpolants()) {
    std::map<PolynomialApproximation*, RealVector2DArray>& prod_t1c
      = prodT1CoeffsIter->second;
    std::map<PolynomialApproximation*, RealMatrix2DArray>& prod_t2c
      = prodT2CoeffsIter->second;
    std::map<PolynomialApproximation*, RealVector2DArray>::iterator e1_it;
    std::map<PolynomialApproximation*, RealMatrix2DArray>::iterator e2_it;
    if (save_data) {
      std::map<PolynomialApproximation*,RealVectorDequeArray>& pop_prod_t1c
	= poppedProdType1Coeffs[key];
      if (pop_prod_t1c.empty()) // sync pointer set with prod_t1c
	for (e1_it = prod_t1c.begin(); e1_it != prod_t1c.end(); ++e1_it)
	  pop_prod_t1c.insert(std::pair<PolynomialApproximation*,
	    RealVectorDequeArray>(e1_it->first, RealVectorDequeArray()));
      std::map<PolynomialApproximation*,RealVectorDequeArray>::iterator p1_it
	= pop_prod_t1c.begin();
      std::map<PolynomialApproximation*,RealMatrixDequeArray>::iterator p2_it;
      if (use_derivs) {
	std::map<PolynomialApproximation*, RealMatrixDequeArray>& pop_prod_t2c
	  = poppedProdType2Coeffs[key];
	if (pop_prod_t2c.empty()) // sync pointer set with prod_t2c
	  for (e2_it = prod_t2c.begin(); e2_it != prod_t2c.end(); ++e2_it)
	    pop_prod_t2c.insert(std::pair<PolynomialApproximation*,
	      RealMatrixDequeArray>(e2_it->first, RealMatrixDequeArray()));
	e2_it = prod_t2c.begin(); p2_it = pop_prod_t2c.begin();
      }
      switch (data_rep->expConfigOptions.refineControl) {
      case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED:
	for (e1_it = prod_t1c.begin(); e1_it != prod_t1c.end(); ++e1_it) {
	  PolynomialApproximation* poly_approx_2 = e1_it->first;
	  RealVectorDequeArray& pop_prod_t1c = p1_it->second;
	  if (pop_prod_t1c.size() <= lev) pop_prod_t1c.resize(lev+1);
	  push_back_to_back(e1_it->second[lev], pop_prod_t1c[lev]);
	  ++p1_it;  //++e1_it;
	  if (use_derivs) {
	    RealMatrixDequeArray& pop_prod_t2c = p2_it->second;
	    if (pop_prod_t2c.size() <= lev) pop_prod_t2c.resize(lev+1);
	    push_back_to_back(e2_it->second[lev], pop_prod_t2c[lev]);
	    ++p2_it;  ++e2_it;
	  }
	}
	break;
      default: {
	const UShortArray& incr_sets = data_rep->hsg_driver()->increment_sets();
	for (e1_it=prod_t1c.begin(); e1_it!=prod_t1c.end(); ++e1_it, ++p1_it) {
	  if (p1_it->second.size() <= num_lev) p1_it->second.resize(num_lev+1);
	  if (use_derivs && p2_it->second.size() <= num_lev)
	    p2_it->second.resize(num_lev+1);
	  for (lev=0; lev<num_lev; ++lev) {
	    start_set = incr_sets[lev];
	    push_range_to_back(e1_it->second[lev], start_set,
			       p1_it->second[lev]);
	    if (use_derivs)
	      push_range_to_back(e2_it->second[lev], start_set,
				 p2_it->second[lev]);
	  }
	  if (use_derivs) { ++p2_it; ++e2_it; }
	}
	break;
      }
      }
    }
    else { // don't save_data
      if (use_derivs) e2_it = prod_t2c.begin();
      switch (data_rep->expConfigOptions.refineControl) {
      case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED:
	for (e1_it = prod_t1c.begin(); e1_it != prod_t1c.end(); ++e1_it) {
	  e1_it->second[lev].pop_back();
	  if (use_derivs) { e2_it->second[lev].pop_back(); ++e2_it; }
	}
	break;
      default: {
	const UShortArray& incr_sets = data_rep->hsg_driver()->increment_sets();
	for (e1_it = prod_t1c.begin(); e1_it != prod_t1c.end(); ++e1_it) {
	  for (lev=0; lev<num_lev; ++lev) {
	    start_set = incr_sets[lev];
	    e1_it->second[lev].resize(start_set);
	    if (use_derivs) { e2_it->second[lev].resize(start_set); ++e2_it; }
	  }
	}
	break;
      }
      }
    }
  }
}


void HierarchInterpPolyApproximation::push_coefficients()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;

  /* Prior to activeKey tracking for computed bits and cached moments:
  bool updated = update_active_iterators(key);
  if (updated) clear_computed_bits();
  else         increment_reference_to_current();
  */
  update_active_iterators(key);
  increment_reference_to_current();

  switch (data_rep->expConfigOptions.refineControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
    // Note: for both a restored or a selected candidate, trial set is pushed
    // to Driver and corresponds to sm_mi[lev].back(), but in the latter case,
    // must locate corresponding popped data within interior to push onto back
    // of coeff arrays --> adopt SharedPolyApproxData::push_index() approach
    // > pushIndex is defined from *SparseGridDriver::poppedLevMultiIndex.
    // > nested maps approach requires a valid map key and trial_set() is
    //   a problem in pop_coefficients() since it is downstream (pushIndex
    //   is defined upstream in shared data) --> would have to retrieve from
    //   poppedLevMultiIndex somehow (too late?)

    //HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
    //const UShortArray& tr_set = hsg_driver->trial_set();
    //size_t lev = l1_norm(tr_set);
    size_t tr_lev = data_rep->hsg_driver()->trial_level(),
          p_index = data_rep->push_index();
    bool use_derivs = data_rep->basisConfigOptions.useDerivs;
    RealVectorDeque::iterator v_it;  RealMatrixDeque::iterator m_it;
    if (expansionCoeffFlag) {
      push_index_to_back(poppedExpT1Coeffs[key][tr_lev], p_index,
			 expT1CoeffsIter->second[tr_lev]);
      SharedHierarchInterpPolyApproxData* data_rep
	= (SharedHierarchInterpPolyApproxData*)sharedDataRep;
      if (use_derivs)
	push_index_to_back(poppedExpT2Coeffs[key][tr_lev], p_index,
			   expT2CoeffsIter->second[tr_lev]);
    }
    if (expansionCoeffGradFlag)
      push_index_to_back(poppedExpT1CoeffGrads[key][tr_lev], p_index,
			 expT1CoeffGradsIter->second[tr_lev]);
#ifdef DEBUG
    PCout << "Hierarchical expansion coeffs (T1 in push_coefficients()):\n"
	  << expT1CoeffsIter->second;
#endif

    // update productType{1,2}Coeffs if in use
    if (product_interpolants()) {
      // Note: we define all entries in productType{1,2}Coeffs for convenience
      // even if use_derivs is inactive
      std::map<PolynomialApproximation*, RealVector2DArray>& prod_t1c
	= prodT1CoeffsIter->second;
      std::map<PolynomialApproximation*, RealVectorDequeArray>& pop_prod_t1c
	= poppedProdType1Coeffs[key];
      std::map<PolynomialApproximation*, RealVector2DArray>::iterator    e1_it;
      std::map<PolynomialApproximation*, RealMatrix2DArray>::iterator    e2_it;
      std::map<PolynomialApproximation*, RealVectorDequeArray>::iterator p1_it;
      std::map<PolynomialApproximation*, RealMatrixDequeArray>::iterator p2_it;
      if (use_derivs) {
	e2_it = prodT2CoeffsIter->second.begin();
	p2_it = poppedProdType2Coeffs[key].begin();
      }
      for (e1_it  = prod_t1c.begin(), p1_it  = pop_prod_t1c.begin();
	   e1_it != prod_t1c.end() && p1_it != pop_prod_t1c.end();
	   ++e1_it, ++p1_it) {
	push_index_to_back(p1_it->second[tr_lev], p_index,
			   e1_it->second[tr_lev]);
	if (use_derivs) {
	  push_index_to_back(p2_it->second[tr_lev], p_index,
			     e2_it->second[tr_lev]);
	  ++e2_it; ++p2_it;
	}
      }
    }
    break;
  }
  default: // only one candidate for iso/aniso refinement: push all of its sets
    promote_all_popped_coefficients();
    break;
  }
}


void HierarchInterpPolyApproximation::finalize_coefficients()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;

  /* Prior to activeKey tracking for computed bits and cached moments:
  bool updated = update_active_iterators(data_rep->activeKey);
  if (updated) clear_computed_bits();
  else         clear_current_computed_bits(); //clear_delta_computed_bits();
  */
  // synchronize expansionCoeff{s,Grads} and approxData
  update_active_iterators(data_rep->activeKey);
  clear_current_computed_bits(); //clear_delta_computed_bits();

  promote_all_popped_coefficients();
}


void HierarchInterpPolyApproximation::promote_all_popped_coefficients()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;

  // Code supports iso/aniso refinement (one candidate with multiple index sets)
  // or generalized refinement (multiple candidates, each with one index set)

  RealVector2DArray&    exp_t1c  = expT1CoeffsIter->second;
  RealMatrix2DArray&    exp_t2c  = expT2CoeffsIter->second;
  RealMatrix2DArray&    exp_t1cg = expT1CoeffGradsIter->second;
  RealVectorDequeArray& pop_t1c  = poppedExpT1Coeffs[key];
  RealMatrixDequeArray& pop_t2c  = poppedExpT2Coeffs[key];
  RealMatrixDequeArray& pop_t1cg = poppedExpT1CoeffGrads[key];
  size_t lev, num_lev = std::max(exp_t1c.size(), exp_t1cg.size());
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;
  for (lev=0; lev<num_lev; ++lev) {
    if (expansionCoeffFlag) {
      push_range_to_back(pop_t1c[lev], 0, exp_t1c[lev]);
      if (use_derivs) push_range_to_back(pop_t2c[lev], 0, exp_t2c[lev]);
    }
    if (expansionCoeffGradFlag)
      push_range_to_back(pop_t1cg[lev], 0, exp_t1cg[lev]);
  }
#ifdef DEBUG
  PCout << "Hierarchical expansion coeffs (T1 in promote_all_popped_"
	<< "coefficients()):\n"	<< expT1CoeffsIter->second;
#endif

  // update productType{1,2}Coeffs if in use
  if (product_interpolants()) {
    std::map<PolynomialApproximation*, RealVector2DArray>&        prod_t1c
      = prodT1CoeffsIter->second;
    std::map<PolynomialApproximation*, RealVectorDequeArray>& pop_prod_t1c
      = poppedProdType1Coeffs[key];
    std::map<PolynomialApproximation*, RealVector2DArray>::iterator    e1_it;
    std::map<PolynomialApproximation*, RealVectorDequeArray>::iterator p1_it;
    std::map<PolynomialApproximation*, RealMatrix2DArray>::iterator    e2_it;
    std::map<PolynomialApproximation*, RealMatrixDequeArray>::iterator p2_it;
    if (use_derivs) {
      e2_it = prodT2CoeffsIter->second.begin();
      p2_it = poppedProdType2Coeffs[key].begin();
    }
    for (e1_it  = prod_t1c.begin(), p1_it  = pop_prod_t1c.begin();
	 e1_it != prod_t1c.end() && p1_it != pop_prod_t1c.end();
	 ++e1_it, ++p1_it) {
      for (lev=0; lev<num_lev; ++lev) {
	push_range_to_back(p1_it->second[lev], 0, e1_it->second[lev]);
	if (use_derivs)
	  push_range_to_back(p2_it->second[lev], 0, e2_it->second[lev]);
      }
      if (use_derivs) { ++e2_it; ++p2_it; }
    }
  }
}


void HierarchInterpPolyApproximation::clear_inactive()
{
  std::map<UShortArray, RealVector2DArray>::iterator e1c_it
    = expansionType1Coeffs.begin();
  std::map<UShortArray, RealMatrix2DArray>::iterator e2c_it
    = expansionType2Coeffs.begin();
  std::map<UShortArray, RealMatrix2DArray>::iterator e1g_it
    = expansionType1CoeffGrads.begin();
  while (e1c_it != expansionType1Coeffs.end())
    if (e1c_it == expT1CoeffsIter) // preserve active
      { ++e1c_it; ++e2c_it; ++e1g_it; }
    else { // clear inactive: postfix increments manage iterator invalidations
      expansionType1Coeffs.erase(e1c_it++);
      expansionType2Coeffs.erase(e2c_it++);
      expansionType1CoeffGrads.erase(e1g_it++);
    }
}


void HierarchInterpPolyApproximation::combine_coefficients()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();

  // Coefficient combination is not dependent on active state
  //update_active_iterators(data_rep->activeKey);

  // Note: computed bits are also cleared when refineStatsType is changed
  if (data_rep->expConfigOptions.refineStatsType == COMBINED_EXPANSION_STATS)
    clear_computed_bits(); // increment: clear_current_computed_bits() ?

  allocate_component_sobol(); // size sobolIndices from shared sobolIndexMap

  const UShort4DArray&   comb_key = hsg_driver->combined_collocation_key();
  const Sizet3DArray& comb_sm_map
    = hsg_driver->combined_smolyak_multi_index_map();
  size_t i, lev, set, pt, num_lev = comb_key.size(), num_sets, num_tp_pts,
    num_v = modSurrData.num_derivative_variables();
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;

  // Resize combined expansion arrays and initialize to zero.
  // Note: STL resize() is no-op if already correct size, but Teuchos
  // resize()/reshape() is not; therefore protect the latter.
  combinedExpT1Coeffs.resize(num_lev); combinedExpT2Coeffs.resize(num_lev);
  combinedExpT1CoeffGrads.resize(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    const UShort3DArray& comb_key_l = comb_key[lev];
    num_sets = comb_key_l.size();
    RealVectorArray& comb_t1c_l = combinedExpT1Coeffs[lev];
    RealMatrixArray& comb_t2c_l = combinedExpT2Coeffs[lev];
    RealMatrixArray& comb_t1g_l = combinedExpT1CoeffGrads[lev];
    comb_t1c_l.resize(num_sets);  comb_t2c_l.resize(num_sets);
    comb_t1g_l.resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = comb_key_l[set].size();
      RealVector& comb_t1c_ls = comb_t1c_l[set];
      RealMatrix& comb_t2c_ls = comb_t2c_l[set];
      RealMatrix& comb_t1g_ls = comb_t1g_l[set];
      if (expansionCoeffFlag) {
	if (comb_t1c_ls.length() != num_tp_pts)
	  comb_t1c_ls.size(num_tp_pts);         // init to 0
	else comb_t1c_ls = 0.;
      }
      else if (!comb_t1c_ls.empty())
	comb_t1c_ls.size(0);
      if (use_derivs && expansionCoeffFlag) {
	if (comb_t2c_ls.numRows()!=num_v || comb_t2c_ls.numCols()!=num_tp_pts)
	  comb_t2c_ls.shape(num_v, num_tp_pts); // init to 0
	else comb_t2c_ls = 0.;
      }
      else if (!comb_t2c_ls.empty())
	comb_t2c_ls.shape(0, 0);
      if (expansionCoeffGradFlag) {
	if (comb_t1g_ls.numRows()!=num_v || comb_t1g_ls.numCols()!=num_tp_pts)
	  comb_t1g_ls.shape(num_v, num_tp_pts); // init to 0
	else comb_t1g_ls = 0.;
      }
      else if (!comb_t1g_ls.empty())
	comb_t1g_ls.shape(0, 0);
    }
  }

  // Roll up hierarchical surplus increments for each level-set-point
  if (expansionCoeffFlag) {
    std::map<UShortArray, RealVector2DArray>::iterator ec1_it;
    std::map<UShortArray, RealMatrix2DArray>::iterator ec2_it;
    for (ec1_it  = expansionType1Coeffs.begin(),
	 ec2_it  = expansionType2Coeffs.begin(), i=0;
	 ec1_it != expansionType1Coeffs.end() &&
	 ec2_it != expansionType2Coeffs.end(); ++ec1_it, ++ec2_it, ++i) {
      const RealVector2DArray& t1c = ec1_it->second;
      const RealMatrix2DArray& t2c = ec2_it->second;
      num_lev = t1c.size();
      for (lev=0; lev<num_lev; ++lev) {
	const RealVectorArray& t1c_l = t1c[lev];
	const RealMatrixArray& t2c_l = t2c[lev];
	RealVectorArray&  comb_t1c_l = combinedExpT1Coeffs[lev];
	RealMatrixArray&  comb_t2c_l = combinedExpT2Coeffs[lev];
	const SizetArray& comb_sm_map_il = comb_sm_map[i][lev];
	num_sets = t1c_l.size();
	for (set=0; set<num_sets; ++set) {
	  // map from smolyakMultiIndex for this key to corresponding set
	  // in combinedSmolyakMultiIndex
	  size_t comb_set = comb_sm_map_il[set];
	  comb_t1c_l[comb_set]                 += t1c_l[set];
	  if (use_derivs) comb_t2c_l[comb_set] += t2c_l[set];
	}
      }
    }
  }
  if (expansionCoeffGradFlag) { // sum up expansionType1CoeffGrads
    std::map<UShortArray, RealMatrix2DArray>::iterator eg1_it;
    for (eg1_it  = expansionType1CoeffGrads.begin(), i=0;
	 eg1_it != expansionType1CoeffGrads.end(); ++eg1_it, ++i) {
      const RealMatrix2DArray& t1g = eg1_it->second;
      num_lev = t1g.size();
      for (lev=0; lev<num_lev; ++lev) {
	const RealMatrixArray& t1g_l = t1g[lev];
	RealMatrixArray&  comb_t1g_l = combinedExpT1CoeffGrads[lev];
	const SizetArray& comb_sm_map_il = comb_sm_map[i][lev];
	num_sets = t1g_l.size();
	for (set=0; set<num_sets; ++set) {
	  // map from smolyakMultiIndex for this key to corresponding set
	  // in combinedSmolyakMultiIndex
	  size_t comb_set = comb_sm_map_il[set];
	  comb_t1g_l[comb_set] += t1g_l[set];
	}
      }
    }
  }
}


void HierarchInterpPolyApproximation::combined_to_active(bool clear_combined)
{
  // replace active expansions with combined expansion arrays
  // > clear_inactive() takes care of the auxilliary inactive expansions
  //   that are now assimilated within the active expansion
  // > we reassign T1Coeffs,T2Coeffs,T1CoeffGrads even if not active in order
  //   to preserve hierarchical sizing needed downstream (in value() et al.)

  if (clear_combined) {
    std::swap(expT1CoeffsIter->second,     combinedExpT1Coeffs);
    std::swap(expT2CoeffsIter->second,     combinedExpT2Coeffs);
    std::swap(expT1CoeffGradsIter->second, combinedExpT1CoeffGrads);
    combinedExpT1Coeffs.clear();  combinedExpT2Coeffs.clear();
    combinedExpT1CoeffGrads.clear();
  }
  else { // (redundant) copies
    expT1CoeffsIter->second     = combinedExpT1Coeffs;
    expT2CoeffsIter->second     = combinedExpT2Coeffs;
    expT1CoeffGradsIter->second = combinedExpT1CoeffGrads;
  }

  // clear accumulated (raw) product coefficients
  // (used to accelerate delta covariance adaptations)
  productType1Coeffs.clear();     prodT1CoeffsIter = productType1Coeffs.end();
  productType2Coeffs.clear();     prodT2CoeffsIter = productType2Coeffs.end();
  poppedProdType1Coeffs.clear();  poppedProdType2Coeffs.clear();

  // Create a dummy modSurrData for the combined-now-active coeffs, for
  // accelerating FINAL_RESULTS (integration, VBD processing, etc.)
  synthetic_surrogate_data(modSurrData); // overwrite data for activeKey

  // If outgoing stats type is active (e.g., as in Dakota::NonDExpansion::
  // multifidelity_expansion()), then previous active stats are invalidated.
  // But if outgoing stats type is combined, then can avoid recomputation
  // and carry over current moment stats from combined to active. 
  // Note: due to this carry-over optimization, updating of stats type from
  //       COMBINED to ACTIVE must follow this function
  // Note: reference and delta are less important to reuse in the context of
  //       final active processing, so clear these tracker bits for simplicity
  //       (even though they could be preserved as well with sufficient care).
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  if (data_rep->expConfigOptions.refineStatsType == COMBINED_EXPANSION_STATS)
    { clear_reference_computed_bits(); clear_delta_computed_bits(); }
  else // previous active coeffs overwritten -> all moments invalidated
    clear_computed_bits();
}


void HierarchInterpPolyApproximation::
synthetic_surrogate_data(SurrogateData& surr_data)
{
  // Update the active key of surr_data with synthetic data based on the
  // active grid from hsg_driver

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();

  // HierarchSparseGridDriver::combined_to_active() transfers all data except
  // collocation indices, which are invalidated by the combination.  In support
  // of the synthetic data to be created, new colloc point count and (default)
  // colloc indices are defined at the end of HSGDriver::combined_to_active(),
  // and this default colloc index sequence is employed below. 
  const RealMatrix2DArray& var_sets = hsg_driver->hierarchical_variable_sets();
  const UShort3DArray&        sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray&   colloc_key = hsg_driver->collocation_key();
  const Sizet3DArray&  colloc_index = hsg_driver->collocation_indices();
  size_t             num_colloc_pts = hsg_driver->collocation_points(),
    lev, num_lev = colloc_key.size(), set, num_sets, pt, num_tp_pts,
    num_v = var_sets[0][0].numRows();

  const RealVector2DArray& t1_coeffs = expT1CoeffsIter->second;
  const RealMatrix2DArray& t2_coeffs = expT2CoeffsIter->second;

  // shallow copy vars array data using HierarchSparseGridDriver::variableSets
  surr_data.clear_all_active();
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;
  short bits = (use_derivs) ? 3 : 1;
  surr_data.resize(num_colloc_pts, bits, num_v);

  // use interpolant to produce data values that are being interpolated
  // Note: for a nested hierarchical construction, we know that contributions
  //       from level > l are zero for colloc pts that correspond to level l.
  SDVArray& sdv_array = surr_data.variables_data(); 
  SDRArray& sdr_array = surr_data.response_data();
  size_t c_index = colloc_index[0][0][0];
  RealVector v0(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),(int)num_v);
  sdv_array[c_index].continuous_variables(v0); // DEFAULT_COPY assumes view
  sdr_array[c_index].response_function(
    value(v0, sm_mi, colloc_key, t1_coeffs, t2_coeffs, 0));
  if (use_derivs)
    sdr_array[c_index].response_gradient(
      gradient_basis_variables(v0, sm_mi, colloc_key, t1_coeffs, t2_coeffs, 0));
  for (lev=1; lev<num_lev; ++lev) {
    num_sets = colloc_key[lev].size();
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      const RealMatrix& var_sets_ls = var_sets[lev][set];
      for (pt=0; pt<num_tp_pts; ++pt) {
	c_index = colloc_index[lev][set][pt];
	RealVector v_lsp(Teuchos::View, const_cast<Real*>(var_sets_ls[pt]),
			 (int)num_v);
	sdv_array[c_index].continuous_variables(v_lsp);
	sdr_array[c_index].response_function(
	  value(v_lsp, sm_mi, colloc_key, t1_coeffs, t2_coeffs, lev));
	if (use_derivs)
	  sdr_array[c_index].response_gradient(
	    gradient_basis_variables(v_lsp, sm_mi, colloc_key,
				     t1_coeffs, t2_coeffs, lev));
      }
    }
  }
}


void HierarchInterpPolyApproximation::initialize_products()
{
  std::map<PolynomialApproximation*, RealVector2DArray>& prod_t1c
    = prodT1CoeffsIter->second;
  std::map<PolynomialApproximation*, RealMatrix2DArray>& prod_t2c
    = prodT2CoeffsIter->second;

  size_t num_cov_ptr = covariancePointers.size();
  if (prod_t1c.size() == num_cov_ptr && prod_t2c.size() == num_cov_ptr) {
    // already configured; only need to clear data
    std::map<PolynomialApproximation*, RealVector2DArray>::iterator it1;
    std::map<PolynomialApproximation*, RealMatrix2DArray>::iterator it2;
    for (it1  = prod_t1c.begin(), it2  = prod_t2c.begin();
	 it1 != prod_t1c.end() && it2 != prod_t2c.end(); ++it1, ++it2)
      { it1->second.clear(); it2->second.clear(); }
  }
  else { // build pointer mappings from sratch
    prod_t1c.clear();  prod_t2c.clear();
    RealVector2DArray empty_rv2a;  RealMatrix2DArray empty_rm2a;
    std::set<PolynomialApproximation*>::iterator it;
    for (it=covariancePointers.begin(); it!=covariancePointers.end(); ++it)
      { prod_t1c[*it] = empty_rv2a;  prod_t2c[*it] = empty_rm2a; }
  }

  // synchronize product interpolants with current state of expType{1,2}Coeffs
  increment_products(); // no set_partition -> compute all active terms
}


void HierarchInterpPolyApproximation::
increment_products(const UShort2DArray& set_partition)
{
  // update coefficients of product interpolants needed for efficient delta
  // covariance calculations

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  std::map<PolynomialApproximation*, RealVector2DArray>& prod_t1c
    = prodT1CoeffsIter->second;
  std::map<PolynomialApproximation*, RealMatrix2DArray>& prod_t2c
    = prodT2CoeffsIter->second;
  std::map<PolynomialApproximation*, RealVector2DArray>::iterator it1;
  std::map<PolynomialApproximation*, RealMatrix2DArray>::iterator it2;
  // loop over all PolynomialApproximation* instances previously initialized
  // (including this pointer)
  if (data_rep->expConfigOptions.refineStatsType == COMBINED_EXPANSION_STATS) {
    UShortArray lf_key;  paired_lf_key(data_rep->activeKey, lf_key);
    for (it1  = prod_t1c.begin(), it2  = prod_t2c.begin();
	 it1 != prod_t1c.end() && it2 != prod_t2c.end(); ++it1, ++it2) {
      product_difference_interpolant(
	(HierarchInterpPolyApproximation*)it1->first, it1->second,
	it2->second, lf_key, set_partition);
    }
  }
  else
    for (it1  = prod_t1c.begin(), it2  = prod_t2c.begin();
	 it1 != prod_t1c.end() && it2 != prod_t2c.end(); ++it1, ++it2)
      product_interpolant((HierarchInterpPolyApproximation*)it1->first,
	it1->second, it2->second, set_partition);
}


Real HierarchInterpPolyApproximation::
value(const RealVector& x, const UShort3DArray& sm_mi,
      const UShort4DArray& colloc_key, const RealVector2DArray& t1_coeffs,
      const RealMatrix2DArray& t2_coeffs, unsigned short level,
      const UShort2DArray& set_partition)
{
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real approx_val = 0.;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  size_t lev, set, set_start = 0, set_end;
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = colloc_key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeffs_l.size();
    for (set=set_start; set<set_end; ++set)
      approx_val +=
	data_rep->tensor_product_value(x, t1_coeffs_l[set], t2_coeffs_l[set],
				       sm_mi_l[set], key_l[set], colloc_index);
  }
  return approx_val;
}


/** All variables version. */
Real HierarchInterpPolyApproximation::
value(const RealVector& x, const UShort3DArray& sm_mi,
      const UShort4DArray& colloc_key, const RealVector2DArray& t1_coeffs,
      const RealMatrix2DArray& t2_coeffs, unsigned short level,
      const SizetList& subset_indices, const UShort2DArray& set_partition)
{
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real approx_val = 0.;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  size_t lev, set, set_start = 0, set_end;
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = colloc_key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeffs_l.size();
    for (set=set_start; set<set_end; ++set)
      approx_val +=
	data_rep->tensor_product_value(x, t1_coeffs_l[set], t2_coeffs_l[set],
				       sm_mi_l[set], key_l[set], colloc_index,
				       subset_indices);
  }
  return approx_val;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& colloc_key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 unsigned short level,
			 const UShort2DArray& set_partition)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t num_v = sharedDataRep->numVars;
  if (approxGradient.length() != num_v)
    approxGradient.sizeUninitialized(num_v);
  approxGradient = 0.;

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  size_t lev, set, set_start = 0, set_end;
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = colloc_key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeffs_l.size();
    for (set=set_start; set<set_end; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& colloc_key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 unsigned short level, const SizetList& subset_indices,
			 const UShort2DArray& set_partition)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t num_v = sharedDataRep->numVars;
  if (approxGradient.length() != num_v)
    approxGradient.sizeUninitialized(num_v);
  approxGradient = 0.;

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  size_t lev, set, set_start = 0, set_end;
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = colloc_key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeffs_l.size();
    for (set=set_start; set<set_end; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index,
	  subset_indices);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& colloc_key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 const SizetArray& dvv, unsigned short level,
			 const UShort2DArray& set_partition)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t lev, set, set_start = 0, set_end, num_deriv_v = dvv.size();
  if (approxGradient.length() != num_deriv_v)
    approxGradient.sizeUninitialized(num_deriv_v);
  approxGradient = 0.;

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = colloc_key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeffs_l.size();
    for (set=set_start; set<set_end; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index, dvv);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			    const UShort4DArray& colloc_key,
			    const RealMatrix2DArray& t1_coeff_grads,
			    unsigned short level,
			    const UShort2DArray& set_partition)
{
  // Error check for required data
  size_t lev, set, set_start = 0, set_end, num_deriv_v;
  if (expansionCoeffGradFlag) {
    if (t1_coeff_grads.size() > level && t1_coeff_grads[level].size())
      num_deriv_v = t1_coeff_grads[level][0].numRows();
    else {
      PCerr << "Error: insufficient size in type1 expansion coefficient "
	    << "gradients in\n       HierarchInterpPolyApproximation::"
	    << "gradient_nonbasis_variables()" << std::endl;
      abort_handler(-1);
    }
  }
  else {
    PCerr << "Error: expansion coefficient gradients not defined in Hierarch"
	  << "InterpPolyApproximation::gradient_nonbasis_variables()"
	  << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != num_deriv_v)
    approxGradient.sizeUninitialized(num_deriv_v);
  approxGradient = 0.;

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  bool partial = !set_partition.empty();
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&            sm_mi_l = sm_mi[lev];
    const UShort3DArray&              key_l = colloc_key[lev];
    const RealMatrixArray& t1_coeff_grads_l = t1_coeff_grads[lev];
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else // sm_mi and key include all current index sets, whereas the t1/t2
         // coeffs may refect a partial state derived from a ref_key
      set_end = t1_coeff_grads_l.size();
    for (set=set_start; set<set_end; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_nonbasis_variables(x,
	  t1_coeff_grads_l[set], sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


const RealSymMatrix& HierarchInterpPolyApproximation::
hessian_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			const UShort4DArray& colloc_key,
			const RealVector2DArray& t1_coeffs,
			unsigned short level,
			const UShort2DArray& set_partition)
{
  PCerr << "Error: HierarchInterpPolyApproximation::hessian_basis_variables() "
	<< "not yet implemented." << std::endl;
  abort_handler(-1);

  return approxHessian;
}


Real HierarchInterpPolyApproximation::mean()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compMeanIter->second & 1))
    return numMomentsIter->second[0];

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  Real mean = expectation(expT1CoeffsIter->second, expT2CoeffsIter->second);
  if (std_mode)
    { numMomentsIter->second[0] = mean; compMeanIter->second |= 1; }
  return mean;
}



Real HierarchInterpPolyApproximation::mean(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean[key]))
    return numMomentsIter->second[0];

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real mean = expectation(x, expT1CoeffsIter->second, expT2CoeffsIter->second);
  if (all_mode) {
    numMomentsIter->second[0] = mean;
    compMeanIter->second |= 1;  xPrevMean[key] = x;
  }
  return mean;
}


Real HierarchInterpPolyApproximation::combined_mean()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compMeanIter->second & 1))
    return numMomentsIter->second[0];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real mean =
    expectation(combinedExpT1Coeffs, combinedExpT2Coeffs,
		hsg_driver->combined_type1_hierarchical_weight_sets(),
		hsg_driver->combined_type2_hierarchical_weight_sets());
  if (std_mode)
    { numMomentsIter->second[0] = mean; compMeanIter->second |= 1; }
  return mean;
}



Real HierarchInterpPolyApproximation::combined_mean(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean[key]))
    return numMomentsIter->second[0];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real mean = expectation(x, combinedExpT1Coeffs, combinedExpT2Coeffs,
			  hsg_driver->combined_smolyak_multi_index(),
			  hsg_driver->combined_collocation_key());
  if (all_mode) {
    numMomentsIter->second[0] = mean;
    compMeanIter->second |= 1;  xPrevMean[key] = x;
  }
  return mean;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& HierarchInterpPolyApproximation::mean_gradient()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (std_mode && (compMeanIter->second & 2))
    return momentGradsIter->second[0];

  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Hierarch"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  RealVector& mean_grad = momentGradsIter->second[0];
  mean_grad = expectation_gradient(expT1CoeffGradsIter->second);
  if (std_mode) compMeanIter->second |=  2;//  activate 2-bit
  else          compMeanIter->second &= ~2;//deactivate 2-bit: protect mixed use
  return mean_grad;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  In this case, the mean of the expansion is the
    expectation over the random subset and the derivative of the mean
    is the derivative of the remaining expansion over the non-random
    subset.  This function must handle the mixed case, where some
    design/state variables are augmented (and are part of the
    expansion: derivatives are evaluated as described above) and some
    are inserted (derivatives are obtained from expansionType1CoeffGrads). */
const RealVector& HierarchInterpPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  // if already computed, return previous result
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if ( all_mode && (compMeanIter->second & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevMeanGrad[key]) )
    // && dvv == dvvPrev)
    return momentGradsIter->second[0];

  // ---------------------------------------------------------------------
  // For xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   mu(sa, si)    = Sum_i r_i(sa, si) wt_prod_i
  //   dmu/ds        = Sum_i dr_i/ds wt_prod_i
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   mu(sa, si)    = Sum_i r_i(si) Lsa_i wt_prod_i
  //   dmu/dsa       = Sum_i r_i(si) dLsa_i/dsa wt_prod_i
  //   dmu/dsi       = Sum_i dr_i/dsi Lsa_i wt_prod_i
  // ---------------------------------------------------------------------
  size_t i, deriv_index, cntr = 0, num_deriv_vars = dvv.size();
  RealVector& mean_grad = momentGradsIter->second[0];
  if (mean_grad.length() != num_deriv_vars)
    mean_grad.sizeUninitialized(num_deriv_vars);
  const RealVector2DArray& exp_t1_coeffs = expT1CoeffsIter->second;
  const RealMatrix2DArray& exp_t2_coeffs = expT2CoeffsIter->second;
  const RealMatrix2DArray& exp_t1_coeff_grads = expT1CoeffGradsIter->second;

  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (data_rep->randomVarsKey[deriv_index]) {
      // --------------------------------------------------------------------
      // derivative of All var expansion w.r.t. random var (design insertion)
      // --------------------------------------------------------------------
      if (!expansionCoeffGradFlag) { // required data check
	PCerr << "Error: expansion coefficient gradients not defined in "
	      << "HierarchInterpPolyApproximation::mean_gradient()."
	      << std::endl;
	abort_handler(-1);
      }
      if (data_rep->basisConfigOptions.useDerivs) {
	PCerr << "Error: combination of coefficient gradients and use_"
	      << "derivatives is not supported in HierarchInterpPoly"
	      << "Approximation::mean_gradient()." << std::endl;
	abort_handler(-1);
      }
      mean_grad[i] = expectation_gradient(x, exp_t1_coeff_grads, cntr);
      ++cntr;
    }
    else {
      // ---------------------------------------------------------------------
      // deriv of All var expansion w.r.t. nonrandom var (design augmentation)
      // ---------------------------------------------------------------------
      if (!expansionCoeffFlag) { // required data check
	PCerr << "Error: expansion coefficients not defined in HierarchInterp"
	      << "PolyApproximation::mean_gradient()." << std::endl;
	abort_handler(-1);
      }
      mean_grad[i]
	= expectation_gradient(x, exp_t1_coeffs, exp_t2_coeffs, deriv_index);
    }
  }
  if (all_mode) { compMeanIter->second |=  2; xPrevMeanGrad[key] = x; }
  else compMeanIter->second &= ~2; // deactivate 2-bit: protect mixed use
  return mean_grad;
}


Real HierarchInterpPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool same = (this == hip_approx_2),
    std_mode = data_rep->nonRandomIndices.empty();

  if (same && std_mode && (compVarIter->second & 1))
    return numMomentsIter->second[1];

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  Real covar, mean_1 = mean(), mean_2 = (same) ? mean_1 : hip_approx_2->mean();
  if (speedOverPrecision && product_interpolants())
     // uncentered raw moment available (accept loss of precision)
    covar = expectation(prodT1CoeffsIter->second[poly_approx_2],
			prodT2CoeffsIter->second[poly_approx_2])
          - mean_1 * mean_2;
  else { // form central product interp from scratch and then integrate it
    RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
    central_product_interpolant(hip_approx_2, mean_1, mean_2,
				cov_t1_coeffs, cov_t2_coeffs);
    covar = expectation(cov_t1_coeffs, cov_t2_coeffs);
  }

  if (same && std_mode)
    { numMomentsIter->second[1] = covar; compVarIter->second |= 1; }
  return covar;
}


Real HierarchInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool same = (this == hip_approx_2),
    all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;

  if ( same && all_mode && (compVarIter->second & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar[key]) )
    return numMomentsIter->second[1];

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  Real covar, mean_1 = mean(x),
    mean_2 = (same) ? mean_1 : hip_approx_2->mean(x);
  if (speedOverPrecision && product_interpolants())
    // uncentered raw moment available (accept loss of precision)
    covar = expectation(x, prodT1CoeffsIter->second[poly_approx_2],
			prodT2CoeffsIter->second[poly_approx_2])
          - mean_1 * mean_2;
  else { // form central product interp from scratch and then integrate it
    RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
    central_product_interpolant(hip_approx_2, mean_1, mean_2,
				cov_t1_coeffs, cov_t2_coeffs);
    covar = expectation(x, cov_t1_coeffs, cov_t2_coeffs);
  }

  if (same && all_mode) {
    numMomentsIter->second[1] = covar;
    compVarIter->second |= 1;  xPrevVar[key] = x;
  }
  return covar;
}


Real HierarchInterpPolyApproximation::
combined_covariance(PolynomialApproximation* poly_approx_2)
{
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool  same = (this == hip_approx_2),
    std_mode = data_rep->nonRandomIndices.empty();

  if (same && std_mode && (compVarIter->second & 1))
    return numMomentsIter->second[1];

  Real mean_1 = combined_mean(),
       mean_2 = (same) ? mean_1 : hip_approx_2->combined_mean();

  // Note: an incremental approach is not critical since we mainly use
  //       delta_combined_*() prior to combined_to_active(), and use regular
  //       covariance() after combined_to_active() --> this fn is mainly used
  //       immediately after switch from active to combined stats to provide
  //       reference prior to greedy adaptation.
  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  central_product_interpolant(hsg_driver->combined_hierarchical_variable_sets(),
			      hsg_driver->combined_smolyak_multi_index(),
			      hsg_driver->combined_collocation_key(),
			      combinedExpT1Coeffs, combinedExpT2Coeffs,
			      hip_approx_2->combinedExpT1Coeffs,
			      hip_approx_2->combinedExpT2Coeffs, same,
			      mean_1, mean_2, cov_t1_coeffs, cov_t2_coeffs);

  // evaluate expectation of these central product t1,t2 coefficients
  Real covar =
    expectation(cov_t1_coeffs, cov_t2_coeffs,
		hsg_driver->combined_type1_hierarchical_weight_sets(),
		hsg_driver->combined_type2_hierarchical_weight_sets());

  if (same && std_mode)
    { numMomentsIter->second[1] = covar; compVarIter->second |= 1; }
  return covar;
}


Real HierarchInterpPolyApproximation::
combined_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool  same = (this == hip_approx_2),
    all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if ( same && all_mode && (compVarIter->second & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar[key]) )
    return numMomentsIter->second[1];

  Real mean_1 = combined_mean(x),
       mean_2 = (same) ? mean_1 : hip_approx_2->combined_mean(x);

  // Note: an incremental approach is not critical since we mainly use
  //       delta_combined_*() prior to combined_to_active(), and use regular
  //       covariance() after combined_to_active() --> this fn is mainly used
  //       immediately after switch from active to combined stats to provide
  //       reference prior to greedy adaptation.
  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& comb_sm_mi = hsg_driver->combined_smolyak_multi_index();
  const UShort4DArray&   comb_key = hsg_driver->combined_collocation_key();
  central_product_interpolant(hsg_driver->combined_hierarchical_variable_sets(),
			      comb_sm_mi, comb_key,
			      combinedExpT1Coeffs, combinedExpT2Coeffs,
			      hip_approx_2->combinedExpT1Coeffs,
			      hip_approx_2->combinedExpT2Coeffs, same,
			      mean_1, mean_2, cov_t1_coeffs, cov_t2_coeffs);

  // evaluate expectation of these central product t1,t2 coefficients
  Real covar
    = expectation(x, cov_t1_coeffs, cov_t2_coeffs, comb_sm_mi, comb_key);

  if (same && all_mode) {
    numMomentsIter->second[1] = covar;
    compVarIter->second |= 1;  xPrevVar[key] = x;
  }
  return covar;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design variables are inserted and some
    are augmented) requires no special treatment. */
const RealVector& HierarchInterpPolyApproximation::variance_gradient()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (std_mode && (compVarIter->second & 2))
    return momentGradsIter->second[1];

  // Error check for required data
  if (!expansionCoeffFlag ||
      !expansionCoeffGradFlag) {
    PCerr << "Error: insufficient expansion coefficient data in HierarchInterp"
	  << "PolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  Real mean1 = mean(); const RealVector& mean1_grad = mean_gradient();
  RealMatrix2DArray cov_t1_coeff_grads;
  central_product_gradient_interpolant(this, mean1, mean1, mean1_grad,
				       mean1_grad, cov_t1_coeff_grads);
  RealVector& var_grad = momentGradsIter->second[1];
  var_grad = expectation_gradient(cov_t1_coeff_grads);
  if (std_mode) compVarIter->second |=  2;
  else          compVarIter->second &= ~2;// deactivate 2-bit: protect mixed use
  return var_grad;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionType1CoeffGrads). */
const RealVector& HierarchInterpPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  // if already computed, return previous result
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if ( all_mode && (compVarIter->second & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevVarGrad[key]) )
    // && dvv == dvvPrev)
    return momentGradsIter->second[1];

  // ---------------------------------------------------------------------
  // For xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   mu(sa, si)    = Sum_i r_i(sa, si) wt_prod_i
  //   dmu/ds        = Sum_i dr_i/ds wt_prod_i
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   mu(sa, si)    = Sum_i r_i(si) Lsa_i wt_prod_i
  //   dmu/dsa       = Sum_i r_i(si) dLsa_i/dsa wt_prod_i
  //   dmu/dsi       = Sum_i dr_i/dsi Lsa_i wt_prod_i
  // ---------------------------------------------------------------------
  size_t i, deriv_index, cntr = 0, num_deriv_vars = dvv.size();
  bool insert = false, augment = false;
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    if (data_rep->randomVarsKey[deriv_index]) insert = true;
    else                                      augment = true;
  }

  RealVector2DArray cov_t1c;
  RealMatrix2DArray cov_t1c_grads, cov_t2c;
  Real mean_1 = mean(x);
  if (insert) {
    const RealVector& mean1_grad = mean_gradient(x, dvv);
    central_product_gradient_interpolant(this, mean_1, mean_1, mean1_grad,
					 mean1_grad, cov_t1c_grads);
  }
  if (augment)
    central_product_interpolant(this, mean_1, mean_1, cov_t1c, cov_t2c);

  RealVector& var_grad = momentGradsIter->second[1];
  if (var_grad.length() != num_deriv_vars)
    var_grad.sizeUninitialized(num_deriv_vars);
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    Real& grad_i = var_grad[i];
    if (data_rep->randomVarsKey[deriv_index]) {
      // --------------------------------------------------------------------
      // derivative of All var expansion w.r.t. random var (design insertion)
      // --------------------------------------------------------------------
      if (!expansionCoeffGradFlag) { // required data check
	PCerr << "Error: expansion coefficient gradients not defined in "
	      << "HierarchInterpPolyApproximation::variance_gradient()."
	      << std::endl;
	abort_handler(-1);
      }
      if (data_rep->basisConfigOptions.useDerivs) {
	PCerr << "Error: combination of coefficient gradients and use_"
	      << "derivatives is not supported in HierarchInterpPoly"
	      << "Approximation::variance_gradient()" << std::endl;
	abort_handler(-1);
      }
      var_grad[i] = expectation_gradient(x, cov_t1c_grads, cntr);
      ++cntr;
    }
    else {
      // ---------------------------------------------------------------------
      // deriv of All var expansion w.r.t. nonrandom var (design augmentation)
      // ---------------------------------------------------------------------
      if (!expansionCoeffFlag) { // required data check
	PCerr << "Error: expansion coefficients not defined in Hierarch"
	      << "InterpPolyApproximation::variance_gradient()." << std::endl;
	abort_handler(-1);
      }
      var_grad[i] = expectation_gradient(x, cov_t1c, cov_t2c, deriv_index);
    }
  }
  if (all_mode) { compVarIter->second |=  2; xPrevVarGrad[key] = x; }
  else            compVarIter->second &= ~2;//deactivate 2-bit:protect mixed use
  return var_grad;
}


Real HierarchInterpPolyApproximation::
reference_mean(const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compRefMeanIter->second & 1))
    return refMomentsIter->second[0];

  Real ref_mean
    = expectation(expT1CoeffsIter->second, expT2CoeffsIter->second, ref_key);
  if (std_mode)
    { refMomentsIter->second[0] = ref_mean; compRefMeanIter->second |= 1; }
  return ref_mean;
}


Real HierarchInterpPolyApproximation::
reference_mean(const RealVector& x, const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compRefMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevRefMean[key]))
    return refMomentsIter->second[0];

  Real ref_mean
    = expectation(x, expT1CoeffsIter->second, expT2CoeffsIter->second, ref_key);
  if (all_mode) {
    refMomentsIter->second[0] = ref_mean;
    compRefMeanIter->second |= 1;  xPrevRefMean[key] = x;
  }
  return ref_mean;
}


/** does not require combinedExpT{1,2}Coeffs as works from full maps */
Real HierarchInterpPolyApproximation::
reference_combined_mean(const std::map<UShortArray, UShort2DArray>& ref_key_map)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compRefMeanIter->second & 1))
    return refMomentsIter->second[0];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real ref_mean
    = expectation(expansionType1Coeffs, expansionType2Coeffs,
		  hsg_driver->type1_weight_sets_map(),
		  hsg_driver->type2_weight_sets_map(), ref_key_map);

  if (std_mode)
    { refMomentsIter->second[0] = ref_mean; compRefMeanIter->second |= 1; }
  return ref_mean;
}


/** does not require combinedExpT{1,2}Coeffs as works from full maps */
Real HierarchInterpPolyApproximation::
reference_combined_mean(const RealVector& x,
			const std::map<UShortArray, UShort2DArray>& ref_key_map)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compRefMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevRefMean[key]))
    return refMomentsIter->second[0];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real ref_mean
    = expectation(x, expansionType1Coeffs, expansionType2Coeffs,
		  hsg_driver->smolyak_multi_index_map(),
		  hsg_driver->collocation_key_map(), ref_key_map);
  if (all_mode) {
    refMomentsIter->second[0] = ref_mean;
    compRefMeanIter->second |= 1;  xPrevRefMean[key] = x;
  }
  return ref_mean;
}


Real HierarchInterpPolyApproximation::
reference_variance(const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compRefVarIter->second & 1))
    return refMomentsIter->second[1];

  Real ref_var, ref_mean = reference_mean(ref_key);
  if (speedOverPrecision && product_interpolants())
    // uncentered raw moment available (accept loss of precision)
    ref_var = expectation(prodT1CoeffsIter->second[this],
			  prodT2CoeffsIter->second[this], ref_key)
            - ref_mean * ref_mean;
  else { // form central product interp from scratch and then integrate it
    RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
    central_product_interpolant(this, ref_mean, ref_mean, cov_t1_coeffs,
				cov_t2_coeffs, ref_key);
    ref_var = expectation(cov_t1_coeffs, cov_t2_coeffs, ref_key);
  }

  if (std_mode)
    { refMomentsIter->second[1] = ref_var; compRefVarIter->second |= 1; }
  return ref_var;
}


Real HierarchInterpPolyApproximation::
reference_variance(const RealVector& x, const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compRefVarIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevRefVar[key]))
    return refMomentsIter->second[1];

  Real ref_var, ref_mean = reference_mean(x, ref_key);
  if (speedOverPrecision && product_interpolants())
    // uncentered raw moment available (accept loss of precision)
    ref_var = expectation(x, prodT1CoeffsIter->second[this],
			     prodT2CoeffsIter->second[this], ref_key)
            - ref_mean * ref_mean;
  else { // form central product interp from scratch and then integrate it
    RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
    central_product_interpolant(this, ref_mean, ref_mean, cov_t1_coeffs,
				cov_t2_coeffs, ref_key);
    ref_var = expectation(x, cov_t1_coeffs, cov_t2_coeffs, ref_key);
  }

  if (all_mode) {
    refMomentsIter->second[1] = ref_var;
    compRefVarIter->second |= 1; xPrevRefVar[key] = x;
  }
  return ref_var;
}


Real HierarchInterpPolyApproximation::reference_combined_variance(
  const std::map<UShortArray, UShort2DArray>& ref_key_map)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compRefVarIter->second & 1))
    return refMomentsIter->second[1];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real ref_var, ref_mean = reference_combined_mean(ref_key_map);
  if (speedOverPrecision && product_interpolants())
    // uncentered raw moment available (accept loss of precision)
    ref_var = expectation(productType1Coeffs, productType2Coeffs, this,
			  hsg_driver->type1_weight_sets_map(),
			  hsg_driver->type2_weight_sets_map(), ref_key_map)
            - ref_mean * ref_mean;
  else { // form central product interp from scratch and then integrate it
    std::map<UShortArray, RealVector2DArray> cov_t1c_map;
    std::map<UShortArray, RealMatrix2DArray> cov_t2c_map;
    central_product_interpolant(this, ref_mean, ref_mean, cov_t1c_map,
				cov_t2c_map, ref_key_map);
    ref_var = expectation(cov_t1c_map, cov_t2c_map,
			  hsg_driver->type1_weight_sets_map(),
			  hsg_driver->type2_weight_sets_map(), ref_key_map);
  }

  if (std_mode)
    { refMomentsIter->second[1] = ref_var; compRefVarIter->second |= 1; }
  return ref_var;
}


Real HierarchInterpPolyApproximation::
reference_combined_variance(const RealVector& x,
  const std::map<UShortArray, UShort2DArray>& ref_key_map)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compRefVarIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevRefVar[key]))
    return refMomentsIter->second[1];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real ref_var, ref_mean = reference_combined_mean(x, ref_key_map);
  if (speedOverPrecision && product_interpolants())
    // uncentered raw moment available (accept loss of precision)
    ref_var = expectation(x, productType1Coeffs, productType2Coeffs, this,
			  hsg_driver->smolyak_multi_index_map(),
			  hsg_driver->collocation_key_map(), ref_key_map)
            - ref_mean * ref_mean;
  else { // form central product interp from scratch and then integrate it
    std::map<UShortArray, RealVector2DArray> cov_t1c_map;
    std::map<UShortArray, RealMatrix2DArray> cov_t2c_map;
    central_product_interpolant(this, ref_mean, ref_mean, cov_t1c_map,
				cov_t2c_map, ref_key_map);
    ref_var = expectation(x, cov_t1c_map, cov_t2c_map,
			  hsg_driver->smolyak_multi_index_map(),
			  hsg_driver->collocation_key_map(), ref_key_map);
  }

  if (all_mode) {
    refMomentsIter->second[1] = ref_var;
    compRefVarIter->second |= 1; xPrevRefVar[key] = x;
  }
  return ref_var;
}


Real HierarchInterpPolyApproximation::delta_mean()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compDeltaMeanIter->second & 1))
    return deltaMomentsIter->second[0];

  UShort2DArray incr_key;
  data_rep->hsg_driver()->partition_increment_key(incr_key);

  Real delta_mean
    = expectation(expT1CoeffsIter->second, expT2CoeffsIter->second, incr_key);
  if (std_mode) {
    deltaMomentsIter->second[0] = delta_mean;
    compDeltaMeanIter->second  |= 1;
  }
  return delta_mean;
}


/** This helper avoids recomputing incr_key. */
Real HierarchInterpPolyApproximation::delta_mean(const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compDeltaMeanIter->second & 1))
    return deltaMomentsIter->second[0];

  Real delta_mean
    = expectation(expT1CoeffsIter->second, expT2CoeffsIter->second, incr_key);
  if (std_mode) {
    deltaMomentsIter->second[0] = delta_mean;
    compDeltaMeanIter->second  |= 1;
  }
  return delta_mean;
}


Real HierarchInterpPolyApproximation::delta_mean(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compDeltaMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaMean[key]))
    return deltaMomentsIter->second[0];

  UShort2DArray incr_key;
  data_rep->hsg_driver()->partition_increment_key(incr_key);

  Real delta_mean =
    expectation(x, expT1CoeffsIter->second, expT2CoeffsIter->second, incr_key);
  if (all_mode) {
    deltaMomentsIter->second[0] = delta_mean;
    compDeltaMeanIter->second  |= 1;  xPrevDeltaMean[key] = x;
  }
  return delta_mean;
}


/** This helper avoids recomputing incr_key. */
Real HierarchInterpPolyApproximation::
delta_mean(const RealVector& x, const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compDeltaMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaMean[key]))
    return deltaMomentsIter->second[0];

  Real delta_mean =
    expectation(x, expT1CoeffsIter->second, expT2CoeffsIter->second, incr_key);
  if (all_mode) {
    deltaMomentsIter->second[0] = delta_mean;
    compDeltaMeanIter->second  |= 1;  xPrevDeltaMean[key] = x;
  }
  return delta_mean;
}


Real HierarchInterpPolyApproximation::delta_combined_mean()
{
  // combined increment is the same as the active increment
  //return delta_mean();

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compDeltaMeanIter->second & 1))
    return deltaMomentsIter->second[0];

  // Avoid dependence on metric_roll_up() (combinedExpT{1,2}Coeffs)
  // by employing incr_key_map on expansionType{1,2}Coeffs

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  std::map<UShortArray, UShort2DArray> incr_key_map;
  hsg_driver->partition_increment_key(incr_key_map);

  Real delta_mean
    = expectation(expansionType1Coeffs, expansionType2Coeffs,
		  hsg_driver->type1_weight_sets_map(),
		  hsg_driver->type2_weight_sets_map(), incr_key_map);
  if (std_mode) {
    deltaMomentsIter->second[0] = delta_mean;
    compDeltaMeanIter->second  |= 1;
  }
  return delta_mean;
}


/** This helper avoids recomputing incr_key_map. */
Real HierarchInterpPolyApproximation::
delta_combined_mean(const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  // combined increment is the same as the active increment
  //return delta_mean(incr_key);

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compDeltaMeanIter->second & 1))
    return deltaMomentsIter->second[0];

  // Avoid dependence on metric_roll_up() (combinedExpT{1,2}Coeffs)
  // by employing incr_key_map on expansionType{1,2}Coeffs

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real delta_mean
    = expectation(expansionType1Coeffs, expansionType2Coeffs,
		  hsg_driver->type1_weight_sets_map(),
		  hsg_driver->type2_weight_sets_map(), incr_key_map);
  if (std_mode) {
    deltaMomentsIter->second[0] = delta_mean;
    compDeltaMeanIter->second  |= 1;
  }
  return delta_mean;
}


Real HierarchInterpPolyApproximation::
delta_combined_mean(const RealVector& x)
{
  // Potentially equivalent for all_vars mode as well, but use key maps
  //return delta_mean(x);

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compDeltaMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaMean[key]))
    return deltaMomentsIter->second[0];

  // Avoid dependence on metric_roll_up() (combinedExpT{1,2}Coeffs)
  // by employing incr_key_map on expansionType{1,2}Coeffs

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  std::map<UShortArray, UShort2DArray> incr_key_map;
  hsg_driver->partition_increment_key(incr_key_map);

  Real delta_mean
    = expectation(x, expansionType1Coeffs, expansionType2Coeffs,
		  hsg_driver->smolyak_multi_index_map(),
		  hsg_driver->collocation_key_map(), incr_key_map);
  if (all_mode) {
    deltaMomentsIter->second[0] = delta_mean;
    compDeltaMeanIter->second  |= 1;  xPrevDeltaMean[key] = x;
  }
  return delta_mean;
}


/** This helper avoids recomputing incr_key_map. */
Real HierarchInterpPolyApproximation::
delta_combined_mean(const RealVector& x,
		    const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  // TO DO: equivalence for all_vars mode ?
  //return delta_mean(x, incr_key);

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compDeltaMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaMean[key]))
    return deltaMomentsIter->second[0];

  // Avoid dependence on metric_roll_up() (combinedExpT{1,2}Coeffs)
  // by employing incr_key_map on expansionType{1,2}Coeffs

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real delta_mean
    = expectation(x, expansionType1Coeffs, expansionType2Coeffs,
		  hsg_driver->smolyak_multi_index_map(),
		  hsg_driver->collocation_key_map(), incr_key_map);
  if (all_mode) {
    deltaMomentsIter->second[0] = delta_mean;
    compDeltaMeanIter->second  |= 1;  xPrevDeltaMean[key] = x;
  }
  return delta_mean;
}


/** This helper avoids recomputing ref_key and incr_key (and also
    eliminates poly_approx_2 from delta_covariance()). */
Real HierarchInterpPolyApproximation::
delta_variance(const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compDeltaVarIter->second & 1))
    return deltaMomentsIter->second[1];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real delta_var;
  if (product_interpolants())
    delta_var = delta_covariance(expT1CoeffsIter->second,
      expT2CoeffsIter->second, expT1CoeffsIter->second, expT2CoeffsIter->second,
      true, prodT1CoeffsIter->second[this], prodT2CoeffsIter->second[this],
      hsg_driver->type1_hierarchical_weight_sets(),
      hsg_driver->type2_hierarchical_weight_sets(), ref_key, incr_key);
  else { // rebuild product interpolant from scratch
    RealVector2DArray prod_t1c; RealMatrix2DArray prod_t2c;
    product_interpolant(this, prod_t1c, prod_t2c);
    delta_var = delta_covariance(expT1CoeffsIter->second,
      expT2CoeffsIter->second, expT1CoeffsIter->second, expT2CoeffsIter->second,
      true, prod_t1c, prod_t2c, hsg_driver->type1_hierarchical_weight_sets(),
      hsg_driver->type2_hierarchical_weight_sets(), ref_key, incr_key);
  }

  if (std_mode)
    { deltaMomentsIter->second[1] = delta_var; compDeltaVarIter->second |= 1; }
  return delta_var;
}


/** This helper avoids recomputing ref_key and incr_key (and also
    eliminates poly_approx_2 from delta_covariance()). */
Real HierarchInterpPolyApproximation::
delta_variance(const RealVector& x, const UShort2DArray& ref_key,
	       const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compDeltaVarIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaVar[key]))
    return deltaMomentsIter->second[1];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real delta_var;
  if (product_interpolants())
    delta_var = delta_covariance(x, expT1CoeffsIter->second,
      expT2CoeffsIter->second, expT1CoeffsIter->second, expT2CoeffsIter->second,
      true, prodT1CoeffsIter->second[this], prodT2CoeffsIter->second[this],
      hsg_driver->smolyak_multi_index(), hsg_driver->collocation_key(),
      ref_key, incr_key);
  else { // rebuild product interpolant from scratch
    RealVector2DArray prod_t1c; RealMatrix2DArray prod_t2c;
    product_interpolant(this, prod_t1c, prod_t2c);
    delta_var = delta_covariance(x, expT1CoeffsIter->second,
      expT2CoeffsIter->second, expT1CoeffsIter->second, expT2CoeffsIter->second,
      true, prod_t1c, prod_t2c, hsg_driver->smolyak_multi_index(),
      hsg_driver->collocation_key(), ref_key, incr_key);
  }
  
  if (all_mode) {
    deltaMomentsIter->second[1] = delta_var;
    compDeltaVarIter->second   |= 1;  xPrevDeltaVar[key] = x;
  }
  return delta_var;
}


/** This helper avoids recomputing ref_key_map and incr_key_map (and
    also eliminates poly_approx_2 from delta_combined_covariance()). */
Real HierarchInterpPolyApproximation::
delta_combined_variance(
  const std::map<UShortArray, UShort2DArray>& ref_key_map,
  const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compDeltaVarIter->second & 1))
    return deltaMomentsIter->second[1];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real delta_var;
  if (product_interpolants())
    delta_var = delta_covariance(expansionType1Coeffs, expansionType2Coeffs,
      expansionType1Coeffs, expansionType2Coeffs, true,
      prodT1CoeffsIter->second[this], prodT2CoeffsIter->second[this],
      hsg_driver->type1_weight_sets_map(), hsg_driver->type2_weight_sets_map(),
      data_rep->activeKey, ref_key_map, incr_key_map);
  else {
    // rebuild product interpolant from scratch for all active sets (only need
    // expectation of R1R2 on active incr_key, but depends on active ref_key
    // since product interp is constructed bottom-up from hierarch surpluses)
    RealVector2DArray prod_t1c; RealMatrix2DArray prod_t2c;
    product_interpolant(this, prod_t1c, prod_t2c);
    delta_var = delta_covariance(expansionType1Coeffs, expansionType2Coeffs,
      expansionType1Coeffs, expansionType2Coeffs, true, prod_t1c, prod_t2c,
      hsg_driver->type1_weight_sets_map(), hsg_driver->type2_weight_sets_map(),
      data_rep->activeKey, ref_key_map, incr_key_map);
  }

  if (std_mode)
    { deltaMomentsIter->second[1] = delta_var; compDeltaVarIter->second |= 1; }
  return delta_var;
}


/** This helper avoids recomputing ref_key_map and incr_key_map (and
    also eliminates poly_approx_2 from delta_combined_covariance()). */
Real HierarchInterpPolyApproximation::
delta_combined_variance(const RealVector& x,
  const std::map<UShortArray, UShort2DArray>& ref_key_map,
  const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compDeltaVarIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaVar[key]))
    return deltaMomentsIter->second[1];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  Real delta_var;
  if (product_interpolants())
    delta_var = delta_covariance(x, expansionType1Coeffs, expansionType2Coeffs,
      expansionType1Coeffs, expansionType2Coeffs, true,
      prodT1CoeffsIter->second[this], prodT2CoeffsIter->second[this],
      hsg_driver->smolyak_multi_index_map(), hsg_driver->collocation_key_map(),
      data_rep->activeKey, ref_key_map, incr_key_map);
  else {
    // rebuild product interpolant from scratch for all active sets (only need
    // expectation of R1R2 on active incr_key, but depends on active ref_key
    // since product interp is constructed bottom-up from hierarch surpluses)
    RealVector2DArray prod_t1c; RealMatrix2DArray prod_t2c;
    product_interpolant(this, prod_t1c, prod_t2c);
    delta_var = delta_covariance(x, expansionType1Coeffs, expansionType2Coeffs,
      expansionType1Coeffs, expansionType2Coeffs, true, prod_t1c, prod_t2c,
      hsg_driver->smolyak_multi_index_map(), hsg_driver->collocation_key_map(),
      data_rep->activeKey, ref_key_map, incr_key_map);
  }

  if (all_mode) {
    deltaMomentsIter->second[1] = delta_var;
    compDeltaVarIter->second   |= 1;  xPrevDeltaVar[key] = x;
  }
  return delta_var;
}


Real HierarchInterpPolyApproximation::
delta_std_deviation(const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  // Preserve precision by avoiding subtractive cancellation
  // delta-sigma = sqrt( var0 + delta-var ) - sigma0
  //             = [ sqrt(1 + delta_var / var0) - 1 ] * sigma0
  //             = sqrt1pm1(delta_var / var0) * sigma0
  // where sqrt1pm1(x) = sqrt(1+x) - 1 = expm1[ log1p(x) / 2 ]

  Real delta_var = delta_variance(ref_key, incr_key),
       var0      = reference_variance(ref_key),
       sigma0    = (var0 > 0.) ? std::sqrt(var0) : 0.;

  // usual case: var0 and sigma0 are both positive, prefer sqrt1pm1
  // unless 1+x could go negative
  if ( sigma0 > 0. &&
       ( delta_var >= 0. || std::abs(delta_var) < var0 / 2.) )
    return bmth::sqrt1pm1(delta_var / var0) * sigma0;
  // negative var0 only supported if var1 recovers to positive
  else {
    Real var1 = var0 + delta_var;
    return (var1 > 0.) ? std::sqrt(var1) - sigma0 : 0.;
  }
}


Real HierarchInterpPolyApproximation::
delta_std_deviation(const RealVector& x, const UShort2DArray& ref_key,
		    const UShort2DArray& incr_key)
{
  // delta-sigma = sqrt( var0 + delta-var ) - sigma0
  //             = [ sqrt(1 + delta_var / var0) - 1 ] * sigma0
  //             = sqrt1pm1(delta_var / var0) * sigma0
  // where sqrt1pm1(x) = expm1[ log1p(x) / 2 ]

  Real delta_var = delta_variance(x, ref_key, incr_key),
       var0      = reference_variance(x, ref_key),
       sigma0    = (var0 > 0.) ? std::sqrt(var0) : 0.;

  // usual case: var0 and sigma0 are both positive, prefer sqrt1pm1
  // unless 1+x could go negative
  if ( sigma0 > 0. &&
       ( delta_var >= 0. || std::abs(delta_var) < var0 / 2.) )
    return bmth::sqrt1pm1(delta_var / var0) * sigma0;
  // negative var0 only supported if var1 recovers to positive
  else {
    Real var1 = var0 + delta_var;
    return (var1 > 0.) ? std::sqrt(var1) - sigma0 : 0.;
  }
}


Real HierarchInterpPolyApproximation::
delta_combined_std_deviation(
  const std::map<UShortArray, UShort2DArray>& ref_key_map,
  const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  // Preserve precision by avoiding subtractive cancellation
  // delta-sigma = sqrt( var0 + delta-var ) - sigma0
  //             = [ sqrt(1 + delta_var / var0) - 1 ] * sigma0
  //             = sqrt1pm1(delta_var / var0) * sigma0
  // where sqrt1pm1(x) = sqrt(1+x) - 1 = expm1[ log1p(x) / 2 ]

  Real delta_var = delta_combined_variance(ref_key_map, incr_key_map),
       var0      = reference_combined_variance(ref_key_map),
       sigma0    = (var0 > 0.) ? std::sqrt(var0) : 0.;

  // usual case: var0 and sigma0 are both positive, prefer sqrt1pm1
  // unless 1+x could go negative
  if ( sigma0 > 0. &&
       ( delta_var >= 0. || std::abs(delta_var) < var0 / 2.) )
    return bmth::sqrt1pm1(delta_var / var0) * sigma0;
  // negative var0 only supported if var1 recovers to positive
  else {
    Real var1 = var0 + delta_var;
    return (var1 > 0.) ? std::sqrt(var1) - sigma0 : 0.;
  }
}


Real HierarchInterpPolyApproximation::
delta_combined_std_deviation(const RealVector& x,
  const std::map<UShortArray, UShort2DArray>& ref_key_map,
  const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  // delta-sigma = sqrt( var0 + delta-var ) - sigma0
  //             = [ sqrt(1 + delta_var / var0) - 1 ] * sigma0
  //             = sqrt1pm1(delta_var / var0) * sigma0
  // where sqrt1pm1(x) = expm1[ log1p(x) / 2 ]

  Real delta_var = delta_combined_variance(x, ref_key_map, incr_key_map),
       var0      = reference_combined_variance(x, ref_key_map),
       sigma0    = (var0 > 0.) ? std::sqrt(var0) : 0.;

  // usual case: var0 and sigma0 are both positive, prefer sqrt1pm1
  // unless 1+x could go negative
  if ( sigma0 > 0. &&
       ( delta_var >= 0. || std::abs(delta_var) < var0 / 2.) )
    return bmth::sqrt1pm1(delta_var / var0) * sigma0;
  // negative var0 only supported if var1 recovers to positive
  else {
    Real var1 = var0 + delta_var;
    return (var1 > 0.) ? std::sqrt(var1) - sigma0 : 0.;
  }
}


Real HierarchInterpPolyApproximation::
delta_beta(bool cdf_flag, Real z_bar, const UShort2DArray& ref_key,
	   const UShort2DArray& incr_key)
{
  Real mu0      = reference_mean(ref_key),
    delta_mu    = delta_mean(incr_key),
    var0        = reference_variance(ref_key),
    delta_sigma = delta_std_deviation(ref_key, incr_key);
  return delta_beta_map(mu0, delta_mu, var0, delta_sigma, cdf_flag, z_bar);
}


Real HierarchInterpPolyApproximation::
delta_beta(const RealVector& x, bool cdf_flag, Real z_bar,
	   const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  Real mu0      = reference_mean(x, ref_key),
    delta_mu    = delta_mean(x, incr_key),
    var0        = reference_variance(x, ref_key),
    delta_sigma = delta_std_deviation(x, ref_key, incr_key);
  return delta_beta_map(mu0, delta_mu, var0, delta_sigma, cdf_flag, z_bar);
}


Real HierarchInterpPolyApproximation::
delta_combined_beta(bool cdf_flag, Real z_bar,
		    const std::map<UShortArray, UShort2DArray>& ref_key_map,
		    const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  Real mu0      = reference_combined_mean(ref_key_map),
    delta_mu    = delta_combined_mean(incr_key_map),
    var0        = reference_combined_variance(ref_key_map),
    delta_sigma = delta_combined_std_deviation(ref_key_map, incr_key_map);
  return delta_beta_map(mu0, delta_mu, var0, delta_sigma, cdf_flag, z_bar);
}


Real HierarchInterpPolyApproximation::
delta_combined_beta(const RealVector& x, bool cdf_flag, Real z_bar,
		    const std::map<UShortArray, UShort2DArray>& ref_key_map,
		    const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  Real mu0      = reference_combined_mean(x, ref_key_map),
    delta_mu    = delta_combined_mean(x, incr_key_map),
    var0        = reference_combined_variance(x, ref_key_map),
    delta_sigma = delta_combined_std_deviation(x, ref_key_map, incr_key_map);
  return delta_beta_map(mu0, delta_mu, var0, delta_sigma, cdf_flag, z_bar);
}


Real HierarchInterpPolyApproximation::
delta_beta_map(Real mu0, Real delta_mu, Real var0, Real delta_sigma,
	       bool cdf_flag, Real z_bar)
{
  //  CDF delta-beta = (mu1 - z-bar)/sigma1 - (mu0 - z-bar)/sigma0
  //    = (mu1 sigma0 - z-bar sigma0 - mu0 sigma1 + z-bar sigma1)/sigma1/sigma0
  //    = (delta-mu sigma0 - mu0 delta-sigma + z-bar delta-sigma)/sigma1/sigma0
  //    = (delta-mu - delta-sigma beta0)/sigma1
  // CCDF delta-beta = (z-bar - mu1)/sigma1 - (z-bar - mu0)/sigma0
  //    = (z-bar sigma0 - mu1 sigma0 - z-bar sigma1 + mu0 sigma1)/sigma1/sigma0
  //    = (mu0 delta-sigma - delta-mu sigma0 - z-bar delta-sigma)/sigma1/sigma0
  //    = -delta-mu/sigma1 - delta_sigma (z-bar - mu0) / sigma0 / sigma1
  //    = (-delta-mu - delta-sigma beta0)/sigma1

  // Error traps are needed for zero variance: a single point ref grid
  // (level=0 sparse or m=1 tensor) has zero variance.  Unchanged response
  // values along an index set could then cause sigma1 also = 0.
  Real beta0, sigma0 = (var0 > 0.) ? std::sqrt(var0) : 0.,
    sigma1  = sigma0 + delta_sigma;
  if (cdf_flag) {
    if (sigma0 > SMALL_NUMBER && sigma1 > SMALL_NUMBER) {
      beta0 = (mu0 - z_bar) / sigma0;
      return ( delta_mu - delta_sigma * beta0) / sigma1;
    }
    else if (sigma1 > SMALL_NUMBER)// neglect beta0 term (zero init reliability)
      return delta_mu / sigma1; // or delta = beta1 = (mu1 - z_bar) / sigma1 ?
    else if (sigma0 > SMALL_NUMBER)// assume beta1 = 0 -> delta = -beta0
      return (z_bar - mu0) / sigma0;
    else                           // assume beta0 = beta1 = 0
      return 0;
  }
  else {
    if (sigma0 > SMALL_NUMBER && sigma1 > SMALL_NUMBER) {
      beta0 = (z_bar - mu0) / sigma0;
      return (-delta_mu - delta_sigma * beta0) / sigma1;
    }
    else if (sigma1 > SMALL_NUMBER)// neglect beta0 term (zero init reliability)
      return -delta_mu / sigma1;
    else if (sigma0 > SMALL_NUMBER)// assume beta1 = 0 -> delta = -beta0
      return (mu0 - z_bar) / sigma0;
    else                           // assume beta0 = beta1 = 0
      return 0;
  }
}


Real HierarchInterpPolyApproximation::
delta_z(bool cdf_flag, Real beta_bar, const UShort2DArray& ref_key,
	const UShort2DArray& incr_key)
{
  //  CDF delta-z = (mu1 - sigma1 beta-bar) - (mu0 - sigma0 beta-bar)
  //              = delta-mu - delta-sigma * beta-bar
  // CCDF delta-z = (mu1 + sigma1 beta-bar) - (mu0 + sigma0 beta-bar)
  //              = delta-mu + delta-sigma * beta-bar

  Real delta_mu = delta_mean(incr_key),
    delta_sigma = delta_std_deviation(ref_key, incr_key);
  return (cdf_flag) ? delta_mu - delta_sigma * beta_bar :
                      delta_mu + delta_sigma * beta_bar;
}


Real HierarchInterpPolyApproximation::
delta_z(const RealVector& x, bool cdf_flag, Real beta_bar,
	const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  //  CDF delta-z = (mu1 - sigma1 beta-bar) - (mu0 - sigma0 beta-bar)
  //              = delta-mu - delta-sigma * beta-bar
  // CCDF delta-z = (mu1 + sigma1 beta-bar) - (mu0 + sigma0 beta-bar)
  //              = delta-mu + delta-sigma * beta-bar

  Real delta_mu = delta_mean(x, incr_key),
    delta_sigma = delta_std_deviation(x, ref_key, incr_key);
  return (cdf_flag) ? delta_mu - delta_sigma * beta_bar :
                      delta_mu + delta_sigma * beta_bar;
}


Real HierarchInterpPolyApproximation::
delta_combined_z(bool cdf_flag, Real beta_bar,
		 const std::map<UShortArray, UShort2DArray>& ref_key_map,
		 const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  //  CDF delta-z = (mu1 - sigma1 beta-bar) - (mu0 - sigma0 beta-bar)
  //              = delta-mu - delta-sigma * beta-bar
  // CCDF delta-z = (mu1 + sigma1 beta-bar) - (mu0 + sigma0 beta-bar)
  //              = delta-mu + delta-sigma * beta-bar

  Real delta_mu = delta_combined_mean(incr_key_map),
    delta_sigma = delta_combined_std_deviation(ref_key_map, incr_key_map);
  return (cdf_flag) ? delta_mu - delta_sigma * beta_bar :
                      delta_mu + delta_sigma * beta_bar;
}


Real HierarchInterpPolyApproximation::
delta_combined_z(const RealVector& x, bool cdf_flag, Real beta_bar,
		 const std::map<UShortArray, UShort2DArray>& ref_key_map,
		 const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  //  CDF delta-z = (mu1 - sigma1 beta-bar) - (mu0 - sigma0 beta-bar)
  //              = delta-mu - delta-sigma * beta-bar
  // CCDF delta-z = (mu1 + sigma1 beta-bar) - (mu0 + sigma0 beta-bar)
  //              = delta-mu + delta-sigma * beta-bar

  Real delta_mu = delta_combined_mean(x, incr_key_map),
    delta_sigma = delta_combined_std_deviation(x, ref_key_map, incr_key_map);
  return (cdf_flag) ? delta_mu - delta_sigma * beta_bar :
                      delta_mu + delta_sigma * beta_bar;
}


Real HierarchInterpPolyApproximation::
delta_covariance(PolynomialApproximation* poly_approx_2)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;

  // Supports multiple grid increments in discerning nominal from delta based
  // on isotropic/anisotropic/generalized index set increments.  In current
  // use, 2D keys with set ranges are sufficient: level -> {start,end} set.
  // In the future, may need 3D keys for level/set/point.
  bool  same = (this == hip_approx_2),
    std_mode = data_rep->nonRandomIndices.empty();
  if (same && std_mode && (compDeltaVarIter->second & 1))
    return deltaMomentsIter->second[1];

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::delta_covariance()" << std::endl;
    abort_handler(-1);
  }

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  UShort2DArray ref_key, incr_key;
  hsg_driver->partition_keys(ref_key, incr_key);

  Real delta_covar;
  if (product_interpolants())
    delta_covar = delta_covariance(expT1CoeffsIter->second,
      expT2CoeffsIter->second, hip_approx_2->expT1CoeffsIter->second,
      hip_approx_2->expT2CoeffsIter->second, same,
      prodT1CoeffsIter->second[poly_approx_2],
      prodT2CoeffsIter->second[poly_approx_2],
      hsg_driver->type1_hierarchical_weight_sets(),
      hsg_driver->type2_hierarchical_weight_sets(), ref_key, incr_key);
  else {
    // rebuild product interpolant from scratch for all sets (only need
    // expectation of R1R2 on incr_key, but depends on ref_key since product
    // interpolant must be constructed bottom-up from hierarchical surpluses)
    RealVector2DArray prod_t1c; RealMatrix2DArray prod_t2c;
    product_interpolant(hip_approx_2, prod_t1c, prod_t2c);
    delta_covar = delta_covariance(expT1CoeffsIter->second,
      expT2CoeffsIter->second, hip_approx_2->expT1CoeffsIter->second,
      hip_approx_2->expT2CoeffsIter->second, same, prod_t1c, prod_t2c,
      hsg_driver->type1_hierarchical_weight_sets(),
      hsg_driver->type2_hierarchical_weight_sets(), ref_key, incr_key);
  }

  if (same && std_mode) {
    deltaMomentsIter->second[1] = delta_covar;
    compDeltaVarIter->second |= 1;
  }
  return delta_covar;
}


Real HierarchInterpPolyApproximation::
delta_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;

  bool  same = (this == hip_approx_2),
    all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (same && all_mode && (compDeltaVarIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaVar[key]))
    return deltaMomentsIter->second[1];

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::delta_covariance()" << std::endl;
    abort_handler(-1);
  }

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  UShort2DArray ref_key, incr_key;
  hsg_driver->partition_keys(ref_key, incr_key);

  Real delta_covar;
  if (product_interpolants())
    delta_covar = delta_covariance(x, expT1CoeffsIter->second,
      expT2CoeffsIter->second, hip_approx_2->expT1CoeffsIter->second,
      hip_approx_2->expT2CoeffsIter->second, same,
      prodT1CoeffsIter->second[poly_approx_2],
      prodT2CoeffsIter->second[poly_approx_2],
      hsg_driver->smolyak_multi_index(), hsg_driver->collocation_key(),
      ref_key, incr_key);
  else {
    // rebuild product interpolant from scratch for all sets (only need
    // expectation of R1R2 on incr_key, but depends on ref_key since product
    // interpolant must be constructed bottom-up from hierarchical surpluses)
    RealVector2DArray prod_t1c; RealMatrix2DArray prod_t2c;
    product_interpolant(hip_approx_2, prod_t1c, prod_t2c);
    delta_covar = delta_covariance(x, expT1CoeffsIter->second,
      expT2CoeffsIter->second, hip_approx_2->expT1CoeffsIter->second,
      hip_approx_2->expT2CoeffsIter->second, same, prod_t1c, prod_t2c,
      hsg_driver->smolyak_multi_index(), hsg_driver->collocation_key(),
      ref_key, incr_key);
  }

  if (same && all_mode) {
    deltaMomentsIter->second[1] = delta_covar;
    compDeltaVarIter->second |= 1; xPrevDeltaVar[key] = x;
  }
  return delta_covar;
}


Real HierarchInterpPolyApproximation::
delta_combined_covariance(PolynomialApproximation* poly_approx_2)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  bool  same = (this == hip_approx_2),
    std_mode = data_rep->nonRandomIndices.empty();
  if (same && std_mode && (compDeltaVarIter->second & 1))
    return deltaMomentsIter->second[1];
  
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  std::map<UShortArray, UShort2DArray> ref_key_map, incr_key_map;
  hsg_driver->partition_keys(ref_key_map, incr_key_map);
  
  // For combined statistics, we utilize the full coefficient maps
  // (expansionType{1,2}Coeffs) for the expected values using ref_key_map,
  // but only require the active coefficients (r1r2_t{1,2}_coeffs and active
  // expansionType{1,2}Coeffs) for the deltas using the active incr_key.
  Real delta_covar;
  if (product_interpolants())
    delta_covar = delta_covariance(expansionType1Coeffs, expansionType2Coeffs,
      hip_approx_2->expansionType1Coeffs, hip_approx_2->expansionType2Coeffs,
      same, prodT1CoeffsIter->second[poly_approx_2],
      prodT2CoeffsIter->second[poly_approx_2],
      hsg_driver->type1_weight_sets_map(), hsg_driver->type2_weight_sets_map(),
      data_rep->activeKey, ref_key_map, incr_key_map);
  else {
    // rebuild product interpolant from scratch for all active sets (only need
    // expectation of R1R2 on active incr_key, but depends on active ref_key
    // since product interp is constructed bottom-up from hierarch surpluses)
    RealVector2DArray prod_t1c; RealMatrix2DArray prod_t2c;
    product_interpolant(hip_approx_2, prod_t1c, prod_t2c);
    delta_covar = delta_covariance(expansionType1Coeffs, expansionType2Coeffs,
      hip_approx_2->expansionType1Coeffs, hip_approx_2->expansionType2Coeffs,
      same, prod_t1c, prod_t2c, hsg_driver->type1_weight_sets_map(),
      hsg_driver->type2_weight_sets_map(), data_rep->activeKey,
      ref_key_map, incr_key_map);
  }

  if (same && std_mode) {
    deltaMomentsIter->second[1] = delta_covar;
    compDeltaVarIter->second |= 1;
  }
  return delta_covar;
}


Real HierarchInterpPolyApproximation::
delta_combined_covariance(const RealVector& x,
			  PolynomialApproximation* poly_approx_2)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;

  bool  same = (this == hip_approx_2),
    all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (same && all_mode && (compDeltaVarIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaVar[key]))
    return deltaMomentsIter->second[1];

  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  std::map<UShortArray, UShort2DArray> ref_key_map, incr_key_map;
  hsg_driver->partition_keys(ref_key_map, incr_key_map);

  Real delta_covar;
  if (product_interpolants())
    delta_covar = delta_covariance(x, expansionType1Coeffs,
      expansionType2Coeffs, hip_approx_2->expansionType1Coeffs,
      hip_approx_2->expansionType2Coeffs, same,
      prodT1CoeffsIter->second[poly_approx_2],
      prodT2CoeffsIter->second[poly_approx_2],
      hsg_driver->smolyak_multi_index_map(), hsg_driver->collocation_key_map(),
      data_rep->activeKey, ref_key_map, incr_key_map);
  else {
    // rebuild product interpolant from scratch for all active sets (only need
    // expectation of R1R2 on active incr_key, but depends on active ref_key
    // since product interp is constructed bottom-up from hierarch surpluses)
    RealVector2DArray prod_t1c; RealMatrix2DArray prod_t2c;
    product_interpolant(hip_approx_2, prod_t1c, prod_t2c);
    delta_covar = delta_covariance(x, expansionType1Coeffs,
      expansionType2Coeffs, hip_approx_2->expansionType1Coeffs,
      hip_approx_2->expansionType2Coeffs, same, prod_t1c, prod_t2c,
      hsg_driver->smolyak_multi_index_map(), hsg_driver->collocation_key_map(),
      data_rep->activeKey, ref_key_map, incr_key_map);
  }

  if (same && all_mode) {
    deltaMomentsIter->second[1] = delta_covar;
    compDeltaVarIter->second |= 1; xPrevDeltaVar[key] = x;
  }
  return delta_covar;
}


Real HierarchInterpPolyApproximation::
delta_covariance(const RealVector2DArray& r1_t1_coeffs,
		 const RealMatrix2DArray& r1_t2_coeffs,
		 const RealVector2DArray& r2_t1_coeffs,
		 const RealMatrix2DArray& r2_t2_coeffs, bool same,
		 const RealVector2DArray& r1r2_t1_coeffs,
		 const RealMatrix2DArray& r1r2_t2_coeffs,
		 const RealVector2DArray& t1_wts,
		 const RealMatrix2DArray& t2_wts, const UShort2DArray& ref_key,
		 const UShort2DArray& incr_key)
{
  // Compute surplus for r1, r2, and r1r2 and retrieve reference mean values
  Real  ref_mean_r1 =
    expectation(r1_t1_coeffs, r1_t2_coeffs, t1_wts, t2_wts, ref_key),
      delta_mean_r1 =
    expectation(r1_t1_coeffs, r1_t2_coeffs, t1_wts, t2_wts, incr_key),
        ref_mean_r2 = (same) ? ref_mean_r1 :
    expectation(r2_t1_coeffs, r2_t2_coeffs, t1_wts, t2_wts, ref_key),
      delta_mean_r2 = (same) ? delta_mean_r1 :
    expectation(r2_t1_coeffs, r2_t2_coeffs, t1_wts, t2_wts, incr_key),
    delta_mean_r1r2 =
    expectation(r1r2_t1_coeffs, r1r2_t2_coeffs, t1_wts, t2_wts, incr_key);

  // Hierarchical increment to covariance:
  // \Delta\Sigma_ij = \Sigma^1_ij - \Sigma^0_ij
  //   = ( E[Ri Rj]^1 - \mu_i^1 \mu_j^1 ) - ( E[Ri Rj]^0 - \mu_i^0 \mu_j^0 )
  //   = E[Ri Rj]^0 + \DeltaE[Ri Rj] - E[Ri Rj]^0 
  //     - (\mu_i^0 + \Delta\mu_i) (\mu_j^0 + \Delta\mu_j) + \mu_i^0 \mu_j^0
  //   = \DeltaE[Ri Rj] - \Delta\mu_i \mu_j^0 - \mu_i^0 \Delta\mu_j
  //     - \Delta\mu_i \Delta\mu_j
  // Note: \DeltaE[Ri Rj] can be expanded to E[\DeltaR_i R_j^0]
  //       + E[\DeltaR_j R_i^0] + E[\DeltaR_i \DeltaR_j].  While this would
  //       reduce subtractive cancellation from (dominant) term E[R_i^0 R_j^0],
  //       it would incur additional product bookkeeping.
#ifdef DEBUG
  PCout << "\n  delta_mean_r1r2 = " << delta_mean_r1r2
	<< "\n  ref_mean_r1 = " << ref_mean_r1
	<< "  delta_mean_r1 = " << delta_mean_r1
	<< "\n  ref_mean_r2 = " << ref_mean_r2
	<< "  delta_mean_r2 = " << delta_mean_r2 << std::endl;
#endif // DEBUG
  return delta_mean_r1r2 - ref_mean_r1 * delta_mean_r2
    - ref_mean_r2 * delta_mean_r1 - delta_mean_r1 * delta_mean_r2;
}


Real HierarchInterpPolyApproximation::
delta_covariance(const std::map<UShortArray, RealVector2DArray>& r1_t1c_map,
		 const std::map<UShortArray, RealMatrix2DArray>& r1_t2c_map,
		 const std::map<UShortArray, RealVector2DArray>& r2_t1c_map,
		 const std::map<UShortArray, RealMatrix2DArray>& r2_t2c_map,
		 bool same, const RealVector2DArray& r1r2_t1c,
		 const RealMatrix2DArray& r1r2_t2c,
		 const std::map<UShortArray, RealVector2DArray>& t1_wts_map,
		 const std::map<UShortArray, RealMatrix2DArray>& t2_wts_map,
		 const UShortArray& active_key,
		 const std::map<UShortArray, UShort2DArray>& ref_key_map,
		 const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  // Compute surplus for r1, r2, and r1r2 and retrieve reference mean values
  std::map<UShortArray, RealVector2DArray>::const_iterator r1t1_cit
    = r1_t1c_map.find(active_key), t1w_cit = t1_wts_map.find(active_key);
  std::map<UShortArray, RealMatrix2DArray>::const_iterator r1t2_cit
    = r1_t2c_map.find(active_key), t2w_cit = t2_wts_map.find(active_key);
  std::map<UShortArray, UShort2DArray>::const_iterator
    //ref_cit = ref_key_map.find(active_key),
    incr_cit = incr_key_map.find(active_key);
  if (r1t1_cit == r1_t1c_map.end() || t1w_cit == t1_wts_map.end() ||
      r1t2_cit == r1_t2c_map.end() || t2w_cit == t2_wts_map.end() ||
      incr_cit == incr_key_map.end()) { //|| ref_cit == ref_key_map.end()) {
    PCerr << "Error: failure in key lookup in HierarchInterpPolyApproximation"
	  << "::delta_covariance()" << std::endl;
    abort_handler(-1);
  }
  const RealVector2DArray& active_t1_wts   =  t1w_cit->second;
  const RealMatrix2DArray& active_t2_wts   =  t2w_cit->second;
  const UShort2DArray&     active_incr_key = incr_cit->second;
  //const UShort2DArray&   active_ref_key  =  ref_cit->second;

  // increments could in general be defined by incr_key_map, but currently
  // are limited to the active incr_key (enabling short cuts below)
  Real  ref_mean_r1 =
    expectation(r1_t1c_map, r1_t2c_map, t1_wts_map, t2_wts_map, ref_key_map),
    //  ref_mean_r1_active =
    //expectation(r1t1_cit->second, r1t2_cit->second, active_t1_wts,
    //            active_t2_wts, active_ref_key),
    //  delta_mean_r1_all =
    //expectation(r1_t1c_map, r1_t2c_map, t1_wts_map, t2_wts_map, incr_key_map),
      delta_mean_r1 =
    expectation(r1t1_cit->second, r1t2_cit->second, active_t1_wts,
		active_t2_wts, active_incr_key); // shortcut
  Real ref_mean_r2, delta_mean_r2;
  if (same)
    { ref_mean_r2 = ref_mean_r1; delta_mean_r2 = delta_mean_r1; }
  else {
    std::map<UShortArray, RealVector2DArray>::const_iterator r2t1_cit
      = r2_t1c_map.find(active_key);
    std::map<UShortArray, RealMatrix2DArray>::const_iterator r2t2_cit
      = r2_t2c_map.find(active_key);
    ref_mean_r2 =
      expectation(r2_t1c_map, r2_t2c_map, t1_wts_map, t2_wts_map, ref_key_map);
    //expectation(r2t1_cit->second, r2t2_cit->second, active_t1_wts,
    //		  active_t2_wts, active_ref_key);
    delta_mean_r2 =
    //expectation(r2_t1c_map, r2_t2c_map, t1_wts_map, t2_wts_map, incr_key_map),
      expectation(r2t1_cit->second, r2t2_cit->second, active_t1_wts,
		  active_t2_wts, active_incr_key); // shortcut
  }
  Real delta_mean_r1r2 =
    expectation(r1r2_t1c, r1r2_t2c, active_t1_wts, active_t2_wts,
		active_incr_key); // shortcut

  // same expression as standard expansion mode case above
#ifdef DEBUG
  PCout << "\n  delta_mean_r1r2 (active) = " << delta_mean_r1r2
        << "\n  ref_mean_r1 (all) = "        << ref_mean_r1
      //<<   "  ref_mean_r1 (active) = "     << ref_mean_r1_active
      //<< "  delta_mean_r1 (all) = "        << delta_mean_r1_all
	<< "  delta_mean_r1 (active) = "     << delta_mean_r1
	<< "\n  ref_mean_r2 (all) = "        << ref_mean_r2
	<< "  delta_mean_r2 (active) = "     << delta_mean_r2 << std::endl;
#endif // DEBUG
  return delta_mean_r1r2 - ref_mean_r1 * delta_mean_r2
    - ref_mean_r2 * delta_mean_r1 - delta_mean_r1 * delta_mean_r2;
}


Real HierarchInterpPolyApproximation::
delta_covariance(const RealVector& x, const RealVector2DArray& r1_t1_coeffs,
		 const RealMatrix2DArray& r1_t2_coeffs,
		 const RealVector2DArray& r2_t1_coeffs,
		 const RealMatrix2DArray& r2_t2_coeffs, bool same,
		 const RealVector2DArray& r1r2_t1_coeffs,
		 const RealMatrix2DArray& r1r2_t2_coeffs,
		 const UShort3DArray& sm_mi,   const UShort4DArray& colloc_key,
		 const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  // Compute surplus for r1, r2, and r1r2 and retrieve reference mean values
  Real  ref_mean_r1 =
    expectation(x, r1_t1_coeffs, r1_t2_coeffs, sm_mi, colloc_key, ref_key),
      delta_mean_r1 =
    expectation(x, r1_t1_coeffs, r1_t2_coeffs, sm_mi, colloc_key, incr_key),
        ref_mean_r2 = (same) ? ref_mean_r1 :
    expectation(x, r2_t1_coeffs, r2_t2_coeffs, sm_mi, colloc_key, ref_key),
      delta_mean_r2 = (same) ? delta_mean_r1 :
    expectation(x, r2_t1_coeffs, r2_t2_coeffs, sm_mi, colloc_key, incr_key),
    delta_mean_r1r2 =
    expectation(x, r1r2_t1_coeffs, r1r2_t2_coeffs, sm_mi, colloc_key, incr_key);

  // same expression as standard expansion mode case above
#ifdef DEBUG
  PCout << "\n  delta_mean_r1r2 = " << delta_mean_r1r2
	<< "\n  ref_mean_r1 = " << ref_mean_r1
	<< "  delta_mean_r1 = " << delta_mean_r1
	<< "\n  ref_mean_r2 = " << ref_mean_r2
	<< "  delta_mean_r2 = " << delta_mean_r2 << std::endl;
#endif // DEBUG
  return delta_mean_r1r2 - ref_mean_r1 * delta_mean_r2
    - ref_mean_r2 * delta_mean_r1 - delta_mean_r1 * delta_mean_r2;
}


Real HierarchInterpPolyApproximation::
delta_covariance(const RealVector& x,
		 const std::map<UShortArray, RealVector2DArray>& r1_t1c_map,
		 const std::map<UShortArray, RealMatrix2DArray>& r1_t2c_map,
		 const std::map<UShortArray, RealVector2DArray>& r2_t1c_map,
		 const std::map<UShortArray, RealMatrix2DArray>& r2_t2c_map,
		 bool same, const RealVector2DArray& r1r2_t1c,
		 const RealMatrix2DArray& r1r2_t2c,
		 const std::map<UShortArray, UShort3DArray>& sm_mi_map,
		 const std::map<UShortArray, UShort4DArray>& colloc_key_map,
		 const UShortArray& active_key,
		 const std::map<UShortArray, UShort2DArray>& ref_key_map,
		 const std::map<UShortArray, UShort2DArray>& incr_key_map)
{
  // Compute surplus for r1, r2, and r1r2 and retrieve reference mean values
  std::map<UShortArray, RealVector2DArray>::const_iterator r1t1_cit
    = r1_t1c_map.find(active_key);
  std::map<UShortArray, RealMatrix2DArray>::const_iterator r1t2_cit
    = r1_t2c_map.find(active_key);
  std::map<UShortArray, UShort3DArray>::const_iterator sm_cit
    = sm_mi_map.find(active_key);
  std::map<UShortArray, UShort4DArray>::const_iterator ck_cit
    = colloc_key_map.find(active_key);
  std::map<UShortArray, UShort2DArray>::const_iterator incr_cit
    = incr_key_map.find(active_key);
  if (r1t1_cit == r1_t1c_map.end() || r1t2_cit == r1_t2c_map.end() ||
      sm_cit   == sm_mi_map.end()  || ck_cit   == colloc_key_map.end() ||
      incr_cit == incr_key_map.end()) {
    PCerr << "Error: failure in key lookup in HierarchInterpPolyApproximation"
	  << "::delta_covariance()" << std::endl;
    abort_handler(-1);
  }
  const UShort3DArray& active_sm_mi      =   sm_cit->second;
  const UShort4DArray& active_colloc_key =   ck_cit->second;
  const UShort2DArray& active_incr_key   = incr_cit->second;

  // increments could in general be defined by incr_key_map, but currently
  // are limited to the active incr_key (enabling short cuts below)
  Real  ref_mean_r1 =
    expectation(x, r1_t1c_map, r1_t2c_map,sm_mi_map,colloc_key_map,ref_key_map),
      delta_mean_r1 =
  //expectation(x, r1_t1c_map, r1_t2c_map, sm_mi_map, colloc_key_map,
  //            incr_key_map),
    expectation(x, r1t1_cit->second, r1t2_cit->second, active_sm_mi,
		active_colloc_key, active_incr_key); // shortcut
  Real ref_mean_r2, delta_mean_r2;
  if (same)
    { ref_mean_r2 = ref_mean_r1; delta_mean_r2 = delta_mean_r1; }
  else {
    std::map<UShortArray, RealVector2DArray>::const_iterator r2t1_cit
      = r2_t1c_map.find(active_key);
    std::map<UShortArray, RealMatrix2DArray>::const_iterator r2t2_cit
      = r2_t2c_map.find(active_key);
    ref_mean_r2 =
      expectation(x, r2_t1c_map, r2_t2c_map, sm_mi_map, colloc_key_map,
		  ref_key_map);
    delta_mean_r2 =
    //expectation(x, r2_t1c_map, r2_t2c_map, sm_mi_map, colloc_key_map,
    //            incr_key_map),
      expectation(x, r2t1_cit->second, r2t2_cit->second, active_sm_mi,
		  active_colloc_key, active_incr_key); // shortcut
  }
  Real delta_mean_r1r2 =
    expectation(x, r1r2_t1c, r1r2_t2c, active_sm_mi, active_colloc_key,
		active_incr_key);

  // same expression as standard expansion mode case above
#ifdef DEBUG
  PCout << "\n  delta_mean_r1r2 = " << delta_mean_r1r2
	<< "\n  ref_mean_r1 = " << ref_mean_r1
	<< "  delta_mean_r1 = " << delta_mean_r1
	<< "\n  ref_mean_r2 = " << ref_mean_r2
	<< "  delta_mean_r2 = " << delta_mean_r2 << std::endl;
#endif // DEBUG
  return delta_mean_r1r2 - ref_mean_r1 * delta_mean_r2
    - ref_mean_r2 * delta_mean_r1 - delta_mean_r1 * delta_mean_r2;
}


Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs, const RealVector2DArray& t1_wts,
	    const RealMatrix2DArray& t2_wts, const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real integral = 0.;
  size_t lev, set, pt, num_lev = t1_coeffs.size(),
    set_start = 0, set_end, num_tp_pts;
  bool partial = !set_partition.empty();
  if (data_rep->basisConfigOptions.useDerivs) {
    size_t v, num_v = sharedDataRep->numVars;
    for (lev=0; lev<num_lev; ++lev) {
      const RealVectorArray& t1c_l = t1_coeffs[lev];
      const RealMatrixArray& t2c_l = t2_coeffs[lev];
      const RealVectorArray& t1w_l = t1_wts[lev];
      const RealMatrixArray& t2w_l = t2_wts[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1c_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1c_ls = t1c_l[set];
	const RealMatrix& t2c_ls = t2c_l[set];
	const RealVector& t1w_ls = t1w_l[set];
	const RealMatrix& t2w_ls = t2w_l[set];
	num_tp_pts = t1c_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	  integral += t1c_ls[pt] * t1w_ls[pt];
	  const Real* t2c_lsp = t2c_ls[pt];
	  const Real* t2w_lsp = t2w_ls[pt];
	  for (v=0; v<num_v; ++v)
	    integral += t2c_lsp[v] * t2w_lsp[v];
	}
      }
    }
  }
  else {
    for (lev=0; lev<num_lev; ++lev) {
      const RealVectorArray& t1c_l = t1_coeffs[lev];
      const RealVectorArray& t1w_l = t1_wts[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1c_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1c_ls = t1c_l[set];
	const RealVector& t1w_ls = t1w_l[set];
	num_tp_pts = t1c_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) // omitted if empty surplus vector
	  integral += t1c_ls[pt] * t1w_ls[pt];
      }
    }
  }
  return integral;
}


/*
Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs,
	    const UShort3DArray& pt_partition)
{
  Real integral = 0.;
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const RealVector2DArray& t1_wts
    = hsg_driver->type1_hierarchical_weight_sets();
  size_t lev, set, pt, num_lev = t1_coeffs.size(), num_sets,
    tp_pt_start = 0, tp_pt_end;
  bool partial = !pt_partition.empty();
  switch (data_rep->basisConfigOptions.useDerivs) {
  case false:
    for (lev=0; lev<num_lev; ++lev) {
      const RealVectorArray& t1c_l = t1_coeffs[lev];
      num_sets = t1c_l.size();
      for (set=0; set<num_sets; ++set) {
	const RealVector& t1c_ls = t1c_l[set];
	const RealVector& t1w_ls = t1_wts[lev][set];
	if (partial) {
	  tp_pt_start = pt_partition[lev][set][0];
	  tp_pt_end   = pt_partition[lev][set][1];
	}
	else
	  tp_pt_end   = t1c_ls.length();
	for (pt=tp_pt_start; pt<tp_pt_end; ++pt)
	  integral += t1c_ls[pt] * t1w_ls[pt];
      }
    }
    break;
  case true: {
    const RealMatrix2DArray& t2_wts
      = hsg_driver->type2_hierarchical_weight_sets();
    size_t v, num_v = sharedDataRep->numVars;
    for (lev=0; lev<num_lev; ++lev) {
      const RealVectorArray& t1c_l = t1_coeffs[lev];
      num_sets = t1c_l.size();
      for (set=0; set<num_sets; ++set) {
	const RealVector& t1c_ls = t1c_l[set];
	const RealMatrix& t2c_ls = t2_coeffs[lev][set];
	const RealVector& t1w_ls = t1_wts[lev][set];
	const RealMatrix& t2w_ls = t2_wts[lev][set];
	if (partial) {
	  tp_pt_start = pt_partition[lev][set][0];
	  tp_pt_end   = pt_partition[lev][set][1];
	}
	else
	  tp_pt_end   = t1c_ls.length();
	for (pt=tp_pt_start; pt<tp_pt_end; ++pt) {
	  integral += t1c_ls[pt] * t1w_ls[pt];
	  const Real* t2c_lsp = t2c_ls[pt];
	  const Real* t2w_lsp = t2w_ls[pt];
	  for (v=0; v<num_v; ++v)
	    integral += t2c_lsp[v] * t2w_lsp[v];
	}
      }
    }
    break;
  }
  }

  return integral;
}
*/


Real HierarchInterpPolyApproximation::
expectation(const std::map<UShortArray, RealVector2DArray>& t1_coeffs_map,
	    const std::map<UShortArray, RealMatrix2DArray>& t2_coeffs_map,
	    const std::map<UShortArray, RealVector2DArray>& t1_wts_map,
	    const std::map<UShortArray, RealMatrix2DArray>& t2_wts_map,
	    const std::map<UShortArray, UShort2DArray>& set_partition_map)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  std::map<UShortArray, RealVector2DArray>::const_iterator t1c_cit, t1w_cit;
  std::map<UShortArray, RealMatrix2DArray>::const_iterator t2c_cit, t2w_cit;
  std::map<UShortArray, UShort2DArray>::const_iterator p_cit;

  Real integral = 0.;
  for (t1c_cit = t1_coeffs_map.begin(),      t2c_cit  = t2_coeffs_map.begin(),
       t1w_cit = t1_wts_map.begin(),         t2w_cit  = t2_wts_map.begin(),
       p_cit   = set_partition_map.begin();  t1c_cit != t1_coeffs_map.end();
       ++t1c_cit, ++t2c_cit, ++t1w_cit, ++t2w_cit, ++p_cit)
    integral += expectation(t1c_cit->second, t2c_cit->second, t1w_cit->second,
			    t2w_cit->second, p_cit->second);
  return integral;
}


Real HierarchInterpPolyApproximation::
expectation(const std::map<UShortArray, std::map<PolynomialApproximation*,
	      RealVector2DArray> >& prod_t1c_map,
	    const std::map<UShortArray, std::map<PolynomialApproximation*,
	      RealMatrix2DArray> >& prod_t2c_map,
	    PolynomialApproximation* poly_approx_2,
	    const std::map<UShortArray, RealVector2DArray>& t1_wts_map,
	    const std::map<UShortArray, RealMatrix2DArray>& t2_wts_map,
	    const std::map<UShortArray, UShort2DArray>& set_partition_map)
{
  std::map<UShortArray, std::map<PolynomialApproximation*, RealVector2DArray> >
    ::const_iterator p1c_cit;
  std::map<UShortArray, std::map<PolynomialApproximation*, RealMatrix2DArray> >
    ::const_iterator p2c_cit;
  std::map<PolynomialApproximation*, RealVector2DArray>::const_iterator t1c_cit;
  std::map<PolynomialApproximation*, RealMatrix2DArray>::const_iterator t2c_cit;
  std::map<UShortArray, RealVector2DArray>::const_iterator t1w_cit;
  std::map<UShortArray, RealMatrix2DArray>::const_iterator t2w_cit;
  std::map<UShortArray, UShort2DArray>::const_iterator p_cit;

  Real integral = 0.;
  for (p1c_cit = prod_t1c_map.begin(),       p2c_cit  = prod_t2c_map.begin(),
       t1w_cit = t1_wts_map.begin(),         t2w_cit  = t2_wts_map.begin(),
       p_cit   = set_partition_map.begin();  p1c_cit != prod_t1c_map.end();
       ++p1c_cit, ++p2c_cit, ++t1w_cit, ++t2w_cit, ++p_cit) {
    t1c_cit = p1c_cit->second.find(poly_approx_2);
    t2c_cit = p2c_cit->second.find(poly_approx_2);
    integral += expectation(t1c_cit->second, t2c_cit->second, t1w_cit->second,
			    t2w_cit->second, p_cit->second);
  }
  return integral;
}


Real HierarchInterpPolyApproximation::
expectation(const RealVector& x, const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs, const UShort3DArray& sm_mi,
	    const UShort4DArray& colloc_key, const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real integral = 0.;
  size_t lev, set, pt, num_lev = t1_coeffs.size(), set_start = 0, set_end,
    num_tp_pts;
  bool partial = !set_partition.empty();
  const SizetList&  r_ind = data_rep->randomIndices;
  const SizetList& nr_ind = data_rep->nonRandomIndices;
  if (data_rep->basisConfigOptions.useDerivs) {
    size_t v, num_v = sharedDataRep->numVars;
    for (lev=0; lev<num_lev; ++lev) {
      const RealVectorArray& t1c_l = t1_coeffs[lev];
      const RealMatrixArray& t2c_l = t2_coeffs[lev];
     if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1c_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector&    t1c_ls = t1c_l[set];
	const RealMatrix&    t2c_ls = t2c_l[set];
	const UShortArray& sm_mi_ls = sm_mi[lev][set];
	num_tp_pts = t1c_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	  const UShortArray& key_lsp = colloc_key[lev][set][pt];
	  integral += t1c_ls[pt]
	    * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls, nr_ind)
	    * data_rep->type1_weight(key_lsp, sm_mi_ls, r_ind);
	  const Real* t2c_lsp = t2c_ls[pt];
	  for (v=0; v<num_v; ++v)
	    integral += t2c_lsp[v]
	      * data_rep->type2_interpolant_value(x, v, key_lsp,sm_mi_ls,nr_ind)
	      * data_rep->type2_weight(v, key_lsp, sm_mi_ls, r_ind);
	}
      }
    }
  }
  else {
    for (lev=0; lev<num_lev; ++lev) {
      const RealVectorArray& t1c_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1c_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector&    t1c_ls = t1c_l[set];
	const UShortArray& sm_mi_ls = sm_mi[lev][set];
	num_tp_pts = t1c_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	  const UShortArray& key_lsp = colloc_key[lev][set][pt];
	  integral += t1c_ls[pt]
	    * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls, nr_ind)
	    * data_rep->type1_weight(key_lsp, sm_mi_ls, r_ind);
	}
      }
    }
  }
  return integral;
}


Real HierarchInterpPolyApproximation::
expectation(const RealVector& x,
	    const std::map<UShortArray, RealVector2DArray>& t1_coeffs_map,
	    const std::map<UShortArray, RealMatrix2DArray>& t2_coeffs_map,
	    const std::map<UShortArray, UShort3DArray>& sm_mi_map,
	    const std::map<UShortArray, UShort4DArray>& colloc_key_map,
	    const std::map<UShortArray, UShort2DArray>& set_partition_map)
{
  std::map<UShortArray, RealVector2DArray>::const_iterator t1c_cit;
  std::map<UShortArray, RealMatrix2DArray>::const_iterator t2c_cit;
  std::map<UShortArray, UShort3DArray>::const_iterator sm_cit;
  std::map<UShortArray, UShort4DArray>::const_iterator ck_cit;
  std::map<UShortArray, UShort2DArray>::const_iterator p_cit;

  Real integral = 0.;
  for (t1c_cit = t1_coeffs_map.begin(),     t2c_cit  = t2_coeffs_map.begin(),
       sm_cit  = sm_mi_map.begin(),          ck_cit  = colloc_key_map.begin(),
       p_cit   = set_partition_map.begin(); t1c_cit != t1_coeffs_map.end();
       ++t1c_cit, ++t2c_cit, ++sm_cit, ++ck_cit, ++p_cit)
    integral += expectation(x, t1c_cit->second, t2c_cit->second, sm_cit->second,
			    ck_cit->second, p_cit->second);
  return integral;
}


Real HierarchInterpPolyApproximation::
expectation(const RealVector& x,
	    const std::map<UShortArray, std::map<PolynomialApproximation*,
	      RealVector2DArray> >& prod_t1c_map,
	    const std::map<UShortArray, std::map<PolynomialApproximation*,
	      RealMatrix2DArray> >& prod_t2c_map,
	    PolynomialApproximation* poly_approx_2,
	    const std::map<UShortArray, UShort3DArray>& sm_mi_map,
	    const std::map<UShortArray, UShort4DArray>& colloc_key_map,
	    const std::map<UShortArray, UShort2DArray>& set_partition_map)
{
  std::map<UShortArray, std::map<PolynomialApproximation*, RealVector2DArray> >
    ::const_iterator p1c_cit;
  std::map<UShortArray, std::map<PolynomialApproximation*, RealMatrix2DArray> >
    ::const_iterator p2c_cit;
  std::map<PolynomialApproximation*, RealVector2DArray>::const_iterator t1c_cit;
  std::map<PolynomialApproximation*, RealMatrix2DArray>::const_iterator t2c_cit;
  std::map<UShortArray, UShort3DArray>::const_iterator sm_cit;
  std::map<UShortArray, UShort4DArray>::const_iterator ck_cit;
  std::map<UShortArray, UShort2DArray>::const_iterator p_cit;
  
  Real integral = 0.;
  for (p1c_cit = prod_t1c_map.begin(),      p2c_cit  = prod_t2c_map.begin(),
       sm_cit  = sm_mi_map.begin(),          ck_cit  = colloc_key_map.begin(),
       p_cit   = set_partition_map.begin(); p1c_cit != prod_t1c_map.end();
       ++p1c_cit, ++p2c_cit, ++sm_cit, ++ck_cit, ++p_cit) {
    t1c_cit = p1c_cit->second.find(poly_approx_2);
    t2c_cit = p2c_cit->second.find(poly_approx_2);
    integral += expectation(x, t1c_cit->second, t2c_cit->second, sm_cit->second,
			    ck_cit->second, p_cit->second);
  }
  return integral;
}


/** For inserted/augmented design/epistemic variables in standard mode. */
const RealVector& HierarchInterpPolyApproximation::
expectation_gradient(const RealMatrix2DArray& t1_coeff_grads,
		     const RealVector2DArray& t1_wts)
{
  size_t lev, num_lev = t1_coeff_grads.size(), set, num_sets, pt, num_tp_pts,
    v, num_deriv_vars = t1_coeff_grads[0][0].numRows();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.;

  for (lev=0; lev<num_lev; ++lev) {
    const RealMatrixArray& t1_coeff_grads_l = t1_coeff_grads[lev];
    num_sets = t1_coeff_grads_l.size();
    for (set=0; set<num_sets; ++set) {
      const RealMatrix& t1_coeff_grads_ls = t1_coeff_grads_l[set];
      num_tp_pts = t1_coeff_grads_ls.numCols();
      for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	const Real* t1_coeff_grads_lsp = t1_coeff_grads_ls[pt];
	Real t1_wt_lsp = t1_wts[lev][set][pt];
	for (v=0; v<num_deriv_vars; ++v)
	  approxGradient[v] += t1_coeff_grads_lsp[v] * t1_wt_lsp;
      }
    }
  }

  return approxGradient;
}


/** For inserted design/epistemic variables in all_variables mode. */
Real HierarchInterpPolyApproximation::
expectation_gradient(const RealVector& x,
		     const RealMatrix2DArray& t1_coeff_grads,
		     const UShort3DArray& sm_mi,
		     const UShort4DArray& colloc_key, size_t t1cg_index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, num_lev = t1_coeff_grads.size(), set, num_sets, pt, num_tp_pts;
  const SizetList&  r_ind = data_rep->randomIndices;
  const SizetList& nr_ind = data_rep->nonRandomIndices;

  Real grad = 0.;
  for (lev=0; lev<num_lev; ++lev) {
    const RealMatrixArray& t1_coeff_grads_l = t1_coeff_grads[lev];
    num_sets = t1_coeff_grads_l.size();
    for (set=0; set<num_sets; ++set) {
      const RealMatrix& t1_coeff_grads_ls = t1_coeff_grads_l[set];
      const UShortArray&         sm_mi_ls = sm_mi[lev][set];
      num_tp_pts = t1_coeff_grads_ls.numCols();
      for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	const UShortArray& key_lsp = colloc_key[lev][set][pt];
	grad += t1_coeff_grads_ls(t1cg_index, pt)
	  * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls, nr_ind)
	  * data_rep->type1_weight(key_lsp, sm_mi_ls, r_ind);
      }
    }
  }
  return grad;
}


/** For augmented design/epistemic variables in all_variables mode. */
Real HierarchInterpPolyApproximation::
expectation_gradient(const RealVector& x, const RealVector2DArray& t1_coeffs,
		     const RealMatrix2DArray& t2_coeffs,
		     const UShort3DArray& sm_mi,
		     const UShort4DArray& colloc_key, size_t deriv_index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, num_lev = t1_coeffs.size(), set, num_sets, pt, num_tp_pts, v,
    num_v = sharedDataRep->numVars;
  const SizetList&  r_ind = data_rep->randomIndices;
  const SizetList& nr_ind = data_rep->nonRandomIndices;

  Real grad = 0.;
  for (lev=0; lev<num_lev; ++lev) {
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    num_sets = t1_coeffs_l.size();
    for (set=0; set<num_sets; ++set) {
      const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
      const UShortArray&    sm_mi_ls = sm_mi[lev][set];
      num_tp_pts = t1_coeffs_ls.length();
      for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	const UShortArray& key_lsp = colloc_key[lev][set][pt];
	grad += t1_coeffs_ls[pt] * data_rep->
	  type1_interpolant_gradient(x, deriv_index, key_lsp, sm_mi_ls,
				     data_rep->nonRandomIndices) *
	  data_rep->type1_weight(key_lsp, sm_mi_ls, data_rep->randomIndices);
	if (data_rep->basisConfigOptions.useDerivs) {
	  const Real *t2_coeff_lsp = t2_coeffs[lev][set][pt];
	  for (v=0; v<num_v; ++v)
	    grad += t2_coeff_lsp[v] *
	      data_rep->type2_interpolant_gradient(x, deriv_index, v, key_lsp,
						   sm_mi_ls, nr_ind) *
	      data_rep->type2_weight(v, key_lsp, sm_mi_ls, r_ind);
	}
      }
    }
  }

  return grad;
}


/** This version accesses the surrogate data for efficiency in cases where
    the r1 and r2 interpolants correspond to the data and do not need to be
    evaluated.  Note: could consider computing deltaR1R2 = R1 deltaR2 +
    R2 deltaR1 + deltaR1 deltaR2. */
void HierarchInterpPolyApproximation::
product_interpolant(const SDVArray& sdv_array, const SDRArray& sdr_array_1,
		    const SDRArray& sdr_array_2, const UShort3DArray& sm_mi,
		    const UShort4DArray& colloc_key,
		    const Sizet3DArray& colloc_index,
		    RealVector2DArray& prod_t1c, RealMatrix2DArray& prod_t2c,
		    const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_lev = colloc_key.size(), set_start = 0, set_end,
    num_sets, num_tp_pts, cntr = 0, c_index, v, num_v = sharedDataRep->numVars;
  bool partial = !set_partition.empty(), empty_c_index = colloc_index.empty(),
    use_derivs = data_rep->basisConfigOptions.useDerivs;
  Real data_fn1, data_fn2;

  // form hierarchical t1/t2 coeffs for raw moment R1 R2

  // level 0 (no surplus)
  prod_t1c.resize(num_lev);  prod_t2c.resize(num_lev);
  if (!partial || set_partition[0][0] == 0) { // partition includes level 0
    prod_t1c[0].resize(1);  prod_t2c[0].resize(1);
    prod_t1c[0][0].sizeUninitialized(1);
    c_index = (empty_c_index) ? cntr : colloc_index[0][0][0];
    const SurrogateDataResp& sdr0_1 = sdr_array_1[c_index];
    const SurrogateDataResp& sdr0_2 = sdr_array_2[c_index];
    data_fn1 = sdr0_1.response_function();
    data_fn2 = sdr0_2.response_function();
    prod_t1c[0][0][0] = data_fn1 * data_fn2;
    if (use_derivs) {
      prod_t2c[0][0].shapeUninitialized(num_v, 1);
      Real *prod_t2c_000 = prod_t2c[0][0][0];
      const RealVector& data_grad1 = sdr0_1.response_gradient();
      const RealVector& data_grad2 = sdr0_2.response_gradient();
      for (v=0; v<num_v; ++v)
	prod_t2c_000[v]
	  = data_fn1 * data_grad2[v] + data_fn2 * data_grad1[v];
    }
#ifdef DEBUG
    PCout << "Surplus components l0 s0 p0: r1r2 = " << prod_t1c[0][0][0]
	  << std::endl;
#endif // DEBUG
  }
  ++cntr;

  // levels 1:w (subtract lev-1 for surplus estimation)
  for (lev=1; lev<num_lev; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];  num_sets = key_l.size();
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else
      set_end = num_sets;
    prod_t1c[lev].resize(num_sets);  prod_t2c[lev].resize(num_sets);
    if (empty_c_index)
      for (set=0; set<set_start; ++set)
	cntr += key_l[set].size();
    for (set=set_start; set<set_end; ++set) {
      num_tp_pts = key_l[set].size();
      RealVector& prod_t1c_ls = prod_t1c[lev][set];
      prod_t1c_ls.sizeUninitialized(num_tp_pts);
      if (use_derivs) { // type{1,2} interpolation
	RealMatrix& prod_t2c_ls = prod_t2c[lev][set];
	prod_t2c_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt) {
	  c_index = (empty_c_index) ? cntr++ : colloc_index[lev][set][pt];
	  const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	  const SurrogateDataResp& sdr_1 = sdr_array_1[c_index];
	  const SurrogateDataResp& sdr_2 = sdr_array_2[c_index];
	  // type1 hierarchical interpolation of R1 R2
	  data_fn1 = sdr_1.response_function();
	  data_fn2 = sdr_2.response_function();
	  prod_t1c_ls[pt] = data_fn1 * data_fn2
	    - value(c_vars, sm_mi, colloc_key, prod_t1c, prod_t2c, lev-1);
	  // type2 hierarchical interpolation of R1 R2
	  // --> interpolated grads are R1 * R2' + R2 * R1'
	  Real* prod_t2c_lsp = prod_t2c_ls[pt];
	  const RealVector& data_grad1 = sdr_1.response_gradient();
	  const RealVector& data_grad2 = sdr_2.response_gradient();
	  const RealVector& prev_grad  = gradient_basis_variables(c_vars,
	    sm_mi, colloc_key, prod_t1c, prod_t2c, lev-1);
	  for (v=0; v<num_v; ++v)
	    prod_t2c_lsp[v] = data_fn1 * data_grad2[v]
	      + data_fn2 * data_grad1[v] - prev_grad[v];
	}
      }
      else // type1 hierarchical interpolation of R1 R2
	for (pt=0; pt<num_tp_pts; ++pt) {
	  c_index = (empty_c_index) ? cntr++ : colloc_index[lev][set][pt];
#ifdef DEBUG
	  PCout << "c_index = " << c_index << std::endl;
	  Real r1 = sdr_array_1[c_index].response_function(),
	       r2 = sdr_array_2[c_index].response_function(), r1r2 = r1*r2,
	       r1r2_lm1 = value(sdv_array[c_index].continuous_variables(),
		 sm_mi, colloc_key, prod_t1c, prod_t2c, lev-1),
	       surplus = r1r2 - r1r2_lm1;
	  PCout << "Surplus components l" << lev << " s" << set << " p" << pt
		<< ": r1r2 = " << r1r2 << " r1r2_lm1 = " << r1r2_lm1
		<< " surplus = " << surplus << std::endl;
	  prod_t1c_ls[pt] = surplus;
#else
	  prod_t1c_ls[pt] = sdr_array_1[c_index].response_function()
	    * sdr_array_2[c_index].response_function()
	    - value(sdv_array[c_index].continuous_variables(), sm_mi,
		    colloc_key, prod_t1c, prod_t2c, lev-1);
#endif // DEBUG
	}
    }
  }

#ifdef DEBUG
  PCout << "Product interpolant type1 coeffs:\n" << prod_t1c;
  if (use_derivs)
    PCout << "Product interpolant type2 coeffs:\n" << prod_t2c;
#endif // DEBUG
}


/** This version evaluates the r1 and r2 interpolants for cases where
    a particular key of the surrogate data does not correspond (e.g.,
    after combine_coefficients).  Note: consider computing deltaR1R2 =
    R1 deltaR2 + R2 deltaR1 + deltaR1 deltaR2. */
void HierarchInterpPolyApproximation::
product_interpolant(const RealMatrix2DArray& var_sets,
		    const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
		    const RealVector2DArray& r1_t1c,
		    const RealMatrix2DArray& r1_t2c,
		    const RealVector2DArray& r2_t1c,
		    const RealMatrix2DArray& r2_t2c, bool same,
		    RealVector2DArray& prod_t1c, RealMatrix2DArray& prod_t2c,
		    const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_lev = r1_t1c.size(), set_start = 0, set_end,
    num_sets, num_tp_pts, v, num_v = sharedDataRep->numVars;
  bool partial = !set_partition.empty(),
    use_derivs = data_rep->basisConfigOptions.useDerivs;
  Real r1_val, r2_val;

  // form hierarchical t1/t2 coeffs for raw moment R1 R2

  // level 0 (no surplus)
  prod_t1c.resize(num_lev);  prod_t2c.resize(num_lev);
  if (!partial || set_partition[0][0] == 0) { // partition includes level 0
    prod_t1c[0].resize(1); prod_t2c[0].resize(1);
    prod_t1c[0][0].sizeUninitialized(1);
    RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		      (int)num_v);
    // type1 hierarchical interpolation of R1 R2
    r1_val = value(c_vars, sm_mi, colloc_key, r1_t1c, r1_t2c, 0);
    r2_val = (same) ? r1_val :
      value(c_vars, sm_mi, colloc_key, r2_t1c, r2_t2c, 0);
    prod_t1c[0][0][0] = r1_val * r2_val;
    if (use_derivs) {
      const RealVector& r1_grad = gradient_basis_variables(c_vars, sm_mi,
        colloc_key, r1_t1c, r1_t2c, 0);
      const RealVector& r2_grad = (same) ? r1_grad : gradient_basis_variables(
        c_vars, sm_mi, colloc_key, r2_t1c, r2_t2c, 0);
      prod_t2c[0][0].shapeUninitialized(num_v, 1);
      Real *prod_t2c_000 = prod_t2c[0][0][0];
      for (v=0; v<num_v; ++v)
	prod_t2c_000[v] = r1_val * r2_grad[v] + r2_val * r1_grad[v];
    }
  }

  // levels 1:w (subtract lev-1 for surplus estimation)
  for (lev=1; lev<num_lev; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];  num_sets = key_l.size();
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else
      set_end = num_sets;
    prod_t1c[lev].resize(num_sets);  prod_t2c[lev].resize(num_sets);
    for (set=set_start; set<set_end; ++set) {
      num_tp_pts = key_l[set].size();
      RealVector& prod_t1c_ls = prod_t1c[lev][set];
      prod_t1c_ls.sizeUninitialized(num_tp_pts);
      if (use_derivs) {
	RealMatrix& prod_t2c_ls = prod_t2c[lev][set];
	prod_t2c_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt) {
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[lev][set][pt]), (int)num_v);
	  // type1 hierarchical interpolation of R1 R2
	  r1_val = value(c_vars, sm_mi, colloc_key, r1_t1c, r1_t2c, lev);
	  r2_val = (same) ? r1_val :
	    value(c_vars, sm_mi, colloc_key, r2_t1c, r2_t2c, lev);
	  prod_t1c_ls[pt] = r1_val * r2_val -
	    value(c_vars, sm_mi, colloc_key, prod_t1c, prod_t2c, lev-1);
	  // type2 hierarchical interpolation of R1 R2
	  // --> interpolated grads are R1 * R2' + R2 * R1'
	  Real* prod_t2c_lsp = prod_t2c_ls[pt];
	  const RealVector& r1_grad = gradient_basis_variables(c_vars, sm_mi,
	    colloc_key, r1_t1c, r1_t2c, lev);
	  const RealVector& r2_grad = (same) ? r1_grad :
	    gradient_basis_variables(c_vars, sm_mi, colloc_key, r2_t1c,
	    r2_t2c, lev);
	  const RealVector& prev_grad = gradient_basis_variables(c_vars,
	    sm_mi, colloc_key, prod_t1c, prod_t2c, lev-1);
	  for (v=0; v<num_v; ++v)
	    prod_t2c_lsp[v] = r1_val * r2_grad[v] + r2_val * r1_grad[v]
	                    - prev_grad[v];
	}
      }
      else // type1 hierarchical interpolation of R1 R2
	for (pt=0; pt<num_tp_pts; ++pt) {
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[lev][set][pt]), (int)num_v);
	  // type1 hierarchical interpolation of R1 R2
	  r1_val = value(c_vars, sm_mi, colloc_key, r1_t1c, r1_t2c, lev);
	  r2_val = (same) ? r1_val :
	    value(c_vars, sm_mi, colloc_key, r2_t1c, r2_t2c, lev);
	  prod_t1c_ls[pt] = r1_val * r2_val -
	    value(c_vars, sm_mi, colloc_key, prod_t1c, prod_t2c, lev-1);
	}
    }
  }

#ifdef DEBUG
  PCout << "Product interpolant type1 coeffs:\n" << prod_t1c;
  if (use_derivs)
    PCout << "Product interpolant type2 coeffs:\n" << prod_t2c;
#endif // DEBUG
}


/** This version accesses the raw surrogate data for Q^H and Q^L for two 
    QoI in order to hierarchically interpolate the difference of products
    Q^H_1 Q^H_2 - Q^L_1 Q^L_2.  The goal is to have the process for 
    interpolating Q^2 mirror the process for interpolating Q, but with 
    focus on combined covariance (we care about the effect of a level 
    increment on the combined QoI). */
void HierarchInterpPolyApproximation::
product_difference_interpolant(const SurrogateData& surr_data_1,
			       const SurrogateData& surr_data_2,
			       const UShort3DArray& sm_mi,
			       const UShort4DArray& colloc_key,
			       const Sizet3DArray&  colloc_index,
			       RealVector2DArray&   prod_t1c,
			       RealMatrix2DArray&   prod_t2c,
			       const UShortArray&   lf_key,
			       const UShort2DArray& set_partition)
{
  // Hierarchically interpolate R_1 * R_2 exactly as for R, i.e., hierarch
  // interp of discrepancies in R_1 * R_2 for model indices > first
  // > For R, the differences are already in modSurrData as computed by
  //   PolynomialApproximation::response_data_to_discrepancy_data()
  // **********************************************************************
  // > TO DO: carrying delta(R_1 R_2) and delta((R-1-mu_1)(R_2-mu_2))
  //          within additional modSurrData keys could eliminate need to
  //          recompute, but unlike R, all covariance pairs must be tracked.
  // **********************************************************************
  // > only need to integrate terms in incr_key, but need to build up to them
  //   by evaluating all underlying contribs to the hierarchical interpolant

  // This case is _not_ modular on SurrogateData instance:
  // it must use surrData to access lower level data
  const SDVArray& sdv_array      = surr_data_1.variables_data();
  const SDRArray& hf_sdr_array_1 = surr_data_1.response_data();
  const SDRArray& hf_sdr_array_2 = surr_data_2.response_data();

  // Accommodate level 0 --> lf_key is empty
  if (lf_key.empty()) {
    product_interpolant(sdv_array, hf_sdr_array_1, hf_sdr_array_2, sm_mi,
			colloc_key, colloc_index, prod_t1c, prod_t2c,
			set_partition);
    return;
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_lev = colloc_key.size(), set_start = 0, set_end,
    num_sets, num_tp_pts, cntr = 0, c_index, v, num_v = sharedDataRep->numVars;
  bool partial = !set_partition.empty(), empty_c_index = colloc_index.empty(),
    use_derivs = data_rep->basisConfigOptions.useDerivs,
    same = (surr_data_1.data_rep() == surr_data_2.data_rep());
  Real hf_fn1, lf_fn1, hf_fn2, lf_fn2, r1r2_l;

  // Support original (R) data for computing increments in R^2
  std::map<UShortArray, SDRArray>::const_iterator
    r_cit_1 = surr_data_1.response_data_map().find(lf_key),
    r_cit_2 = (same) ? r_cit_1 : surr_data_2.response_data_map().find(lf_key);
  const SDRArray& lf_sdr_array_1 = r_cit_1->second;
  const SDRArray& lf_sdr_array_2 = r_cit_2->second;

  // form hierarchical t1/t2 coeffs for raw moment R1 R2

  // level 0 (no surplus)
  prod_t1c.resize(num_lev);  prod_t2c.resize(num_lev);
  if (!partial || set_partition[0][0] == 0) { // partition includes level 0
    prod_t1c[0].resize(1);  prod_t2c[0].resize(1);
    prod_t1c[0][0].sizeUninitialized(1);
    // delta [ R_1 * R_2 ], not deltaR_1 * deltaR_2
    c_index = (empty_c_index) ? cntr : colloc_index[0][0][0];
    hf_fn1 = hf_sdr_array_1[c_index].response_function();
    lf_fn1 = lf_sdr_array_1[c_index].response_function();
    if (same)
      prod_t1c[0][0][0] = hf_fn1 * hf_fn1 - lf_fn1 * lf_fn1;
    else {
      hf_fn2 = hf_sdr_array_2[c_index].response_function();
      lf_fn2 = lf_sdr_array_2[c_index].response_function();
      prod_t1c[0][0][0] = hf_fn1 * hf_fn2 - lf_fn1 * lf_fn2;
    }
    if (use_derivs) {
      prod_t2c[0][0].shapeUninitialized(num_v, 1);
      Real *prod_t2c_000 = prod_t2c[0][0][0];
      const RealVector& hf_grad1 = hf_sdr_array_1[c_index].response_gradient();
      const RealVector& lf_grad1 = lf_sdr_array_1[c_index].response_gradient();
      if (same)
	for (v=0; v<num_v; ++v)
	  prod_t2c_000[v] = 2. * (hf_fn1 * hf_grad1[v] - lf_fn1 * lf_grad1[v]);
      else {
	const RealVector& hf_grad2 =hf_sdr_array_2[c_index].response_gradient();
	const RealVector& lf_grad2 =lf_sdr_array_2[c_index].response_gradient();
	for (v=0; v<num_v; ++v)
	  prod_t2c_000[v] = hf_fn1 * hf_grad2[v] + hf_fn2 * hf_grad1[v]
	                  - lf_fn1 * lf_grad2[v] - lf_fn2 * lf_grad1[v];
      }
    }
#ifdef DEBUG
    PCout << "Surplus components l0 s0 p0: r1r2 = " << prod_t1c[0][0][0]
	  << std::endl;
#endif // DEBUG
  }
  ++cntr;

  // levels 1:w (subtract lev-1 for surplus estimation)
  for (lev=1; lev<num_lev; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];  num_sets = key_l.size();
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else
      set_end = num_sets;
    prod_t1c[lev].resize(num_sets);  prod_t2c[lev].resize(num_sets);
    if (empty_c_index)
      for (set=0; set<set_start; ++set)
	cntr += key_l[set].size();
    for (set=set_start; set<set_end; ++set) {
      num_tp_pts = key_l[set].size();
      RealVector& prod_t1c_ls = prod_t1c[lev][set];
      prod_t1c_ls.sizeUninitialized(num_tp_pts);
      if (use_derivs) {
	RealMatrix& prod_t2c_ls = prod_t2c[lev][set];
	prod_t2c_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt) {
	  c_index = (empty_c_index) ? cntr++ : colloc_index[lev][set][pt];
	  const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	  // type1 hierarchical interpolation of R1 R2
	  hf_fn1 = hf_sdr_array_1[c_index].response_function();
	  lf_fn1 = lf_sdr_array_1[c_index].response_function();
	  if (same)
	    r1r2_l = hf_fn1 * hf_fn1 - lf_fn1 * lf_fn1;
	  else {
	    hf_fn2 = hf_sdr_array_2[c_index].response_function();
	    lf_fn2 = lf_sdr_array_2[c_index].response_function();
	    r1r2_l = hf_fn1 * hf_fn2 - lf_fn1 * lf_fn2;
	  }
	  prod_t1c_ls[pt] = r1r2_l
	    - value(c_vars, sm_mi, colloc_key, prod_t1c, prod_t2c, lev-1);
	  // type2 hierarchical interpolation of R1 R2
	  // --> interpolated grads are R1 * R2' + R2 * R1'
	  Real* prod_t2c_lsp = prod_t2c_ls[pt];
	  const RealVector& hf_grad1
	    = hf_sdr_array_1[c_index].response_gradient();
	  const RealVector& lf_grad1
	    = lf_sdr_array_1[c_index].response_gradient();
	  const RealVector& prev_grad = gradient_basis_variables(c_vars,
	    sm_mi, colloc_key, prod_t1c, prod_t2c, lev-1);
	  if (same)
	    for (v=0; v<num_v; ++v)
	      prod_t2c_lsp[v] = 2. *
		(hf_fn1 * hf_grad1[v] - lf_fn1 * lf_grad1[v]) - prev_grad[v];
	  else {
	    const RealVector& hf_grad2
	      = hf_sdr_array_2[c_index].response_gradient();
	    const RealVector& lf_grad2
	      = lf_sdr_array_2[c_index].response_gradient();
	    for (v=0; v<num_v; ++v)
	      prod_t2c_lsp[v]
		= hf_fn1 * hf_grad2[v] + hf_fn2 * hf_grad1[v]
		- lf_fn1 * lf_grad2[v] - lf_fn2 * lf_grad1[v] - prev_grad[v];
	  }
	}
      }
      else // type1 hierarchical interpolation of R1 R2
	for (pt=0; pt<num_tp_pts; ++pt) {
	  c_index = (empty_c_index) ? cntr++ : colloc_index[lev][set][pt];
	  hf_fn1 = hf_sdr_array_1[c_index].response_function();
	  lf_fn1 = lf_sdr_array_1[c_index].response_function();
	  if (same)
	    r1r2_l = hf_fn1 * hf_fn1 - lf_fn1 * lf_fn1;
	  else {
	    hf_fn2 = hf_sdr_array_2[c_index].response_function();
	    lf_fn2 = lf_sdr_array_2[c_index].response_function();
	    r1r2_l = hf_fn1 * hf_fn2 - lf_fn1 * lf_fn2;
	  }
#ifdef DEBUG
	  Real r1r2_lm1 =
	    value(sdv_array[c_index].continuous_variables(), sm_mi, colloc_key,
		  prod_t1c, prod_t2c, lev-1),
	        surplus = r1r2_l - r1r2_lm1;
	  PCout << "Surplus components l" << lev << " s" << set << " p" << pt
		<< ": r1r2_l = " << r1r2_l << " r1r2_lm1 = " << r1r2_lm1
		<< " surplus = " << surplus << std::endl;
	  prod_t1c_ls[pt] = surplus;
#else
	  prod_t1c_ls[pt] = r1r2_l -
	    value(sdv_array[c_index].continuous_variables(), sm_mi, colloc_key,
		  prod_t1c, prod_t2c, lev-1);
#endif // DEBUG
	}
    }
  }

#ifdef DEBUG
  PCout << "Product interpolant type1 coeffs:\n" << prod_t1c;
  if (use_derivs)
    PCout << "Product interpolant type2 coeffs:\n" << prod_t2c;
#endif // DEBUG
}


/** This version accesses the surrogate data for efficiency in cases where
    the r1 and r2 interpolants correspond to the data and do not need to be
    evaluated.  Note: whereas expectation() supports either a reference or
    increment key (passed as generic set_partition), functions forming
    hierarchical product interpolants only support a reference key due to
    nonlinearity (only end point can be controlled).  Note: could consider 
    computing deltaR1R2 = R1 deltaR2 + R2 deltaR1 + deltaR1 deltaR2. */
void HierarchInterpPolyApproximation::
central_product_interpolant(const SDVArray& sdv_array,
			    const SDRArray& sdr_array_1,
			    const SDRArray& sdr_array_2, Real mean_1,
			    Real mean_2, const UShort3DArray& sm_mi,
			    const UShort4DArray& colloc_key,
			    const Sizet3DArray& colloc_index,
			    RealVector2DArray& cov_t1c,
			    RealMatrix2DArray& cov_t2c,
			    const UShort2DArray& set_partition)
{
  // form hierarchical t1/t2 coeffs for (R_1 - \mu_1) (R_2 - \mu_2)

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;

  size_t lev, set, pt, num_lev = colloc_key.size(), set_start = 0, set_end,
    num_sets, num_tp_pts, cntr = 0, c_index, v, num_v = sharedDataRep->numVars;
  bool partial = !set_partition.empty(), empty_c_index = colloc_index.empty(),
    use_derivs = data_rep->basisConfigOptions.useDerivs;
  Real data_fn1_mm1, data_fn2_mm2;

  // level 0 (no surplus)
  cov_t1c.resize(num_lev);  cov_t2c.resize(num_lev);
  if (!partial || set_partition[0][0] == 0) { // partition includes level 0
    cov_t1c[0].resize(1);  cov_t2c[0].resize(1);
    cov_t1c[0][0].sizeUninitialized(1);
    c_index = (empty_c_index) ? cntr : colloc_index[0][0][0];
    const SurrogateDataResp& sdr0_1 = sdr_array_1[c_index];
    const SurrogateDataResp& sdr0_2 = sdr_array_2[c_index];
    data_fn1_mm1 = sdr0_1.response_function() - mean_1;
    data_fn2_mm2 = sdr0_2.response_function() - mean_2;
    cov_t1c[0][0][0] = data_fn1_mm1 * data_fn2_mm2;
    if (use_derivs) {
      cov_t2c[0][0].shapeUninitialized(num_v, 1);
      Real *cov_t2c_000 = cov_t2c[0][0][0];
      const RealVector& data_grad1 = sdr0_1.response_gradient();
      const RealVector& data_grad2 = sdr0_2.response_gradient();
      // Note: mean is only a function of nonbasis variables, so
      // basis vars grad of data_fn{1,2}_mm{1,2} is only data_grad{1,2}
      for (v=0; v<num_v; ++v)
	cov_t2c_000[v]
	  = data_fn1_mm1 * data_grad2[v] + data_fn2_mm2 * data_grad1[v];
    }
  }
  ++cntr;

  // levels 1:w (subtract lev-1 for surplus estimation)
  for (lev=1; lev<num_lev; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];  num_sets = key_l.size();
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else
      set_end = num_sets;
    cov_t1c[lev].resize(num_sets); cov_t2c[lev].resize(num_sets);
    if (empty_c_index)
      for (set=0; set<set_start; ++set)
	cntr += key_l[set].size();
    for (set=set_start; set<set_end; ++set) {
      num_tp_pts = key_l[set].size();
      RealVector& cov_t1c_ls = cov_t1c[lev][set];
      cov_t1c_ls.sizeUninitialized(num_tp_pts);
      if (use_derivs) {
	RealMatrix& cov_t2c_ls = cov_t2c[lev][set];
	cov_t2c_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt) {
	  c_index = (empty_c_index) ? cntr++ : colloc_index[lev][set][pt];
	  const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	  const SurrogateDataResp& sdr_1 = sdr_array_1[c_index];
	  const SurrogateDataResp& sdr_2 = sdr_array_2[c_index];
	  // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  data_fn1_mm1 = sdr_1.response_function() - mean_1;
	  data_fn2_mm2 = sdr_2.response_function() - mean_2;
	  cov_t1c_ls[pt] = data_fn1_mm1 * data_fn2_mm2 -
	    value(c_vars, sm_mi, colloc_key, cov_t1c, cov_t2c, lev-1);
	  // type2 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  // --> interpolated grads are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
	  Real* cov_t2c_lsp = cov_t2c_ls[pt];
	  const RealVector& data_grad1 = sdr_1.response_gradient();
	  const RealVector& data_grad2 = sdr_2.response_gradient();
	  const RealVector& prev_grad  = gradient_basis_variables(c_vars,
	    sm_mi, colloc_key, cov_t1c, cov_t2c, lev-1);
	  // Note: mean is only a function of nonbasis variables, so
	  // basis vars grad of data_fn{1,2}_mm{1,2} is only data_grad{1,2}
	  for (v=0; v<num_v; ++v)
	    cov_t2c_lsp[v] = data_fn1_mm1 * data_grad2[v]
	      + data_fn2_mm2 * data_grad1[v] - prev_grad[v];
	}
      }
      else // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	for (pt=0; pt<num_tp_pts; ++pt) {
	  c_index = (empty_c_index) ? cntr++ : colloc_index[lev][set][pt];
	  cov_t1c_ls[pt]
	    = (sdr_array_1[c_index].response_function() - mean_1)
	    * (sdr_array_2[c_index].response_function() - mean_2)
	    - value(sdv_array[c_index].continuous_variables(), sm_mi,
		    colloc_key, cov_t1c, cov_t2c, lev-1);
	}
    }
  }

#ifdef DEBUG
  PCout << "Central product interpolant type1 coeffs:\n" << cov_t1c;
  if (use_derivs)
    PCout << "Central product interpolant type2 coeffs:\n" << cov_t2c;
#endif // DEBUG
}


/** This version evaluates the r1 and r2 interpolants for cases where a
    particular key of the surrogate data does not correspond (e.g., after
    combine_coefficients).  Note: whereas expectation() supports either a
    reference or increment key (passed as generic set_partition), functions
    forming hierarchical product interpolants only support a reference key
    due to nonlinearity (only end point can be controlled).  Note: consider 
    computing deltaR1R2 = R1 deltaR2 + R2 deltaR1 + deltaR1 deltaR2. */
void HierarchInterpPolyApproximation::
central_product_interpolant(const RealMatrix2DArray& var_sets,
			    const UShort3DArray&     sm_mi,
			    const UShort4DArray&     colloc_key,
			    const RealVector2DArray& r1_t1c,
			    const RealMatrix2DArray& r1_t2c,
			    const RealVector2DArray& r2_t1c,
			    const RealMatrix2DArray& r2_t2c,
			    bool same, Real mean_r1, Real mean_r2,
			    RealVector2DArray& cov_t1c,
			    RealMatrix2DArray& cov_t2c,
			    const UShort2DArray& set_partition)
{
  // form hierarchical t1/t2 coeffs for (R_1 - \mu_1) (R_2 - \mu_2)
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_lev = r1_t1c.size(), set_start = 0, set_end,
    num_sets, num_tp_pts, v, num_v = sharedDataRep->numVars;
  bool partial = !set_partition.empty(),
    use_derivs = data_rep->basisConfigOptions.useDerivs;
  Real r1_val_mm1, r2_val_mm2;

  // level 0 (no surplus)
  cov_t1c.resize(num_lev);  cov_t2c.resize(num_lev);
  if (!partial || set_partition[0][0] == 0) { // partition includes level 0
    cov_t1c[0].resize(1);  cov_t2c[0].resize(1);
    cov_t1c[0][0].sizeUninitialized(1);
    RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		      (int)num_v);
    r1_val_mm1 = value(c_vars, sm_mi, colloc_key, r1_t1c, r1_t2c, 0) - mean_r1;
    r2_val_mm2 = (same) ? r1_val_mm1 :
      value(c_vars, sm_mi, colloc_key, r2_t1c, r2_t2c, 0) - mean_r2;
    cov_t1c[0][0][0] = r1_val_mm1 * r2_val_mm2;
    if (use_derivs) {
      const RealVector& r1_grad = gradient_basis_variables(c_vars, sm_mi,
        colloc_key, r1_t1c, r1_t2c, 0);
      const RealVector& r2_grad = (same) ? r1_grad : gradient_basis_variables(
        c_vars, sm_mi, colloc_key, r2_t1c, r2_t2c, 0);
      // Note: mean is only a function of nonbasis variables, so
      // basis vars grad of r{1,2}_val_mm{1,2} is only r{1,2}_grad
      cov_t2c[0][0].shapeUninitialized(num_v, 1);
      Real *cov_t2c_000 = cov_t2c[0][0][0];
      for (v=0; v<num_v; ++v)
	cov_t2c_000[v] = r1_val_mm1 * r2_grad[v] + r2_val_mm2 * r1_grad[v];
    }
  }

  // levels 1:w (subtract lev-1 for surplus estimation)
  for (lev=1; lev<num_lev; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];  num_sets = key_l.size();
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else
      set_end = num_sets;
    cov_t1c[lev].resize(num_sets); cov_t2c[lev].resize(num_sets);
    for (set=set_start; set<set_end; ++set) {
      num_tp_pts = key_l[set].size();
      RealVector& cov_t1c_ls = cov_t1c[lev][set];
      cov_t1c_ls.sizeUninitialized(num_tp_pts);
      if (use_derivs) {
	RealMatrix& cov_t2c_ls = cov_t2c[lev][set];
	cov_t2c_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt) {
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[lev][set][pt]), (int)num_v);
	  // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  r1_val_mm1 =
	    value(c_vars, sm_mi, colloc_key, r1_t1c, r1_t2c, lev) - mean_r1;
	  r2_val_mm2 = (same) ? r1_val_mm1 :
	    value(c_vars, sm_mi, colloc_key, r2_t1c, r2_t2c, lev) - mean_r2;
	  cov_t1c_ls[pt] = r1_val_mm1 * r2_val_mm2 -
	    value(c_vars, sm_mi, colloc_key, cov_t1c, cov_t2c, lev-1);
	  // type2 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  // --> interpolated grads are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
	  Real* cov_t2c_lsp = cov_t2c_ls[pt];
	  const RealVector& r1_grad = gradient_basis_variables(c_vars, sm_mi,
	    colloc_key, r1_t1c, r1_t2c, lev);
	  const RealVector& r2_grad = (same) ? r1_grad :
	    gradient_basis_variables(c_vars, sm_mi, colloc_key, r2_t1c,
	    r2_t2c, lev);
	  const RealVector& prev_grad = gradient_basis_variables(c_vars, sm_mi,
	    colloc_key, cov_t1c, cov_t2c, lev-1);
	  // Note: mean is only a function of nonbasis variables, so
	  // basis vars grad of r{1,2}_val_mm{1,2} is only r{1,2}_grad
	  for (v=0; v<num_v; ++v)
	    cov_t2c_lsp[v] = r1_val_mm1 * r2_grad[v] + r2_val_mm2 * r1_grad[v]
	                   - prev_grad[v];
	}
      }
      else // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	for (pt=0; pt<num_tp_pts; ++pt) {
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[lev][set][pt]), (int)num_v);
	  r1_val_mm1 =
	    value(c_vars, sm_mi, colloc_key, r1_t1c, r1_t2c, lev) - mean_r1;
	  r2_val_mm2 = (same) ? r1_val_mm1 :
	    value(c_vars, sm_mi, colloc_key, r2_t1c, r2_t2c, lev) - mean_r2;
	  cov_t1c_ls[pt] = r1_val_mm1 * r2_val_mm2 -
	    value(c_vars, sm_mi, colloc_key, cov_t1c, cov_t2c, lev-1);
	}
    }
  }

#ifdef DEBUG
  PCout << "Central product interpolant type1 coeffs:\n" << cov_t1c;
  if (use_derivs)
    PCout << "Central product interpolant type2 coeffs:\n" << cov_t2c;
#endif // DEBUG
}


void HierarchInterpPolyApproximation::
central_product_gradient_interpolant(const SDVArray& sdv_array,
				     const SDRArray& sdr_array_1,
				     const SDRArray& sdr_array_2, Real mean_r1,
				     Real mean_r2, const RealVector& mean1_grad,
				     const RealVector& mean2_grad,
				     const UShort3DArray& sm_mi,
				     const UShort4DArray& colloc_key,
				     const Sizet3DArray& colloc_index,
				     RealMatrix2DArray& cov_t1_coeff_grads,
				     const UShort2DArray& set_partition)
{
  // form hierarchical t1 coeff grads for (R_1 - \mu_1) (R_2 - \mu_2)
  size_t lev, set, pt, num_lev = colloc_key.size(), set_start = 0, set_end,
    num_sets, num_tp_pts, cntr = 0, c_index, v,
    num_deriv_vars = expT1CoeffGradsIter->second[0][0].numRows();
  bool partial = !set_partition.empty(), empty_c_index = colloc_index.empty();
  Real r1_mm, r2_mm;

  // level 0 (no surplus)
  cov_t1_coeff_grads.resize(num_lev);
  if (!partial || set_partition[0][0] == 0) { // partition includes level 0
    cov_t1_coeff_grads[0].resize(1);
    cov_t1_coeff_grads[0][0].shapeUninitialized(num_deriv_vars, 1);
    Real* cov_t1_coeff_grads_000 = cov_t1_coeff_grads[0][0][0];
    c_index = (empty_c_index) ? cntr : colloc_index[0][0][0];
    const SurrogateDataResp& sdr0_1 = sdr_array_1[c_index];
    const SurrogateDataResp& sdr0_2 = sdr_array_2[c_index];
    r1_mm = sdr0_1.response_function() - mean_r1;
    r2_mm = sdr0_2.response_function() - mean_r2;
    const RealVector& r1_grad = sdr0_1.response_gradient();
    const RealVector& r2_grad = sdr0_2.response_gradient();
    for (v=0; v<num_deriv_vars; ++v)
      cov_t1_coeff_grads_000[v] = r1_mm * (r2_grad[v] - mean2_grad[v])
	                        + r2_mm * (r1_grad[v] - mean1_grad[v]);
  }
  ++cntr;

  // levels 1:w (subtract lev-1 for surplus estimation)
  for (lev=1; lev<num_lev; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];  num_sets = key_l.size();
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else
      set_end = num_sets;
    cov_t1_coeff_grads[lev].resize(num_sets);
    if (empty_c_index)
      for (set=0; set<set_start; ++set)
	cntr += key_l[set].size();
    for (set=set_start; set<set_end; ++set) {
      num_tp_pts = key_l[set].size();
      RealMatrix& cov_t1_coeff_grads_ls = cov_t1_coeff_grads[lev][set];
      cov_t1_coeff_grads_ls.shapeUninitialized(num_deriv_vars, num_tp_pts);
      // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      for (pt=0; pt<num_tp_pts; ++pt) {
	c_index = (empty_c_index) ? cntr++ : colloc_index[lev][set][pt];
	const SurrogateDataResp& sdr_1 = sdr_array_1[c_index];
	const SurrogateDataResp& sdr_2 = sdr_array_2[c_index];
        r1_mm = sdr_1.response_function() - mean_r1;
	r2_mm = sdr_2.response_function() - mean_r2;
	const RealVector& r1_grad = sdr_1.response_gradient();
	const RealVector& r2_grad = sdr_2.response_gradient();
	const RealVector& prev_grad = gradient_nonbasis_variables(
	  sdv_array[c_index].continuous_variables(), sm_mi, colloc_key,
	  cov_t1_coeff_grads, lev-1);
	Real* cov_t1_coeff_grads_lsp = cov_t1_coeff_grads_ls[pt];
	for (v=0; v<num_deriv_vars; ++v)
	  cov_t1_coeff_grads_lsp[v] = r1_mm * (r2_grad[v] - mean2_grad[v])
	    + r2_mm * (r1_grad[v] - mean1_grad[v]) - prev_grad[v];
      }
    }
  }

#ifdef DEBUG
  PCout << "Central product gradient interpolant type1 coeff grads:\n"
	<< cov_t1_coeff_grads << std::endl;
#endif // DEBUG
}


void HierarchInterpPolyApproximation::
central_product_gradient_interpolant(const RealMatrix2DArray& var_sets,
				     const UShort3DArray&     sm_mi,
				     const UShort4DArray&     colloc_key,
				     const RealVector2DArray& r1_t1_coeffs,
				     const RealMatrix2DArray& r1_t2_coeffs,
				     const RealMatrix2DArray& r1_t1_coeff_grads,
				     const RealVector2DArray& r2_t1_coeffs,
				     const RealMatrix2DArray& r2_t2_coeffs,
				     const RealMatrix2DArray& r2_t1_coeff_grads,
				     bool same, Real mean_r1, Real mean_r2,
				     const RealVector& mean1_grad,
				     const RealVector& mean2_grad, 
				     RealMatrix2DArray& cov_t1_coeff_grads,
				     const UShort2DArray& set_partition)
{
  // form hierarchical t1 coeff grads for (R_1 - \mu_1) (R_2 - \mu_2)
  size_t lev, set, pt, num_lev = colloc_key.size(),
    num_sets, set_start = 0, set_end, num_tp_pts,
    v, num_deriv_vars = expT1CoeffGradsIter->second[0][0].numRows();
  bool partial = !set_partition.empty();
  Real r1_mm, r2_mm;

  // level 0 (no surplus)
  cov_t1_coeff_grads.resize(num_lev);
  if (!partial || set_partition[0][0] == 0) { // partition includes level 0
    cov_t1_coeff_grads[0].resize(1);
    cov_t1_coeff_grads[0][0].shapeUninitialized(num_deriv_vars, 1);
    Real* cov_t1_coeff_grads_000 = cov_t1_coeff_grads[0][0][0];
    RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		      (int)num_deriv_vars);
    r1_mm = value(c_vars, sm_mi, colloc_key, r1_t1_coeffs, r1_t2_coeffs, 0) -
      mean_r1;
    r2_mm = (same) ? r1_mm :
      value(c_vars, sm_mi, colloc_key, r2_t1_coeffs, r2_t2_coeffs, 0) - mean_r2;
    const RealVector& r1_grad = gradient_nonbasis_variables(c_vars, sm_mi,
      colloc_key, r1_t1_coeff_grads, 0);
    const RealVector& r2_grad = (same) ? r1_grad : gradient_nonbasis_variables(
      c_vars, sm_mi, colloc_key, r2_t1_coeff_grads, 0);
    for (v=0; v<num_deriv_vars; ++v)
      cov_t1_coeff_grads_000[v] = r1_mm * (r2_grad[v] - mean2_grad[v])
	                        + r2_mm * (r1_grad[v] - mean1_grad[v]);
  }

  // levels 1:w (subtract lev-1 for surplus estimation)
  for (lev=1; lev<num_lev; ++lev) {
    const UShort3DArray& key_l = colloc_key[lev];  num_sets = key_l.size();
    if (partial)
      { set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
    else
      set_end = num_sets;
    cov_t1_coeff_grads[lev].resize(num_sets);
    for (set=set_start; set<set_end; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      RealMatrix& cov_t1_coeff_grads_ls = cov_t1_coeff_grads[lev][set];
      cov_t1_coeff_grads_ls.shapeUninitialized(num_deriv_vars, num_tp_pts);
      // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      for (pt=0; pt<num_tp_pts; ++pt) {
	RealVector c_vars(Teuchos::View,
	  const_cast<Real*>(var_sets[lev][set][pt]), (int)num_deriv_vars);
        r1_mm = value(c_vars, sm_mi, colloc_key, r1_t1_coeffs, r1_t2_coeffs,
	  lev) - mean_r1;
	r2_mm = (same) ? r1_mm : value(c_vars, sm_mi, colloc_key, r2_t1_coeffs,
	  r2_t2_coeffs, lev) - mean_r2;
	const RealVector& r1_grad = gradient_nonbasis_variables(c_vars, sm_mi,
	  colloc_key, r1_t1_coeff_grads, lev);
	const RealVector& r2_grad = gradient_nonbasis_variables(c_vars, sm_mi,
	  colloc_key, r2_t1_coeff_grads, lev);
	const RealVector& prev_grad = gradient_nonbasis_variables(c_vars, sm_mi,
	  colloc_key, cov_t1_coeff_grads, lev-1);
	Real* cov_t1_coeff_grads_lsp = cov_t1_coeff_grads_ls[pt];
	for (v=0; v<num_deriv_vars; ++v)
	  cov_t1_coeff_grads_lsp[v]
	    = r1_mm * (r2_grad[v] - mean2_grad[v])
	    + r2_mm * (r1_grad[v] - mean1_grad[v]) - prev_grad[v];
      }
    }
  }

#ifdef DEBUG
  PCout << "Central product gradient interpolant type1 coeff grads:\n"
	<< cov_t1_coeff_grads << std::endl;
#endif // DEBUG
}


void HierarchInterpPolyApproximation::
central_product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
  Real mean_1, Real mean_2,
  std::map<UShortArray, RealVector2DArray>& cov_t1c_map,
  std::map<UShortArray, RealMatrix2DArray>& cov_t2c_map,
  const std::map<UShortArray, UShort2DArray>& set_partition_map)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();

  std::map<UShortArray, RealVector2DArray>::const_iterator t1c_cit1, t1c_cit2;
  std::map<UShortArray, RealMatrix2DArray>::const_iterator
    t2c_cit1, t2c_cit2, v_cit;
  std::map<UShortArray, UShort3DArray>::const_iterator sm_cit;
  std::map<UShortArray, UShort4DArray>::const_iterator ck_cit;
  std::map<UShortArray,  Sizet3DArray>::const_iterator ci_cit;
  std::map<UShortArray, UShort2DArray>::const_iterator p_cit;
  std::map<UShortArray, SDVArray>::const_iterator sdv_cit;
  std::map<UShortArray, SDRArray>::const_iterator sdr1_cit, sdr2_cit;

  bool same = (this == hip_approx_2),
    track_c_index = hsg_driver->track_collocation_indices();

  for (t1c_cit1  = expansionType1Coeffs.begin(),
       t1c_cit2  = hip_approx_2->expansionType1Coeffs.begin(),
       t2c_cit1  = expansionType2Coeffs.begin(),
       t2c_cit2  = hip_approx_2->expansionType2Coeffs.begin(),
       sdv_cit   = modSurrData.variables_data_map().begin(),
       sdr1_cit  = modSurrData.response_data_map().begin(),
       sdr2_cit  = hip_approx_2->modSurrData.response_data_map().begin(),
       v_cit     = hsg_driver->variable_sets_map().begin(),
       sm_cit    = hsg_driver->smolyak_multi_index_map().begin(),
       ck_cit    = hsg_driver->collocation_key_map().begin(),
       ci_cit    = hsg_driver->collocation_indices_map().begin(),
       p_cit     = set_partition_map.begin();
       t1c_cit1 != expansionType1Coeffs.end();
       ++t1c_cit1, ++t1c_cit2, ++t2c_cit1, ++t2c_cit2, ++sdv_cit, ++sdr1_cit,
       ++sdr2_cit, ++v_cit, ++sm_cit, ++ck_cit, ++ci_cit, ++p_cit) {
    const UShortArray& key     = t1c_cit1->first;
    RealVector2DArray& cov_t1c = cov_t1c_map[key]; // update in place
    RealMatrix2DArray& cov_t2c = cov_t2c_map[key]; // update in place

    if (track_c_index && ci_cit->second.empty())// active invalidated by combine
      central_product_interpolant(v_cit->second, sm_cit->second, ck_cit->second,
	t1c_cit1->second, t2c_cit1->second, t1c_cit2->second, t2c_cit2->second,
	same, mean_1, mean_2, cov_t1c, cov_t2c, p_cit->second);
    else // use modSurrData & colloc_indices for forming central product interp
      central_product_interpolant(sdv_cit->second, sdr1_cit->second,
	sdr2_cit->second, mean_1, mean_2, sm_cit->second, ck_cit->second,
	ci_cit->second, cov_t1c, cov_t2c, p_cit->second);
  }
}


/** See integrate_response_moments(size_t, bool) -- this variant is only
    valid for special case of active coefficients/stats when surrogate
    data is available (collocation indices not invalidated). */
void HierarchInterpPolyApproximation::
integrate_response_moments(size_t num_moments, const UShort3DArray& sm_mi,
			   const UShort4DArray& colloc_key,
			   const Sizet3DArray&  colloc_index,
			   const SDVArray& sdv_array, const SDRArray& sdr_array)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::integrate_response_moments()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_lev = colloc_key.size(), num_sets, num_tp_pts,
    cntr, c_index, v, num_v = sharedDataRep->numVars;
  int m_index, moment;
  bool empty_c_index = colloc_index.empty(),
    use_derivs = data_rep->basisConfigOptions.useDerivs;

  RealVector& numer_mom = numMomentsIter->second;
  if (numer_mom.length() != num_moments)
    numer_mom.sizeUninitialized(num_moments);
  Real& mean = numer_mom[0];
  // this integrate_response_moments() variant supports active exp coeffs:
  mean = expectation(expT1CoeffsIter->second, expT2CoeffsIter->second);

  // size moment coefficient arrays
  RealVector2DArray mom_t1_coeffs(num_lev);
  RealMatrix2DArray mom_t2_coeffs(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    num_sets = colloc_key[lev].size();
    mom_t1_coeffs[lev].resize(num_sets);
    mom_t2_coeffs[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      mom_t1_coeffs[lev][set].sizeUninitialized(num_tp_pts);
      if (use_derivs)
	mom_t2_coeffs[lev][set].shapeUninitialized(num_v, num_tp_pts);
    }
  }

  for (m_index=1; m_index<num_moments; ++m_index) {
    moment = m_index+1; cntr = 0;
    if (use_derivs) {
      // level 0
      c_index = (empty_c_index) ? cntr++ : colloc_index[0][0][0];
      const SurrogateDataResp& sdr0 = sdr_array[c_index];
      Real data_fn_mm         = sdr0.response_function() - mean;
      mom_t1_coeffs[0][0][0]  = std::pow(data_fn_mm, moment);
      Real* mom_t2_coeffs_000 = mom_t2_coeffs[0][0][0];
      Real deriv = moment * std::pow(data_fn_mm, m_index);
      const RealVector& data_grad = sdr0.response_gradient();
      for (v=0; v<num_v; ++v)
	mom_t2_coeffs_000[v] = deriv * data_grad[v];
      // levels 1:w
      for (lev=1; lev<num_lev; ++lev) {
	num_sets = colloc_key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  num_tp_pts = colloc_key[lev][set].size();
	  RealVector& mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  RealMatrix& mom_t2_coeffs_ls = mom_t2_coeffs[lev][set];
	  for (pt=0; pt<num_tp_pts; ++pt) {
	    c_index = (empty_c_index) ? cntr++ : colloc_index[lev][set][pt];
	    const RealVector& c_vars
	      = sdv_array[c_index].continuous_variables();
	    const SurrogateDataResp& sdr = sdr_array[c_index];
	    // type1 hierarchical interpolation of (R - \mu)^moment
	    data_fn_mm = sdr.response_function() - mean;
	    mom_t1_coeffs_ls[pt] = std::pow(data_fn_mm, moment)
	      - value(c_vars, sm_mi, colloc_key, mom_t1_coeffs,
		      mom_t2_coeffs, lev-1);
	    // type2 hierarchical interpolation of (R - \mu)^moment
	    // --> interpolated grads are moment(R-\mu)^{moment-1} R'
	    Real* mom_t2_coeffs_lsp = mom_t2_coeffs_ls[pt];
	    deriv = moment * std::pow(data_fn_mm, m_index);
	    const RealVector& data_grad = sdr.response_gradient();
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, colloc_key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    for (v=0; v<num_v; ++v)
	      mom_t2_coeffs_lsp[v] = deriv * data_grad[v] - prev_grad[v];
	  }
	}
      }
    }
    else {
      // level 0
      c_index = (empty_c_index) ? cntr++ : colloc_index[0][0][0];
      mom_t1_coeffs[0][0][0]
	= std::pow(sdr_array[c_index].response_function() - mean, moment);
      // levels 1:w
      for (lev=1; lev<num_lev; ++lev) {
	num_sets = colloc_key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  RealVector& mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  num_tp_pts = colloc_key[lev][set].size();
	  // type1 hierarchical interpolation of (R - \mu)^moment
	  for (pt=0; pt<num_tp_pts; ++pt) {
	    c_index = (empty_c_index) ? cntr++ : colloc_index[lev][set][pt];
	    mom_t1_coeffs_ls[pt]
	      = std::pow(sdr_array[c_index].response_function() - mean, moment)
	      - value(sdv_array[c_index].continuous_variables(), sm_mi,
		      colloc_key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	  }
	}
      }
    }

    // pass these exp coefficients into a general expectation fn
    numer_mom[m_index] = expectation(mom_t1_coeffs, mom_t2_coeffs);
  }

  // standardize third and higher central moments, if present
  //standardize_moments(numer_mom);

  /*
  if (numer_mom.size() != num_moments)
    numer_mom.size(num_moments);
  if (use_derivs)
    integrate_moments(expT1CoeffsIter->second, expT2CoeffsIter->second,
		      hsg_driver->type1_hierarchical_weight_sets(),
		      hsg_driver->type2_hierarchical_weight_sets(), numer_mom);
  else
    integrate_moments(expT1CoeffsIter->second,
                      hsg_driver->type1_hierarchical_weight_sets(), numer_mom);
  */
}


void HierarchInterpPolyApproximation::
integrate_response_moments(size_t num_moments,
			   const RealMatrix2DArray& var_sets,
			   const UShort3DArray&     sm_mi,
			   const UShort4DArray&     colloc_key,
			   const RealVector2DArray& t1_coeffs,
			   const RealMatrix2DArray& t2_coeffs,
			   const RealVector2DArray& t1_wts,
			   const RealMatrix2DArray& t2_wts)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::integrate_response_moments()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  size_t lev, set, pt, num_lev = colloc_key.size(), num_sets, num_tp_pts,
    v, num_v = sharedDataRep->numVars;
  int m_index, moment;
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;

  RealVector& numer_mom = numMomentsIter->second;
  if (numer_mom.length() != num_moments)
    numer_mom.sizeUninitialized(num_moments);
  Real& mean = numer_mom[0];
  mean = expectation(t1_coeffs, t2_coeffs, t1_wts, t2_wts);

  // size moment coefficient arrays
  RealVector2DArray mom_t1_coeffs(num_lev),  r_vals(num_lev);
  RealMatrix2DArray mom_t2_coeffs(num_lev), r_grads(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    num_sets = colloc_key[lev].size();
    mom_t1_coeffs[lev].resize(num_sets);  r_vals[lev].resize(num_sets);
    mom_t2_coeffs[lev].resize(num_sets);  r_grads[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      mom_t1_coeffs[lev][set].sizeUninitialized(num_tp_pts);
      r_vals[lev][set].sizeUninitialized(num_tp_pts);
      if (use_derivs) {
	mom_t2_coeffs[lev][set].shapeUninitialized(num_v, num_tp_pts);
	r_grads[lev][set].shapeUninitialized(num_v, num_tp_pts);
      }
    }
  }

  // precompute all value() and gradient_basis_variables() for {t1,t2}_coeffs
  // Note: this is essentially what synthetic_surrogate_data() does, following
  //       combined_to_active() (eliminating the need to exercise this
  //       overloaded form of integrate_response_moments() for FINAL_RESULTS).
  RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		    (int)num_v);
  r_vals[0][0][0] =
    value(c_vars, sm_mi, colloc_key, t1_coeffs, t2_coeffs, 0) - mean;
  if (use_derivs)
    Teuchos::setCol(gradient_basis_variables(c_vars, sm_mi, colloc_key,
      t1_coeffs, t2_coeffs, 0), 0, r_grads[0][0]);
  for (lev=1; lev<num_lev; ++lev) {
    num_sets = colloc_key[lev].size();
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = colloc_key[lev][set].size();
      RealVector&         r_vals_ls =   r_vals[lev][set];
      RealMatrix&        r_grads_ls =  r_grads[lev][set];
      const RealMatrix& var_sets_ls = var_sets[lev][set];
      for (pt=0; pt<num_tp_pts; ++pt) {
	RealVector
	  c_vars(Teuchos::View, const_cast<Real*>(var_sets_ls[pt]), (int)num_v);
	r_vals_ls[pt] =
	  value(c_vars, sm_mi, colloc_key, t1_coeffs, t2_coeffs, lev) - mean;
	if (use_derivs)
	  Teuchos::setCol(gradient_basis_variables(c_vars, sm_mi, colloc_key,
	    t1_coeffs, t2_coeffs, lev), (int)pt, r_grads_ls);
      }
    }
  }

  for (m_index=1; m_index<num_moments; ++m_index) {
    moment = m_index+1;

    // level 0
    RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[0][0][0]),
		      (int)num_v);
    Real r_val_mm = r_vals[0][0][0];
    mom_t1_coeffs[0][0][0] = std::pow(r_val_mm, moment);
    if (use_derivs) {
      // level 0 continued
      Real *mom_t2_coeffs_000 = mom_t2_coeffs[0][0][0],
	   *r_grad = r_grads[0][0][0],
	   deriv = moment * std::pow(r_val_mm, m_index);
      for (v=0; v<num_v; ++v)
	mom_t2_coeffs_000[v] = deriv * r_grad[v];
    
      // levels 1:w
      for (lev=1; lev<num_lev; ++lev) {
	num_sets = colloc_key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  num_tp_pts = colloc_key[lev][set].size();
	  RealVector&  mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  RealMatrix&  mom_t2_coeffs_ls = mom_t2_coeffs[lev][set];
	  RealVector&         r_vals_ls =        r_vals[lev][set];
	  RealMatrix&        r_grads_ls =       r_grads[lev][set];
	  const RealMatrix& var_sets_ls =      var_sets[lev][set];
	  for (pt=0; pt<num_tp_pts; ++pt) {
	    RealVector c_vars(Teuchos::View,
	      const_cast<Real*>(var_sets_ls[pt]), (int)num_v);
	    // type1 hierarchical interpolation of (R - \mu)^moment
	    r_val_mm = r_vals_ls[pt];
	    mom_t1_coeffs_ls[pt] = std::pow(r_val_mm, moment) - value(c_vars,
	      sm_mi, colloc_key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    // type2 hierarchical interpolation of (R - \mu)^moment
	    // --> interpolated grads are moment(R-\mu)^{moment-1} R'
	    Real *mom_t2_coeffs_lsp = mom_t2_coeffs_ls[pt],
	         *r_grad = r_grads_ls[pt],
	         deriv = moment * std::pow(r_val_mm, m_index);
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, colloc_key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    for (v=0; v<num_v; ++v)
	      mom_t2_coeffs_lsp[v] = deriv * r_grad[v] - prev_grad[v];
	  }
	}
      }
    }
    else {
      // levels 1:w
      for (lev=1; lev<num_lev; ++lev) {
	num_sets = colloc_key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  RealVector&  mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  RealVector&         r_vals_ls =        r_vals[lev][set];
	  const RealMatrix& var_sets_ls =      var_sets[lev][set];
	  num_tp_pts = colloc_key[lev][set].size();
	  // type1 hierarchical interpolation of (R - \mu)^moment
	  for (pt=0; pt<num_tp_pts; ++pt) {
	    RealVector c_vars(Teuchos::View,
	      const_cast<Real*>(var_sets_ls[pt]), (int)num_v);
	    mom_t1_coeffs_ls[pt] = std::pow(r_vals_ls[pt], moment) - value(
	      c_vars, sm_mi, colloc_key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	  }
	}
      }
    }

    // pass these exp coefficients into a general expectation fn
    numer_mom[m_index]
      = expectation(mom_t1_coeffs, mom_t2_coeffs, t1_wts, t2_wts);
  }

  // standardize third and higher central moments, if present
  //standardize_moments(numer_mom);
}


/** Computes the variance of component functions.  Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value. */
void HierarchInterpPolyApproximation::
compute_partial_variance(const BitArray& set_value)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  // Compute the partial integral corresponding to set_value
  Real& variance = partialVariance[data_rep->sobolIndexMap[set_value]];

  // Follows Tang, Iaccarino, Eldred (conference paper AIAA-2010-2922),
  // but adapted for hierarchical interpolation.

  // Perform inner integral over complementary set u' to form new weighted
  // coefficients h (stored as member_coeffs), stored hierarchically over
  // the member dimensions.
  RealVector2DArray member_t1_coeffs, member_t1_wts;
  RealMatrix2DArray member_t2_coeffs, member_t2_wts;
  UShort4DArray member_colloc_key; Sizet3DArray member_colloc_index;
  member_coefficients_weights(set_value, member_t1_coeffs, member_t1_wts,
			      member_t2_coeffs, member_t2_wts,
			      member_colloc_key, member_colloc_index);

  // re-interpolate h^2 (0 passed for mean to central product) using a
  // hierarchical formulation over the reduced member dimensions
  RealVector2DArray prod_member_t1_coeffs;
  RealMatrix2DArray prod_member_t2_coeffs;
  central_product_member_coefficients(set_value, member_t1_coeffs,
				      member_t2_coeffs, member_colloc_key,
				      member_colloc_index, 0.,
				      prod_member_t1_coeffs,
				      prod_member_t2_coeffs);

  // integrate h^2 over the reduced member dimensions
  variance = expectation(prod_member_t1_coeffs, prod_member_t2_coeffs,
			 member_t1_wts,         member_t2_wts);

#ifdef VBD_DEBUG
  PCout << "Partial variance = " << variance;
#endif // VBD_DEBUG
  // compute proper subsets and subtract their contributions
  InterpPolyApproximation::compute_partial_variance(set_value);
#ifdef VBD_DEBUG
  PCout << " (raw) " << variance << " (minus subsets)\n";
#endif // VBD_DEBUG
}


void HierarchInterpPolyApproximation::compute_total_sobol_indices()
{
  // Compute the total expansion mean and variance.  For standard mode, these
  // are likely already available, as managed by computed{Mean,Variance} data
  // reuse trackers.  For all vars mode, they are computed without partial
  // integration restricted to the random indices.  If negligible variance
  // (deterministic test fn) or negative variance (poor sparse grid resolution),
  // then attribution of this variance is suspect.  Note: zero is a good choice
  // since it drops out from anisotropic refinement based on a response-average
  // of total Sobol' indices.
  Real total_variance = variance();
  if (total_variance <= SMALL_NUMBER)
    { totalSobolIndices = 0.; return; }
  Real total_mean = mean();

  UShortArray quad_order; Real complement_variance;
  size_t v, num_v = sharedDataRep->numVars;
  BitArray complement_set(num_v);
  RealVector2DArray member_t1_coeffs, member_t1_wts, cprod_member_t1_coeffs;
  RealMatrix2DArray member_t2_coeffs, member_t2_wts, cprod_member_t2_coeffs;
  UShort4DArray member_colloc_key; Sizet3DArray member_colloc_index;
  // iterate each variable 
  for (v=0; v<num_v; ++v) {
    // define complement_set that includes all but index of interest
    complement_set.set(); complement_set.flip(v);

    // Perform inner integral over complementary set u' to form new weighted
    // coefficients h (stored as member_coeffs), stored hierarchically over
    // the member dimensions.
    member_coefficients_weights(complement_set, member_t1_coeffs, member_t1_wts,
				member_t2_coeffs, member_t2_wts,
				member_colloc_key, member_colloc_index);

    // re-interpolate (h-\mu)^2 using a hierarchical formulation over the
    // reduced member dimensions
    central_product_member_coefficients(complement_set, member_t1_coeffs,
					member_t2_coeffs, member_colloc_key,
					member_colloc_index, total_mean,
					cprod_member_t1_coeffs,
					cprod_member_t2_coeffs);

    // integrate (h-\mu)^2 over the reduced member dimensions using member wts
    complement_variance = expectation(cprod_member_t1_coeffs,
				      cprod_member_t2_coeffs, member_t1_wts,
				      member_t2_wts);

    // define total Sobol' index for this var from complement & total variance
    totalSobolIndices[v] = 1. - complement_variance / total_variance;
  }
}


/** Forms a lower dimensional interpolant over variables that are
    members of the given set. */
void HierarchInterpPolyApproximation::
member_coefficients_weights(const BitArray& member_bits,
  RealVector2DArray& member_t1_coeffs,  RealVector2DArray& member_t1_wts,
  RealMatrix2DArray& member_t2_coeffs,  RealMatrix2DArray& member_t2_wts,
  UShort4DArray&     member_colloc_key, Sizet3DArray&      member_colloc_index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key   = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();

  const RealVector2DArray& exp_t1_coeffs = expT1CoeffsIter->second;
  const RealMatrix2DArray& exp_t2_coeffs = expT2CoeffsIter->second;

  // Perform inner integral over complementary set u' (non-member vars) to
  // form new weighted expansion coefficients h (stored as member_t1_coeffs)
  size_t v, num_v = sharedDataRep->numVars, lev, set, pt,
    num_lev = exp_t1_coeffs.size(), num_sets, num_tp_pts,
    num_member_coeffs, v_cntr, p_cntr = 0, m_index;
  bool empty_c_index = colloc_index.empty(),
    use_derivs = data_rep->basisConfigOptions.useDerivs;
  Real member_wt, nonmember_wt;
  SizetArray indexing_factor; // for indexing member coeffs,wts
  UShortArray delta_sizes; UShort2DArray delta_keys;
  member_t1_coeffs.resize(num_lev);  member_t1_wts.resize(num_lev);
  member_t2_coeffs.resize(num_lev);  member_t2_wts.resize(num_lev);
  member_colloc_key.resize(num_lev); member_colloc_index.resize(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    const RealVectorArray& t1_coeffs_l = exp_t1_coeffs[lev];
    num_sets = t1_coeffs_l.size();
    member_t1_coeffs[lev].resize(num_sets); member_t1_wts[lev].resize(num_sets);
    member_t2_coeffs[lev].resize(num_sets); member_t2_wts[lev].resize(num_sets);
    member_colloc_key[lev].resize(num_sets);
    member_colloc_index[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
      const RealMatrix& t2_coeffs_ls = exp_t2_coeffs[lev][set];
      const UShortArray&    sm_mi_ls = sm_mi[lev][set];
      RealVector&     m_t1_coeffs_ls = member_t1_coeffs[lev][set];
      RealVector&        m_t1_wts_ls = member_t1_wts[lev][set];
      RealMatrix&     m_t2_coeffs_ls = member_t2_coeffs[lev][set];
      RealMatrix&        m_t2_wts_ls = member_t2_wts[lev][set];
      UShort2DArray&        m_key_ls = member_colloc_key[lev][set];
      SizetArray&         m_index_ls = member_colloc_index[lev][set];

      // precompute sizes and indexing factors and init member arrays
      num_tp_pts = t1_coeffs_ls.length();
      if (num_tp_pts) { // can be zero in case of growth restriction
	hsg_driver->levels_to_delta_keys(sm_mi_ls, delta_keys);
	hsg_driver->levels_to_delta_sizes(sm_mi_ls, delta_sizes);
	indexing_factor.clear();
	for (v=0, num_member_coeffs=1; v<num_v; ++v)
	  if (member_bits[v]) {
	    indexing_factor.push_back(num_member_coeffs); // for m_index below
	    num_member_coeffs *= delta_sizes[v];
	  }
	m_t1_coeffs_ls.size(num_member_coeffs);
	m_t1_wts_ls.size(num_member_coeffs);
	if (use_derivs) {
	  m_t2_coeffs_ls.shape(num_v, num_member_coeffs);
	  m_t2_wts_ls.shape(num_v, num_member_coeffs);
	}
	m_key_ls.resize(num_member_coeffs);
	m_index_ls.resize(num_member_coeffs);
      }

      for (pt=0; pt<num_tp_pts; ++pt) {
	// convert key_lsp to a corresponding index on member_t1_{coeffs,wts}.
	// We must also define a mapping to a member-variable collocation
	// key/index (from collapsing non-members) for use downstream.
	const UShortArray& key_lsp = colloc_key[lev][set][pt];
	for (v=0, m_index=0, v_cntr=0; v<num_v; ++v)
	  if (member_bits[v]) // key_lsp spans all pts, not deltas
	    m_index += find_index(delta_keys[v], key_lsp[v])
	            *  indexing_factor[v_cntr++];

	// integrate out nonmember dimensions and aggregate with type1 coeffs
	data_rep->type1_weight(key_lsp, sm_mi_ls, member_bits,
			       member_wt, nonmember_wt);
	m_t1_coeffs_ls[m_index] += nonmember_wt * t1_coeffs_ls[pt];
	// reduced dimension data updated more times than necessary, but
	// tracking this redundancy would be more expensive/complicated.
	// Note: key and index->c_vars components are the same only for the
	// member dimensions later used in value()/gradient_basis_variables().
	m_t1_wts_ls[m_index] = member_wt;
	m_index_ls[m_index]  = (empty_c_index) ? p_cntr++ :
	  colloc_index[lev][set][pt];   // links back to modSurrData c_vars
	m_key_ls[m_index]    = key_lsp; // links back to interp polynomials

	// now do the same for the type2 coeffs and weights
	if (use_derivs) {
	  Real     *m_t2_coeffs_lsp = m_t2_coeffs_ls[m_index],
	              *m_t2_wts_lsp = m_t2_wts_ls[m_index];
	  const Real *t2_coeffs_lsp = t2_coeffs_ls[pt];
	  for (v=0; v<num_v; ++v) {
	    data_rep->type2_weight(v, key_lsp, sm_mi_ls, member_bits,
				   member_wt, nonmember_wt);
	    m_t2_coeffs_lsp[v] += nonmember_wt * t2_coeffs_lsp[v];
	    m_t2_wts_lsp[v]    =  member_wt;
	  }
	}
      }
#ifdef VBD_DEBUG
      PCout << "member_bits: " << member_bits // MSB->LSB: order reversed
	    << "\nmember_t1_coeffs[" << lev << "][" << set << "]:\n"
	    << member_t1_coeffs[lev][set] << "member_t1_wts[" << lev << "]["
	    << set << "]:\n" << member_t1_wts[lev][set];
      if (use_derivs) {
	PCout << "member_t2_coeffs[" << lev << "][" << set << "]:\n";
	write_data(PCout, member_t2_coeffs[lev][set], false, true, true);
	PCout << "member_t2_wts[" << lev << "][" << set << "]:\n";
	write_data(PCout, member_t2_wts[lev][set],    false, true, true);
      }
      PCout << std::endl;
#endif // VBD_DEBUG
    }
  }
}


void HierarchInterpPolyApproximation::
central_product_member_coefficients(const BitArray& m_bits,
  const RealVector2DArray& m_t1_coeffs, const RealMatrix2DArray& m_t2_coeffs,
  const UShort4DArray& m_colloc_key,    const Sizet3DArray& m_colloc_index,
  Real mean, RealVector2DArray& cprod_m_t1_coeffs,
  RealMatrix2DArray& cprod_m_t2_coeffs)
{
  // We now have a lower dimensional hierarchical interpolant, for which
  // (h - mean)^2 must be formed hierarchically.  We use value() to both
  // reconstruct h from its surpluses and to evaluate new surpluses for
  // (h - mean)^2.  Note that we can simply pass mean = 0 for an h^2
  // product interpolant (used by compute_partial_variance()).

  // while colloc_{key,index} are redefined for member vars, sm_mi is not
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShort3DArray& sm_mi = data_rep->hsg_driver()->smolyak_multi_index();
  const SDVArray& sdv_array = modSurrData.variables_data();
  bool empty_c_index = m_colloc_index.empty();

  size_t v, num_v = sharedDataRep->numVars, lev, set, pt, cntr = 0,
    num_lev = m_t1_coeffs.size(), num_sets, num_tp_pts, c_index,
    h_level = num_lev - 1;
  SizetList member_indices;
  for (v=0; v<num_v; ++v)
    if (m_bits[v])
      member_indices.push_back(v);

  // form hierarchical t1/t2 coeffs for (h - mean)^2
  cprod_m_t1_coeffs.resize(num_lev); cprod_m_t1_coeffs[0].resize(1);
  cprod_m_t2_coeffs.resize(num_lev); cprod_m_t2_coeffs[0].resize(1);
  // level 0 type1
  cprod_m_t1_coeffs[0][0].sizeUninitialized(1);
  c_index = (empty_c_index) ? cntr++ : m_colloc_index[0][0][0];
  const RealVector& c_vars0 = sdv_array[c_index].continuous_variables();
  Real h_val_mm = value(c_vars0, sm_mi, m_colloc_key, m_t1_coeffs, m_t2_coeffs,
			h_level, member_indices) - mean;
  cprod_m_t1_coeffs[0][0][0] = h_val_mm * h_val_mm;
  if (data_rep->basisConfigOptions.useDerivs) {
    // level 0 type2
    const RealVector& h_grad_000 = gradient_basis_variables(c_vars0, sm_mi,
      m_colloc_key, m_t1_coeffs, m_t2_coeffs, h_level, member_indices);
    cprod_m_t2_coeffs[0][0].shapeUninitialized(num_v, 1);
    Real *cprod_m_t2_coeffs_000 = cprod_m_t2_coeffs[0][0][0];
    for (v=0; v<num_v; ++v)
      cprod_m_t2_coeffs_000[v] = 2. * h_val_mm * h_grad_000[v];
    // levels 1:w type1 and type2
    for (lev=1; lev<num_lev; ++lev) {
      num_sets = m_colloc_key[lev].size();
      cprod_m_t1_coeffs[lev].resize(num_sets);
      cprod_m_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& cprod_m_t1_coeffs_ls = cprod_m_t1_coeffs[lev][set];
	RealMatrix& cprod_m_t2_coeffs_ls = cprod_m_t2_coeffs[lev][set];
	num_tp_pts = m_colloc_key[lev][set].size();
	cprod_m_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	cprod_m_t2_coeffs_ls.shapeUninitialized(num_v,num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt) {
	  c_index = (empty_c_index) ? cntr++ : m_colloc_index[lev][set][pt];
	  const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	  // type1 hierarchical interpolation of h^2
	  h_val_mm = value(c_vars, sm_mi, m_colloc_key, m_t1_coeffs,
			   m_t2_coeffs, h_level, member_indices) - mean;
	  cprod_m_t1_coeffs_ls[pt] = h_val_mm * h_val_mm -
	    value(c_vars, sm_mi, m_colloc_key, cprod_m_t1_coeffs,
		  cprod_m_t2_coeffs, lev-1, member_indices);
	  // type2 hierarchical interpolation of h^2
	  // --> interpolated grads are 2 h h'
	  // --> make a copy of h_grad since prev_grad reuses approxGradient
	  Real* cprod_m_t2_coeffs_lsp = cprod_m_t2_coeffs_ls[pt];
	  RealVector h_grad = gradient_basis_variables(c_vars, sm_mi,
	    m_colloc_key, m_t1_coeffs, m_t2_coeffs, h_level, member_indices);
	  const RealVector& prev_grad =
	    gradient_basis_variables(c_vars, sm_mi, m_colloc_key,
	    cprod_m_t1_coeffs, cprod_m_t2_coeffs, lev-1, member_indices);
	  for (v=0; v<num_v; ++v)
	    cprod_m_t2_coeffs_lsp[v] = 2. * h_val_mm * h_grad[v] - prev_grad[v];
	}
#ifdef VBD_DEBUG
	PCout << "cprod_m_t1_coeffs[" << lev << "][" << set << "]:\n"
	      << cprod_m_t1_coeffs[lev][set]
	      << "cprod_m_t2_coeffs[" << lev << "][" << set << "]:\n";
	write_data(PCout, cprod_m_t2_coeffs[lev][set], false, true, true);
	PCout << std::endl;
#endif // VBD_DEBUG
      }
    }
  }
  else {
    // levels 1:w type1
    for (lev=1; lev<num_lev; ++lev) {
      num_sets = m_colloc_key[lev].size();
      cprod_m_t1_coeffs[lev].resize(num_sets);
      cprod_m_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = m_colloc_key[lev][set].size();
	RealVector& cprod_m_t1_coeffs_ls = cprod_m_t1_coeffs[lev][set];
	cprod_m_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of (h - mean)^2
	for (pt=0; pt<num_tp_pts; ++pt) {
	  c_index = (empty_c_index) ? cntr++ : m_colloc_index[lev][set][pt];
	  const RealVector& c_vars = sdv_array[c_index].continuous_variables();
	  h_val_mm = value(c_vars, sm_mi, m_colloc_key, m_t1_coeffs,
			   m_t2_coeffs, h_level, member_indices) - mean;
	  cprod_m_t1_coeffs_ls[pt] = h_val_mm * h_val_mm -
	    value(c_vars, sm_mi, m_colloc_key, cprod_m_t1_coeffs,
		  cprod_m_t2_coeffs, lev-1, member_indices);
	}
#ifdef VBD_DEBUG
	PCout << "cprod_m_t1_coeffs[" << lev << "][" << set << "]:\n"
	      << cprod_m_t1_coeffs[lev][set] << std::endl;
#endif // VBD_DEBUG
      }
    }
  }
}

}
