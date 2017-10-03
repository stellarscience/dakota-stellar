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

namespace Pecos {


void HierarchInterpPolyApproximation::allocate_arrays()
{
  InterpPolyApproximation::allocate_arrays();

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const ExpansionConfigOptions& ec_options = data_rep->expConfigOptions;
  const UShort4DArray& key = data_rep->hsg_driver()->collocation_key();
  size_t i, j, k, num_levels = key.size(), num_sets, num_tp_pts,
    num_deriv_vars = surrData.num_derivative_variables();

  if (expansionType1Coeffs.size() != num_levels)
    expansionType1Coeffs.resize(num_levels);
  if ( expansionType2Coeffs.size() != num_levels)
    expansionType2Coeffs.resize(num_levels);
  if (expansionType1CoeffGrads.size() != num_levels)
    expansionType1CoeffGrads.resize(num_levels);
  for (i=0; i<num_levels; ++i) {
    const UShort3DArray& key_i = key[i];
    num_sets = key_i.size();
    if (expansionType1Coeffs[i].size() != num_sets)
      expansionType1Coeffs[i].resize(num_sets);
    if (expansionType2Coeffs[i].size() != num_sets)
      expansionType2Coeffs[i].resize(num_sets);
    if (expansionType1CoeffGrads[i].size() != num_sets)
      expansionType1CoeffGrads[i].resize(num_sets);
    for (j=0; j<num_sets; ++j) {
      num_tp_pts = key_i[j].size();
      for (k=0; k<num_tp_pts; ++k) {
	if (expansionCoeffFlag) {
	  expansionType1Coeffs[i][j].sizeUninitialized(num_tp_pts);
	  if (data_rep->basisConfigOptions.useDerivs)
	    expansionType2Coeffs[i][j].shapeUninitialized(num_deriv_vars,
							  num_tp_pts);
	}
	if (expansionCoeffGradFlag)
	  expansionType1CoeffGrads[i][j].shapeUninitialized(num_deriv_vars,
							    num_tp_pts);
      }
    }
  }

  // checking num_points is insufficient due to anisotropy --> changes in
  // anisotropic weights could move points around without changing the total
  //size_t num_points = surrData.points();
  //bool update_exp_form =
  //  ( (expansionCoeffFlag && expansionType1Coeffs.length() != num_points) ||
  //    (expansionCoeffGradFlag &&
  //     expansionType1CoeffGrads.numCols() != num_points ) );

  if (ec_options.refinementControl) {
    size_t num_moments = (data_rep->nonRandomIndices.empty()) ? 4 : 2;
    if (referenceMoments.empty())
      referenceMoments.sizeUninitialized(num_moments);
    if (deltaMoments.empty())
      deltaMoments.sizeUninitialized(num_moments);
  }
}


void HierarchInterpPolyApproximation::compute_expansion_coefficients()
{
  if (surrData.anchor()) {
    PCerr << "Error: anchor point not supported in HierarchInterpPoly"
	  << "Approximation::compute_expansion_coefficients" << std::endl;
    abort_handler(-1);
    /*
    if (expansionCoeffFlag) {
      expansionType1Coeffs[0][0][0] = surrData.anchor_function();
      if (data_rep->basisConfigOptions.useDerivs)
	Teuchos::setCol(surrData.anchor_gradient(), 0,
			expansionType2Coeffs[0][0]);
    }
    if (expansionCoeffGradFlag)
      Teuchos::setCol(surrData.anchor_gradient(), 0,
		      expansionType1CoeffGrads[0][0]);
    */
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key          = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, v, num_levels = key.size(), num_sets, num_tp_pts,
    cntr = 0, index, num_deriv_vars = surrData.num_derivative_variables();

  // level 0
  index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
  if (expansionCoeffFlag) {
    expansionType1Coeffs[0][0][0] = surrData.response_function(index);
    if (data_rep->basisConfigOptions.useDerivs)
      Teuchos::setCol(surrData.response_gradient(index), 0,
		      expansionType2Coeffs[0][0]);
  }
  if (expansionCoeffGradFlag)
    Teuchos::setCol(surrData.response_gradient(index), 0,
		    expansionType1CoeffGrads[0][0]);
  ++cntr;
  // levels 1 to num_levels
  for (lev=1; lev<num_levels; ++lev) {
    const UShort3DArray& key_l = key[lev];
    num_sets = key_l.size();
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = key_l[set].size();
      for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	const RealVector& c_vars = surrData.continuous_variables(index);
	// coefficients are hierarchical surpluses
	if (expansionCoeffFlag) {
	  expansionType1Coeffs[lev][set][pt] = surrData.response_function(index)
	    - value(c_vars, sm_mi, key, expansionType1Coeffs,
		    expansionType2Coeffs, lev-1);
	  if (data_rep->basisConfigOptions.useDerivs) {
	    const RealVector& data_grad = surrData.response_gradient(index);
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, key, expansionType1Coeffs, expansionType2Coeffs, lev-1);
	    Real* hier_grad = expansionType2Coeffs[lev][set][pt];
	    for (v=0; v<num_deriv_vars; ++v)
	      hier_grad[v] = data_grad[v] - prev_grad[v];
	  }
	}
	if (expansionCoeffGradFlag) {
	  const RealVector& data_grad = surrData.response_gradient(index);
	  const RealVector& prev_grad = gradient_nonbasis_variables(c_vars,
	    sm_mi, key, expansionType1CoeffGrads, lev-1);
	  Real* hier_grad = expansionType1CoeffGrads[lev][set][pt];
	  for (v=0; v<num_deriv_vars; ++v)
	    hier_grad[v] = data_grad[v] - prev_grad[v];
	}
      }
    }
  }

  computedMean = computedVariance
    = computedRefMean = computedDeltaMean
    = computedRefVariance = computedDeltaVariance = 0;
}


void HierarchInterpPolyApproximation::increment_coefficients()
{
  increment_current_from_reference();

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  switch (data_rep->expConfigOptions.refinementControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: // generalized sparse grids
    increment_coefficients(hsg_driver->trial_set());
    break;
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

  // size sobolIndices based on shared sobolIndexMap
  allocate_component_sobol();
}


// ***************************************************************
// TO DO: verify that decrement/push is always valid for surpluses
// ***************************************************************


void HierarchInterpPolyApproximation::decrement_coefficients()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShortArray& trial_set = data_rep->hsg_driver()->trial_set();
  size_t lev = l1_norm(trial_set);

  if (expansionCoeffFlag) {
    poppedExpT1Coeffs[trial_set] = expansionType1Coeffs[lev].back();
    expansionType1Coeffs[lev].pop_back();
    if (data_rep->basisConfigOptions.useDerivs) {
      poppedExpT2Coeffs[trial_set] = expansionType2Coeffs[lev].back();
      expansionType2Coeffs[lev].pop_back();
    }
  }
  if (expansionCoeffGradFlag) {
    poppedExpT1CoeffGrads[trial_set] = expansionType1CoeffGrads[lev].back();
    expansionType1CoeffGrads[lev].pop_back();
  }

  decrement_current_to_reference();
}


void HierarchInterpPolyApproximation::finalize_coefficients()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  const UShort3DArray& sm_mi = data_rep->hsg_driver()->smolyak_multi_index();

  size_t lev, set, num_sets, num_levels = sm_mi.size(), num_smolyak_sets,
    num_coeff_sets;
  for (lev=0; lev<num_levels; ++lev) {
    const UShort2DArray& sm_mi_l = sm_mi[lev];
    num_smolyak_sets = sm_mi_l.size();
    num_coeff_sets = (expansionCoeffFlag) ? expansionType1Coeffs[lev].size() :
      expansionType1CoeffGrads[lev].size();
    for (set=num_coeff_sets; set<num_smolyak_sets; ++set)
      push_coefficients(sm_mi_l[set]);
  }
  poppedExpT1Coeffs.clear(); poppedExpT2Coeffs.clear();
  poppedExpT1CoeffGrads.clear();

  computedMean = computedVariance = 0;
}


void HierarchInterpPolyApproximation::store_coefficients(size_t index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;

  size_t stored_len = storedExpType1Coeffs.size();
  if (index == _NPOS || index == stored_len) { // append
    if (expansionCoeffFlag) {
      storedExpType1Coeffs.push_back(expansionType1Coeffs);
      storedExpType2Coeffs.push_back(expansionType2Coeffs);
    }
    else { // keep indexing consistent
      storedExpType1Coeffs.push_back(RealVector2DArray());
      storedExpType2Coeffs.push_back(RealMatrix2DArray());
    }
    if (expansionCoeffGradFlag)
      storedExpType1CoeffGrads.push_back(expansionType1CoeffGrads);
    else // keep indexing consistent
      storedExpType1CoeffGrads.push_back(RealMatrix2DArray());
  }
  else if (index < stored_len) { // replace
    if (expansionCoeffFlag) {
      storedExpType1Coeffs[index] = expansionType1Coeffs;
      storedExpType2Coeffs[index] = expansionType2Coeffs;
    }
    else { // keep indexing consistent
      storedExpType1Coeffs[index] = RealVector2DArray();
      storedExpType2Coeffs[index] = RealMatrix2DArray();
    }
    if (expansionCoeffGradFlag)
      storedExpType1CoeffGrads[index] = expansionType1CoeffGrads;
    else // keep indexing consistent
      storedExpType1CoeffGrads[index] = RealMatrix2DArray();
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in HierarchInterpPoly"
	  << "Approximation::store_coefficients()" << std::endl;
    abort_handler(-1);
  }
}


void HierarchInterpPolyApproximation::restore_coefficients(size_t index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;

  size_t stored_len = storedExpType1Coeffs.size();
  if (index == _NPOS) {
    expansionType1Coeffs = storedExpType1Coeffs.back();
    expansionType2Coeffs = storedExpType2Coeffs.back();
    expansionType1CoeffGrads = storedExpType1CoeffGrads.back();
  }
  else if (index < stored_len) {
    expansionType1Coeffs = storedExpType1Coeffs[index];
    expansionType2Coeffs = storedExpType2Coeffs[index];
    expansionType1CoeffGrads = storedExpType1CoeffGrads[index];
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in HierarchInterpPoly"
	  << "Approximation::restore_coefficients()" << std::endl;
    abort_handler(-1);
  }
}


void HierarchInterpPolyApproximation::swap_coefficients(size_t index)
{
  if (expansionCoeffFlag) {
    std::swap(expansionType1Coeffs, storedExpType1Coeffs[index]);
    SharedHierarchInterpPolyApproxData* data_rep
      = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
    if (data_rep->basisConfigOptions.useDerivs)
      std::swap(expansionType2Coeffs, storedExpType2Coeffs[index]);
  }
  if (expansionCoeffGradFlag)
    std::swap(expansionType1CoeffGrads, storedExpType1CoeffGrads[index]);
}


void HierarchInterpPolyApproximation::remove_stored_coefficients(size_t index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;

  size_t stored_len = storedExpType1Coeffs.size();
  if (index == _NPOS || index == stored_len) {
    storedExpType1Coeffs.pop_back(); storedExpType2Coeffs.pop_back();
    storedExpType1CoeffGrads.pop_back();
  }
  else if (index < stored_len) {
    RealVector3DArray::iterator vit = storedExpType1Coeffs.begin();
    std::advance(vit, index); storedExpType1Coeffs.erase(vit);
    RealMatrix3DArray::iterator mit = storedExpType2Coeffs.begin();
    std::advance(mit, index); storedExpType2Coeffs.erase(mit);
    mit = storedExpType1CoeffGrads.begin();
    std::advance(mit, index); storedExpType1CoeffGrads.erase(mit);
  }
}


void HierarchInterpPolyApproximation::
combine_coefficients(short combine_type, size_t swap_index)
{
  if (swap_index != _NPOS) {
    swap_coefficients(swap_index);
    allocate_component_sobol(); // size sobolIndices from shared sobolIndexMap
  }

  // update expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads} by adding or
  // multiplying stored expansion evaluated at current collocation points
  size_t i, j, num_pts = surrData.points();
  Real curr_val, stored_val;
  /*
  for (i=0; i<num_pts; ++i) {
    const RealVector& c_vars = (anchor_pt && i == 0) ?
      surrData.anchor_continuous_variables() :
      surrData.continuous_variables(i);
    if (combine_type == MULT_COMBINE) { // eval once for both Coeffs/CoeffGrads
      stored_val = stored_value(c_vars);
      curr_val = expansionType1Coeffs[i]; // copy prior to update
    }
    if (expansionCoeffFlag) {
      // split up type1/type2 contribs so increments are performed properly
      if (combine_type == ADD_COMBINE)
	expansionType1Coeffs[i] += stored_value(c_vars);
      else if (combine_type == MULT_COMBINE)
	expansionType1Coeffs[i] *= stored_val;
      if (data_rep->basisConfigOptions.useDerivs) {
	const RealVector& stored_grad
	  = stored_gradient_basis_variables(c_vars);
	Real* exp_t2_coeffs_i = expansionType2Coeffs[i];
	size_t num_deriv_vars = stored_grad.length();
	if (combine_type == ADD_COMBINE)
	  for (j=0; j<num_deriv_vars; ++j)
	    exp_t2_coeffs_i[j] += stored_grad[j];
	else if (combine_type == MULT_COMBINE)
	  // hf = curr*stored --> dhf/dx = dcurr/dx*stored + curr*dstored/dx
	  for (j=0; j<num_deriv_vars; ++j)
	    exp_t2_coeffs_i[j] = exp_t2_coeffs_i[j] * stored_val
	                       + stored_grad[j]     * curr_val;
      }
    }
    if (expansionCoeffGradFlag) {
      Real* exp_t1_grad_i = expansionType1CoeffGrads[i];
      const RealVector& stored_grad
	= stored_gradient_nonbasis_variables(c_vars);
      size_t num_deriv_vars = stored_grad.length();
      if (combine_type == ADD_COMBINE)
	for (j=0; j<num_deriv_vars; ++j)
	  exp_t1_grad_i[j] += stored_grad[j];
      else if (combine_type == MULT_COMBINE)
	for (j=0; j<num_deriv_vars; ++j)
	  exp_t1_grad_i[j] = exp_t1_grad_i[j] * stored_val
	                   + stored_grad[j]   * curr_val;
    }
  }
  */

  // clear stored data now that it has been combined
  storedExpType1Coeffs.clear(); storedExpType2Coeffs.clear();
  storedExpType1CoeffGrads.clear();

  computedMean = computedVariance = 0;
}


/** Lower level helper function to process a single index set. */
void HierarchInterpPolyApproximation::
increment_coefficients(const UShortArray& index_set)
{
  size_t lev, old_levels = expansionType1Coeffs.size(), set, old_sets,
    pt, old_pts = 0;
  for (lev=0; lev<old_levels; ++lev) {
    old_sets = expansionType1Coeffs[lev].size();
    for (set=0; set<old_sets; ++set)
      old_pts += (expansionCoeffFlag) ?
	expansionType1Coeffs[lev][set].length()	:
	expansionType1CoeffGrads[lev][set].numCols();
  }
  lev = l1_norm(index_set);
  if (lev >= old_levels) {
    expansionType1Coeffs.resize(lev+1);
    expansionType2Coeffs.resize(lev+1);
    expansionType1CoeffGrads.resize(lev+1);
  }
  set = expansionType1Coeffs[lev].size();
  // append empty and update in place
  RealVector fns; RealMatrix grads;
  expansionType1Coeffs[lev].push_back(fns);
  expansionType2Coeffs[lev].push_back(grads);
  expansionType1CoeffGrads[lev].push_back(grads);
  RealVector& t1_coeffs      = expansionType1Coeffs[lev][set];
  RealMatrix& t2_coeffs      = expansionType2Coeffs[lev][set];
  RealMatrix& t1_coeff_grads = expansionType1CoeffGrads[lev][set];

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray& key   = hsg_driver->collocation_key();
  size_t index, num_trial_pts = key[lev][set].size(), v, num_deriv_vars = 0;
  if (expansionCoeffFlag) {
    t1_coeffs.sizeUninitialized(num_trial_pts);
    if (data_rep->basisConfigOptions.useDerivs) {
      num_deriv_vars = expansionType2Coeffs[0][0].numRows();
      t2_coeffs.shapeUninitialized(num_deriv_vars, num_trial_pts);
    }
  }
  if (expansionCoeffGradFlag) {
    num_deriv_vars = expansionType1CoeffGrads[0][0].numRows();
    t1_coeff_grads.shapeUninitialized(num_deriv_vars, num_trial_pts);
  }
 
  for (pt=0, index=old_pts; pt<num_trial_pts; ++pt, ++index) {
    const RealVector& c_vars = surrData.continuous_variables(index);
    if (expansionCoeffFlag) {
      t1_coeffs[pt] = surrData.response_function(index) - value(c_vars, sm_mi,
	key, expansionType1Coeffs, expansionType2Coeffs, lev-1);
      if (data_rep->basisConfigOptions.useDerivs) {
	const RealVector& data_grad = surrData.response_gradient(index);
	const RealVector& prev_grad = gradient_basis_variables(c_vars, sm_mi,
	  key, expansionType1Coeffs, expansionType2Coeffs, lev-1);
	Real* hier_grad = t2_coeffs[pt];
	for (v=0; v<num_deriv_vars; ++v)
	  hier_grad[v] = data_grad[v] - prev_grad[v];
      }
    }
    if (expansionCoeffGradFlag) {
      const RealVector& data_grad = surrData.response_gradient(index);
      const RealVector& prev_grad = gradient_nonbasis_variables(c_vars, sm_mi,
	key, expansionType1CoeffGrads, lev-1);
      Real* hier_grad = t1_coeff_grads[pt];
      for (v=0; v<num_deriv_vars; ++v)
	hier_grad[v] = data_grad[v] - prev_grad[v];
    }
  }
}


/** Lower level helper function to process a single index set. */
void HierarchInterpPolyApproximation::
push_coefficients(const UShortArray& push_set)
{
  size_t lev = l1_norm(push_set);
  if (expansionCoeffFlag) {
    expansionType1Coeffs[lev].push_back(poppedExpT1Coeffs[push_set]);
    poppedExpT1Coeffs.erase(push_set);
    SharedHierarchInterpPolyApproxData* data_rep
      = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
    if (data_rep->basisConfigOptions.useDerivs) {
      expansionType2Coeffs[lev].push_back(poppedExpT2Coeffs[push_set]);
      poppedExpT2Coeffs.erase(push_set);
    }
  }
  if (expansionCoeffGradFlag) {
    expansionType1CoeffGrads[lev].push_back(poppedExpT1CoeffGrads[push_set]);
    poppedExpT1CoeffGrads.erase(push_set);
  }
}


void HierarchInterpPolyApproximation::increment_current_from_reference()
{
  computedRefMean     = computedMean;
  computedRefVariance = computedVariance;

  if ( (computedMean & 1) || (computedVariance & 1) )
    referenceMoments = numericalMoments;
  if (computedMean & 2)
    meanRefGradient = meanGradient;
  if (computedVariance & 2)
    varianceRefGradient = varianceGradient;

  // clear current and delta
  computedMean = computedVariance =
    computedDeltaMean = computedDeltaVariance = 0;
}


void HierarchInterpPolyApproximation::decrement_current_to_reference()
{
  computedMean     = computedRefMean;
  computedVariance = computedRefVariance;

  if ( (computedRefMean & 1) || (computedRefVariance & 1) )
    numericalMoments = referenceMoments;
  if (computedRefMean & 2)
    meanGradient = meanRefGradient;
  if (computedRefVariance & 2)
    varianceGradient = varianceRefGradient;

  // leave reference settings, but clear delta settings
  computedDeltaMean = computedDeltaVariance = 0;
}


Real HierarchInterpPolyApproximation::
value(const RealVector& x, const UShort3DArray& sm_mi, const UShort4DArray& key,
      const RealVector2DArray& t1_coeffs, const RealMatrix2DArray& t2_coeffs,
      unsigned short level)
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
  size_t lev, set, num_sets;
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    // sm_mi and key include all current index sets, whereas the t1/t2
    // coeffs may refect a partial state derived from a reference_key
    num_sets = t1_coeffs_l.size();
    for (set=0; set<num_sets; ++set)
      approx_val +=
	data_rep->tensor_product_value(x, t1_coeffs_l[set], t2_coeffs_l[set],
				       sm_mi_l[set], key_l[set], colloc_index);
  }
  return approx_val;
}


/** All variables version. */
Real HierarchInterpPolyApproximation::
value(const RealVector& x, const UShort3DArray& sm_mi, const UShort4DArray& key,
      const RealVector2DArray& t1_coeffs, const RealMatrix2DArray& t2_coeffs,
      unsigned short level, const SizetList& subset_indices)
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
  size_t lev, set, num_sets;
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    // sm_mi and key include all current index sets, whereas the t1/t2
    // coeffs may refect a partial state derived from a reference_key
    num_sets = t1_coeffs_l.size();
    for (set=0; set<num_sets; ++set)
      approx_val +=
	data_rep->tensor_product_value(x, t1_coeffs_l[set], t2_coeffs_l[set],
				       sm_mi_l[set], key_l[set], colloc_index,
				       subset_indices);
  }
  return approx_val;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 unsigned short level)
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
  size_t lev, set, num_sets;
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    // sm_mi and key include all current index sets, whereas the t1/t2
    // coeffs may refect a partial state derived from a reference_key
    num_sets = t1_coeffs_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 unsigned short level, const SizetList& subset_indices)
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
  size_t lev, set, num_sets;
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    // sm_mi and key include all current index sets, whereas the t1/t2
    // coeffs may refect a partial state derived from a reference_key
    num_sets = t1_coeffs_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index,
	  subset_indices);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			 const UShort4DArray& key,
			 const RealVector2DArray& t1_coeffs,
			 const RealMatrix2DArray& t2_coeffs,
			 const SizetArray& dvv, unsigned short level)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in HierarchInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t lev, set, num_sets, num_deriv_vars = dvv.size();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.;

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&       sm_mi_l = sm_mi[lev];
    const UShort3DArray&         key_l = key[lev];
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    const RealMatrixArray& t2_coeffs_l = t2_coeffs[lev];
    // sm_mi and key include all current index sets, whereas the t1/t2
    // coeffs may refect a partial state derived from a reference_key
    num_sets = t1_coeffs_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_basis_variables(x, t1_coeffs_l[set],
	  t2_coeffs_l[set], sm_mi_l[set], key_l[set], colloc_index, dvv);
  }

  return approxGradient;
}


const RealVector& HierarchInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x, const UShort3DArray& sm_mi,
			    const UShort4DArray& key,
			    const RealMatrix2DArray& t1_coeff_grads,
			    unsigned short level)
{
  // Error check for required data
  size_t lev, set, num_sets, num_deriv_vars;
  if (expansionCoeffGradFlag) {
    if (t1_coeff_grads.size() > level && t1_coeff_grads[level].size())
      num_deriv_vars = t1_coeff_grads[level][0].numRows();
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

  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.;

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> 2DArrays allow default indexing
  for (lev=0; lev<=level; ++lev) {
    const UShort2DArray&            sm_mi_l = sm_mi[lev];
    const UShort3DArray&              key_l = key[lev];
    const RealMatrixArray& t1_coeff_grads_l = t1_coeff_grads[lev];
    // sm_mi and key include all current index sets, whereas the t1/t2
    // coeffs may refect a partial state derived from a reference_key
    num_sets = t1_coeff_grads_l.size();
    for (set=0; set<num_sets; ++set)
      approxGradient +=
	data_rep->tensor_product_gradient_nonbasis_variables(x,
	  t1_coeff_grads_l[set], sm_mi_l[set], key_l[set], colloc_index);
  }

  return approxGradient;
}


const RealSymMatrix& HierarchInterpPolyApproximation::
hessian_basis_variables(const RealVector& x)
{
  PCerr << "Error: HierarchInterpPolyApproximation::hessian_basis_variables() "
	<< "not yet implemented." << std::endl;
  abort_handler(-1);

  return approxHessian;
}


Real HierarchInterpPolyApproximation::mean()
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedMean & 1))
    return numericalMoments[0];

  Real mean = expectation(expansionType1Coeffs, expansionType2Coeffs);
  if (std_mode)
    { numericalMoments[0] = mean; computedMean |= 1; }
  return mean;
}



Real HierarchInterpPolyApproximation::mean(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedMean & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean))
    return numericalMoments[0];

  Real mean = expectation(x, expansionType1Coeffs, expansionType2Coeffs);
  if (all_mode)
    { numericalMoments[0] = mean; computedMean |= 1; xPrevMean = x; }
  return mean;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& HierarchInterpPolyApproximation::mean_gradient()
{
  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Hierarch"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedMean & 2))
    return meanGradient;

  meanGradient = expectation_gradient(expansionType1CoeffGrads);
  if (std_mode) computedMean |=  2; //   activate 2-bit
  else          computedMean &= ~2; // deactivate 2-bit: protect mixed usage
  return meanGradient;
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
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if ( all_mode && (computedMean & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevMeanGrad) ) // && dvv == dvvPrev)
    return meanGradient;

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
  if (meanGradient.length() != num_deriv_vars)
    meanGradient.sizeUninitialized(num_deriv_vars);
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
      meanGradient[i] = expectation_gradient(x, expansionType1CoeffGrads, cntr);
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
      meanGradient[i] = expectation_gradient(x, expansionType1Coeffs,
					     expansionType2Coeffs, deriv_index);
    }
  }
  if (all_mode) { computedMean |=  2; xPrevMeanGrad = x; }
  else            computedMean &= ~2; // deactivate 2-bit: protect mixed usage
  return meanGradient;
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

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if (same && std_mode && (computedVariance & 1))
    return numericalMoments[1];

  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  Real mean_1 = mean(), mean_2 = (same) ? mean_1 : hip_approx_2->mean();
  central_product_interpolant(hip_approx_2, mean_1, mean_2,
			      cov_t1_coeffs, cov_t2_coeffs);

  // evaluate expectation of these t1/t2 coefficients
  Real covar = expectation(cov_t1_coeffs, cov_t2_coeffs);
  // Note: separation of reference and increment using cov_t{1,2}_coeffs
  // with {ref,incr}_key would provide an increment of a central moment
  // around an invariant center.  For hierarchical covariance, one must
  // also account for the change in mean as in delta_covariance().

  if (same && std_mode)
    { numericalMoments[1] = covar; computedVariance |= 1; }
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

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if ( same && all_mode && (computedVariance & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar) )
    return numericalMoments[1];

  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  Real mean_1 = mean(x), mean_2 = (same) ? mean_1 : hip_approx_2->mean(x);
  central_product_interpolant(hip_approx_2, mean_1, mean_2,
			      cov_t1_coeffs, cov_t2_coeffs);

  // evaluate expectation of these t1/t2 coefficients
  Real covar = expectation(x, cov_t1_coeffs, cov_t2_coeffs);

  if (same && all_mode)
    { numericalMoments[1] = covar; computedVariance |= 1; xPrevVar = x; }
  return covar;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design variables are inserted and some
    are augmented) requires no special treatment. */
const RealVector& HierarchInterpPolyApproximation::variance_gradient()
{
  // Error check for required data
  if (!expansionCoeffFlag ||
      !expansionCoeffGradFlag) {
    PCerr << "Error: insufficient expansion coefficient data in HierarchInterp"
	  << "PolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedVariance & 2))
    return varianceGradient;

  Real mean_1 = mean(); const RealVector& mean1_grad = mean_gradient();
  RealMatrix2DArray cov_t1_coeff_grads;
  central_product_gradient_interpolant(this, mean_1, mean_1, mean1_grad,
				       mean1_grad, cov_t1_coeff_grads);
  varianceGradient = expectation_gradient(cov_t1_coeff_grads);
  if (std_mode) computedVariance |=  2;
  else          computedVariance &= ~2; // deactivate 2-bit: protect mixed usage
  return varianceGradient;
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
  if ( all_mode && (computedVariance & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevVarGrad) ) // && dvv == dvvPrev)
    return varianceGradient;

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

  RealVector2DArray cov_t1_coeffs;
  RealMatrix2DArray cov_t1_coeff_grads, cov_t2_coeffs;
  Real mean_1 = mean(x);
  if (insert) {
    const RealVector& mean1_grad = mean_gradient(x, dvv);
    central_product_gradient_interpolant(this, mean_1, mean_1, mean1_grad,
					 mean1_grad, cov_t1_coeff_grads);
  }
  if (augment)
    central_product_interpolant(this, mean_1, mean_1, cov_t1_coeffs,
				cov_t2_coeffs);

  if (varianceGradient.length() != num_deriv_vars)
    varianceGradient.sizeUninitialized(num_deriv_vars);
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    Real& grad_i = varianceGradient[i];
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
      varianceGradient[i] = expectation_gradient(x, cov_t1_coeff_grads, cntr);
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
      varianceGradient[i] = expectation_gradient(x, cov_t1_coeffs,
						 cov_t2_coeffs, deriv_index);
    }
  }
  if (all_mode) { computedVariance |=  2; xPrevVarGrad = x; }
  else            computedVariance &= ~2;//deactivate 2-bit: protect mixed usage
  return varianceGradient;
}


Real HierarchInterpPolyApproximation::
reference_mean(const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedRefMean & 1))
    return referenceMoments[0];

  Real ref_mean
    = expectation(expansionType1Coeffs, expansionType2Coeffs, ref_key);
  if (std_mode)
    { referenceMoments[0] = ref_mean; computedRefMean |= 1; }
  return ref_mean;
}


Real HierarchInterpPolyApproximation::
reference_mean(const RealVector& x, const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedRefMean & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevRefMean))
    return referenceMoments[0];

  Real ref_mean
    = expectation(x, expansionType1Coeffs, expansionType2Coeffs, ref_key);
  if (all_mode)
    { referenceMoments[0] = ref_mean; computedRefMean |= 1; xPrevRefMean = x; }
  return ref_mean;
}


Real HierarchInterpPolyApproximation::
reference_variance(const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedRefVariance & 1))
    return referenceMoments[1];

  Real ref_mean = reference_mean(ref_key);
  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  central_product_interpolant(this, ref_mean, ref_mean, cov_t1_coeffs,
			      cov_t2_coeffs, ref_key);
  Real ref_var = expectation(cov_t1_coeffs, cov_t2_coeffs, ref_key);
  if (std_mode)
    { referenceMoments[1] = ref_var; computedRefVariance |= 1; }
  return ref_var;
}


Real HierarchInterpPolyApproximation::
reference_variance(const RealVector& x, const UShort2DArray& ref_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedRefVariance & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevRefVar))
    return referenceMoments[1];

  Real ref_mean = reference_mean(x, ref_key);
  RealVector2DArray cov_t1_coeffs; RealMatrix2DArray cov_t2_coeffs;
  central_product_interpolant(this, ref_mean, ref_mean, cov_t1_coeffs,
			      cov_t2_coeffs, ref_key);
  Real ref_var = expectation(x, cov_t1_coeffs, cov_t2_coeffs, ref_key);
  if (all_mode) {
    referenceMoments[1] = ref_var;
    computedRefVariance |= 1; xPrevRefVar = x;
  }
  return ref_var;
}


Real HierarchInterpPolyApproximation::delta_mean(const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedDeltaMean & 1))
    return deltaMoments[0];

  Real delta_mean
    = expectation(expansionType1Coeffs, expansionType2Coeffs, incr_key);
  if (std_mode)
    { deltaMoments[0] = delta_mean; computedDeltaMean |= 1; }
  return delta_mean;
}


Real HierarchInterpPolyApproximation::
delta_mean(const RealVector& x, const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedDeltaMean & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaMean))
    return deltaMoments[0];

  Real delta_mean
    = expectation(x, expansionType1Coeffs, expansionType2Coeffs, incr_key);
  if (all_mode) {
    deltaMoments[0] = delta_mean;
    computedDeltaMean |= 1; xPrevDeltaMean = x;
  }
  return delta_mean;
}


Real HierarchInterpPolyApproximation::
delta_variance(const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedDeltaVariance & 1))
    return deltaMoments[1];

  Real delta_var = delta_covariance(this, ref_key, incr_key);
  if (std_mode)
    { deltaMoments[1] = delta_var; computedDeltaVariance |= 1; }
  return delta_var;
}


Real HierarchInterpPolyApproximation::
delta_variance(const RealVector& x, const UShort2DArray& ref_key,
	       const UShort2DArray& incr_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedDeltaVariance & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevDeltaVar))
    return deltaMoments[1];

  Real delta_var = delta_covariance(x, this, ref_key, incr_key);
  if (all_mode) {
    deltaMoments[1] = delta_var;
    computedDeltaVariance |= 1; xPrevDeltaVar = x;
  }
  return delta_var;
}


Real HierarchInterpPolyApproximation::
delta_std_deviation(const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  // delta-sigma = sqrt( var0 + delta-var ) - sigma0
  //             = [ sqrt(1 + delta_var / var0) - 1 ] * sigma0
  //             = sqrt1pm1(delta_var / var0) * sigma0
  // where sqrt1pm1(x) = expm1[ log1p(x) / 2 ]

  Real delta_var = delta_variance(ref_key, incr_key),
       var0      = reference_variance(ref_key),
       sigma0    = std::sqrt(var0);

  return (delta_var < var0) ?
    bmth::sqrt1pm1(delta_var / var0) * sigma0 :            // preserve precision
    std::sqrt(var0 + delta_var) - sigma0; // precision OK; prevent div by var0=0
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
       sigma0    = std::sqrt(var0);

  return (delta_var < var0) ?
    bmth::sqrt1pm1(delta_var / var0) * sigma0 :            // preserve precision
    std::sqrt(var0 + delta_var) - sigma0; // precision OK; prevent div by var0=0
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
  Real beta0, sigma0 = std::sqrt(var0), sigma1 = sigma0 + delta_sigma;
  if (cdf_flag) {
    if (sigma0 > SMALL_NUMBER && sigma1 > SMALL_NUMBER) {
      beta0 = (mu0 - z_bar) / sigma0;
      return ( delta_mu - delta_sigma * beta0) / sigma1;
    }
    else if (sigma1 > SMALL_NUMBER)// neglect beta0 term (zero init reliability)
      return delta_mu / sigma1; // or delta = beta1 = (mu1 - z_bar) / sigma1 ?
    else if (sigma0 > SMALL_NUMBER) // assume beta1 = 0 -> delta = -beta0
      return (z_bar - mu0) / sigma0;
    else                      // assume beta0 = beta1 = 0
      return 0;
  }
  else {
    if (sigma0 > SMALL_NUMBER && sigma1 > SMALL_NUMBER) {
      beta0 = (z_bar - mu0) / sigma0;
      return (-delta_mu - delta_sigma * beta0) / sigma1;
    }
    else if (sigma1 > SMALL_NUMBER)// neglect beta0 term (zero init reliability)
      return -delta_mu / sigma1;
    else if (sigma0 > SMALL_NUMBER) // assume beta1 = 0 -> delta = -beta0
      return (mu0 - z_bar) / sigma0;
    else                      // assume beta0 = beta1 = 0
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
delta_covariance(PolynomialApproximation* poly_approx_2,
		 const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  // Supports multiple grid increments in discerning nominal from delta based
  // on isotropic/anisotropic/generalized index set increments.  In current
  // use, 2D keys with set ranges are sufficient: level -> {start,end} set.
  // In the future, may need 3D keys for level/set/point.
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool same = (this == hip_approx_2);

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::delta_covariance()" << std::endl;
    abort_handler(-1);
  }

  RealVector2DArray r1r2_t1_coeffs; RealMatrix2DArray r1r2_t2_coeffs;
  product_interpolant(hip_approx_2, r1r2_t1_coeffs, r1r2_t2_coeffs);

  // Compute surplus for r1, r2, and r1r2 and retrieve reference mean values
  Real ref_mean_r1  = reference_mean(ref_key),
    delta_mean_r1   = delta_mean(incr_key),
    ref_mean_r2     = (same) ? ref_mean_r1 :
      expectation(hip_approx_2->expansionType1Coeffs,
		  hip_approx_2->expansionType2Coeffs, ref_key),
    delta_mean_r2   = (same) ? delta_mean_r1 :
      expectation(hip_approx_2->expansionType1Coeffs,
		  hip_approx_2->expansionType2Coeffs, incr_key),
    delta_mean_r1r2 = expectation(r1r2_t1_coeffs, r1r2_t2_coeffs, incr_key);

  // Hierarchical increment to covariance:
  // \Delta\Sigma_ij = \Sigma^1_ij - \Sigma^0_ij
  //   = ( E[Ri Rj]^1 - E[Ri]^1 E[Rj]^1 ) - ( E[Ri Rj]^0 - E[Ri]^0 E[Rj]^0 )
  //   = E[Ri Rj]^0 + \DeltaE[Ri Rj]
  //     - (E[Ri]^0 + \DeltaE[Ri]) (E[Rj]^0 + \DeltaE[Rj])
  //     - E[Ri Rj]^0 + E[Ri]^0 E[Rj]^0
  //   = \DeltaE[Ri Rj] - \DeltaE[Ri] E[Rj]^0 - E[Ri]^0 \DeltaE[Rj]
  //     - \DeltaE[Ri] \DeltaE[Rj]
  Real delta_covar = delta_mean_r1r2 - ref_mean_r1 * delta_mean_r2
     - ref_mean_r2 * delta_mean_r1 - delta_mean_r1 * delta_mean_r2;
  if (same && data_rep->nonRandomIndices.empty()) // std mode
    { deltaMoments[1] = delta_covar; computedDeltaVariance |= 1; }
  return delta_covar;
}


Real HierarchInterpPolyApproximation::
delta_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2,
		 const UShort2DArray& ref_key, const UShort2DArray& incr_key)
{
  HierarchInterpPolyApproximation* hip_approx_2 = 
    (HierarchInterpPolyApproximation*)poly_approx_2;
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  bool same = (this == hip_approx_2);

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !hip_approx_2->expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "HierarchInterpPolyApproximation::delta_covariance()" << std::endl;
    abort_handler(-1);
  }

  RealVector2DArray r1r2_t1_coeffs; RealMatrix2DArray r1r2_t2_coeffs;
  product_interpolant(hip_approx_2, r1r2_t1_coeffs, r1r2_t2_coeffs);

  // Compute surplus for r1, r2, and r1r2 and retrieve reference mean values
  Real ref_mean_r1  = reference_mean(x, ref_key),
    delta_mean_r1   = delta_mean(x, incr_key),
    ref_mean_r2     = (same) ? ref_mean_r1 :
      expectation(x, hip_approx_2->expansionType1Coeffs,
		  hip_approx_2->expansionType2Coeffs, ref_key),
    delta_mean_r2   = (same) ? delta_mean_r1 :
      expectation(x, hip_approx_2->expansionType1Coeffs,
		  hip_approx_2->expansionType2Coeffs, incr_key),
    delta_mean_r1r2 = expectation(x, r1r2_t1_coeffs, r1r2_t2_coeffs, incr_key);

  // same expression as standard expansion mode case above
  Real delta_covar = delta_mean_r1r2 - ref_mean_r1 * delta_mean_r2
     - ref_mean_r2 * delta_mean_r1 - delta_mean_r1 * delta_mean_r2;
  if (same && !data_rep->nonRandomIndices.empty()) { // all vars mode
    deltaMoments[1] = delta_covar;
    computedDeltaVariance |= 1; xPrevDeltaVar = x;
  }
  return delta_covar;
}


Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs, const RealVector2DArray& t1_wts,
	    const RealMatrix2DArray& t2_coeffs, const RealMatrix2DArray& t2_wts,
	    const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  Real integral = 0.;
  size_t lev, set, pt, num_levels = t1_coeffs.size(), set_start = 0, set_end,
    num_tp_pts;
  bool partial = !set_partition.empty();
  if (data_rep->basisConfigOptions.useDerivs) {
    size_t v, num_v = sharedDataRep->numVars;
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1_coeffs_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealMatrix& t2_coeffs_ls = t2_coeffs[lev][set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	const RealMatrix&    t2_wts_ls = t2_wts[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
	  const Real* t2_coeffs_lsp = t2_coeffs_ls[pt];
	  const Real* t2_wts_lsp    = t2_wts_ls[pt];
	  for (v=0; v<num_v; ++v)
	    integral += t2_coeffs_lsp[v] * t2_wts_lsp[v];
	}
      }
    }
  }
  else {
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1_coeffs_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) // omitted if empty surplus vector
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
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
  const RealVector2DArray& t1_wts = hsg_driver->type1_weight_set_arrays();
  size_t lev, set, pt, num_levels = t1_coeffs.size(), num_sets,
    tp_pt_start = 0, tp_pt_end;
  bool partial = !pt_partition.empty();
  switch (data_rep->basisConfigOptions.useDerivs) {
  case false:
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      num_sets = t1_coeffs_l.size();
      for (set=0; set<num_sets; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	if (partial) {
	  tp_pt_start = pt_partition[lev][set][0];
	  tp_pt_end   = pt_partition[lev][set][1];
	}
	else
	  tp_pt_end   = t1_coeffs_ls.length();
	for (pt=tp_pt_start; pt<tp_pt_end; ++pt)
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
      }
    }
    break;
  case true: {
    const RealMatrix2DArray& t2_wts = hsg_driver->type2_weight_set_arrays();
    size_t v, num_v = sharedDataRep->numVars;
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      num_sets = t1_coeffs_l.size();
      for (set=0; set<num_sets; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealMatrix& t2_coeffs_ls = t2_coeffs[lev][set];
	const RealVector&    t1_wts_ls = t1_wts[lev][set];
	const RealMatrix&    t2_wts_ls = t2_wts[lev][set];
	if (partial) {
	  tp_pt_start = pt_partition[lev][set][0];
	  tp_pt_end   = pt_partition[lev][set][1];
	}
	else
	  tp_pt_end   = t1_coeffs_ls.length();
	for (pt=tp_pt_start; pt<tp_pt_end; ++pt) {
	  integral += t1_coeffs_ls[pt] * t1_wts_ls[pt];
	  const Real* t2_coeffs_lsp = t2_coeffs_ls[pt];
	  const Real* t2_wts_lsp    = t2_wts_ls[pt];
	  for (v=0; v<num_v; ++v)
	    integral += t2_coeffs_lsp[v] * t2_wts_lsp[v];
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
expectation(const RealVector& x, const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs,
	    const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi      = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key = hsg_driver->collocation_key();

  Real integral = 0.;
  size_t lev, set, pt, num_levels = t1_coeffs.size(), set_start = 0, set_end,
    num_tp_pts;
  bool partial = !set_partition.empty();
  if (data_rep->basisConfigOptions.useDerivs) {
    size_t v, num_v = sharedDataRep->numVars;
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1_coeffs_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const RealMatrix& t2_coeffs_ls = t2_coeffs[lev][set];
	const UShortArray&    sm_mi_ls = sm_mi[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	  const UShortArray& key_lsp = colloc_key[lev][set][pt];
	  integral += t1_coeffs_ls[pt]
	    * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls,
						data_rep->nonRandomIndices)
	    * data_rep->type1_weight(key_lsp, sm_mi_ls,
				     data_rep->randomIndices);
	  const Real* t2_coeffs_lsp = t2_coeffs_ls[pt];
	  for (v=0; v<num_v; ++v)
	    integral += t2_coeffs_lsp[v]
	      * data_rep->type2_interpolant_value(x, v, key_lsp, sm_mi_ls,
						  data_rep->nonRandomIndices)
	      * data_rep->type2_weight(v, key_lsp, sm_mi_ls,
				       data_rep->randomIndices);
	}
      }
    }
  }
  else {
    for (lev=0; lev<num_levels; ++lev) {
      const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
      if (partial)
	{ set_start = set_partition[lev][0]; set_end = set_partition[lev][1]; }
      else
	set_end = t1_coeffs_l.size();
      for (set=set_start; set<set_end; ++set) {
	const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
	const UShortArray&    sm_mi_ls = sm_mi[lev][set];
	num_tp_pts = t1_coeffs_ls.length();
	for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	  const UShortArray& key_lsp = colloc_key[lev][set][pt];
	  integral += t1_coeffs_ls[pt]
	    * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls,
						data_rep->nonRandomIndices)
	    * data_rep->type1_weight(key_lsp, sm_mi_ls,
				     data_rep->randomIndices);
	}
      }
    }
  }
  return integral;
}


/** For inserted/augmented design/epistemic variables in standard mode. */
const RealVector& HierarchInterpPolyApproximation::
expectation_gradient(const RealMatrix2DArray& t1_coeff_grads,
		     const RealVector2DArray& t1_wts)
{
  size_t lev, num_levels = t1_coeff_grads.size(), set, num_sets, pt, num_tp_pts,
    v, num_deriv_vars = t1_coeff_grads[0][0].numRows();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.;

  for (lev=0; lev<num_levels; ++lev) {
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
		     const RealMatrix2DArray& t1_coeff_grads, size_t t1cg_index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi      = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key = hsg_driver->collocation_key();

  size_t lev, num_levels = t1_coeff_grads.size(), set, num_sets, pt, num_tp_pts;

  Real grad = 0.;
  for (lev=0; lev<num_levels; ++lev) {
    const RealMatrixArray& t1_coeff_grads_l = t1_coeff_grads[lev];
    num_sets = t1_coeff_grads_l.size();
    for (set=0; set<num_sets; ++set) {
      const RealMatrix& t1_coeff_grads_ls = t1_coeff_grads_l[set];
      const UShortArray&         sm_mi_ls = sm_mi[lev][set];
      num_tp_pts = t1_coeff_grads_ls.numCols();
      for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	const UShortArray& key_lsp = colloc_key[lev][set][pt];
	grad += t1_coeff_grads_ls(t1cg_index, pt)
	  * data_rep->type1_interpolant_value(x, key_lsp, sm_mi_ls,
					      data_rep->nonRandomIndices)
	  * data_rep->type1_weight(key_lsp, sm_mi_ls, data_rep->randomIndices);
      }
    }
  }
  return grad;
}


/** For augmented design/epistemic variables in all_variables mode. */
Real HierarchInterpPolyApproximation::
expectation_gradient(const RealVector& x, const RealVector2DArray& t1_coeffs,
		     const RealMatrix2DArray& t2_coeffs, size_t deriv_index)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi      = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      colloc_key = hsg_driver->collocation_key();

  size_t lev, num_levels = t1_coeffs.size(), set, num_sets, pt, num_tp_pts, v,
    num_v = sharedDataRep->numVars;

  Real grad = 0.;
  for (lev=0; lev<num_levels; ++lev) {
    const RealVectorArray& t1_coeffs_l = t1_coeffs[lev];
    num_sets = t1_coeffs_l.size();
    for (set=0; set<num_sets; ++set) {
      const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
      const UShortArray&    sm_mi_ls = sm_mi[lev][set];
      num_tp_pts = t1_coeffs_ls.length();
      for (pt=0; pt<num_tp_pts; ++pt) { // omitted if empty surplus vector
	const UShortArray& key_lsp = colloc_key[lev][set][pt];
	grad += t1_coeffs_ls[pt]
	  * data_rep->type1_interpolant_gradient(x, deriv_index, key_lsp,
						 sm_mi_ls,
						 data_rep->nonRandomIndices)
	  * data_rep->type1_weight(key_lsp, sm_mi_ls, data_rep->randomIndices);
	if (data_rep->basisConfigOptions.useDerivs) {
	  const Real *t2_coeff_lsp = t2_coeffs[lev][set][pt];
	  for (v=0; v<num_v; ++v)
	    grad += t2_coeff_lsp[v]
	      * data_rep->type2_interpolant_gradient(x, deriv_index, v, key_lsp,
						     sm_mi_ls,
						     data_rep->nonRandomIndices)
	      * data_rep->type2_weight(v, key_lsp, sm_mi_ls,
				       data_rep->randomIndices);
	}
      }
    }
  }

  return grad;
}


/** Whereas expectation() supports either a reference or increment key
    (passed as generic set_partition), functions forming hierarchical
    interpolant coefficients support only a reference key (starting
    point must be set 0; end point can be controlled). */
void HierarchInterpPolyApproximation::
product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
		    RealVector2DArray& r1r2_t1_coeffs,
		    RealMatrix2DArray& r1r2_t2_coeffs,
		    const UShort2DArray& reference_key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key          = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, num_levels = expansionType1Coeffs.size(),
    num_sets, num_tp_pts, cntr = 0, index;
  bool partial = !reference_key.empty();
  const SurrogateData& s_data_2 = hip_approx_2->surrData;

  // form hierarchical t1/t2 coeffs for raw moment R1 R2
  r1r2_t1_coeffs.resize(num_levels); r1r2_t1_coeffs[0].resize(1);
  r1r2_t2_coeffs.resize(num_levels); r1r2_t2_coeffs[0].resize(1);
  r1r2_t1_coeffs[0][0].sizeUninitialized(1);
  if (data_rep->basisConfigOptions.useDerivs) {
    size_t v, num_v = sharedDataRep->numVars;
    // level 0 (assume this is always contained in a partial reference_key)
    index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
    Real data_fn1 = surrData.response_function(index);
    Real data_fn2 = s_data_2.response_function(index);
    r1r2_t1_coeffs[0][0][0] = data_fn1 * data_fn2;
    r1r2_t2_coeffs[0][0].shapeUninitialized(num_v, 1);
    Real *r1r2_t2_coeffs_000 = r1r2_t2_coeffs[0][0][0];
    const RealVector& data_grad1 = surrData.response_gradient(index);
    const RealVector& data_grad2 = s_data_2.response_gradient(index);
    for (v=0; v<num_v; ++v)
      r1r2_t2_coeffs_000[v]
	= data_fn1 * data_grad2[v] + data_fn2 * data_grad1[v];
    ++cntr;
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? reference_key[lev][1] : key[lev].size();
      r1r2_t1_coeffs[lev].resize(num_sets);
      r1r2_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	RealMatrix& r1r2_t2_coeffs_ls = r1r2_t2_coeffs[lev][set];
	num_tp_pts = key[lev][set].size();
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	r1r2_t2_coeffs_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  const RealVector& c_vars = surrData.continuous_variables(index);
	  // type1 hierarchical interpolation of R1 R2
	  data_fn1 = surrData.response_function(index);
	  data_fn2 = s_data_2.response_function(index);
	  r1r2_t1_coeffs_ls[pt] = data_fn1 * data_fn2 -
	    value(c_vars, sm_mi, key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1);
	  // type2 hierarchical interpolation of R1 R2
	  // --> interpolated grads are R1 * R2' + R2 * R1'
	  Real* r1r2_t2_coeffs_lsp = r1r2_t2_coeffs_ls[pt];
	  const RealVector& data_grad1 = surrData.response_gradient(index);
	  const RealVector& data_grad2 = s_data_2.response_gradient(index);
	  const RealVector& prev_grad  = gradient_basis_variables(c_vars,
	    sm_mi, key, r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1);
	  for (v=0; v<num_v; ++v)
	    r1r2_t2_coeffs_lsp[v] = data_fn1 * data_grad2[v]
	      + data_fn2 * data_grad1[v] - prev_grad[v];
	}
      }
    }
  }
  else {
    // level 0 (assume this is always contained in a partial reference_key)
    index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
    r1r2_t1_coeffs[0][0][0]
      = surrData.response_function(index) * s_data_2.response_function(index);
    ++cntr;
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? reference_key[lev][1] : key[lev].size();
      r1r2_t1_coeffs[lev].resize(num_sets);
      r1r2_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = key[lev][set].size();
	RealVector& r1r2_t1_coeffs_ls = r1r2_t1_coeffs[lev][set];
	r1r2_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of R1 R2
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  r1r2_t1_coeffs_ls[pt] = surrData.response_function(index)
	    * s_data_2.response_function(index)
	    - value(surrData.continuous_variables(index), sm_mi, key,
		    r1r2_t1_coeffs, r1r2_t2_coeffs, lev-1);
	}
      }
    }
  }

}


/** Whereas expectation() supports either a reference or increment key
    (passed as generic set_partition), functions forming hierarchical
    interpolant coefficients support only a reference key (starting
    point must be set 0; end point can be controlled). */
void HierarchInterpPolyApproximation::
central_product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
			    Real mean_1, Real mean_2,
			    RealVector2DArray& cov_t1_coeffs,
			    RealMatrix2DArray& cov_t2_coeffs,
			    const UShort2DArray& reference_key)
{
  // form hierarchical t1/t2 coeffs for (R_1 - \mu_1) (R_2 - \mu_2)
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key          = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, num_levels = expansionType1Coeffs.size(),
    num_sets, num_tp_pts, cntr = 0, index;
  bool partial = !reference_key.empty();
  const SurrogateData& s_data_2 = hip_approx_2->surrData;

  cov_t1_coeffs.resize(num_levels); cov_t1_coeffs[0].resize(1);
  cov_t2_coeffs.resize(num_levels); cov_t2_coeffs[0].resize(1);
  cov_t1_coeffs[0][0].sizeUninitialized(1);
  if (data_rep->basisConfigOptions.useDerivs) {
    size_t v, num_v = sharedDataRep->numVars;
    // level 0 (assume this is always contained in a partial reference_key)
    index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
    Real data_fn1_mm1 = surrData.response_function(index) - mean_1;
    Real data_fn2_mm2 = s_data_2.response_function(index) - mean_2;
    cov_t1_coeffs[0][0][0] = data_fn1_mm1 * data_fn2_mm2;
    cov_t2_coeffs[0][0].shapeUninitialized(num_v, 1);
    Real *cov_t2_coeffs_000 = cov_t2_coeffs[0][0][0];
    const RealVector& data_grad1 = surrData.response_gradient(index);
    const RealVector& data_grad2 = s_data_2.response_gradient(index);
    for (v=0; v<num_v; ++v)
      cov_t2_coeffs_000[v]
	= data_fn1_mm1 * data_grad2[v] + data_fn2_mm2 * data_grad1[v];
    ++cntr;
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? reference_key[lev][1] : key[lev].size();
      cov_t1_coeffs[lev].resize(num_sets); cov_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	RealVector& cov_t1_coeffs_ls = cov_t1_coeffs[lev][set];
	RealMatrix& cov_t2_coeffs_ls = cov_t2_coeffs[lev][set];
	num_tp_pts = key[lev][set].size();
	cov_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	cov_t2_coeffs_ls.shapeUninitialized(num_v, num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  const RealVector& c_vars = surrData.continuous_variables(index);
	  // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  data_fn1_mm1 = surrData.response_function(index) - mean_1;
	  data_fn2_mm2 = s_data_2.response_function(index) - mean_2;
	  cov_t1_coeffs_ls[pt] = data_fn1_mm1 * data_fn2_mm2 -
	    value(c_vars, sm_mi, key, cov_t1_coeffs, cov_t2_coeffs, lev-1);
	  // type2 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	  // --> interpolated grads are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
	  Real* cov_t2_coeffs_lsp = cov_t2_coeffs_ls[pt];
	  const RealVector& data_grad1 = surrData.response_gradient(index);
	  const RealVector& data_grad2 = s_data_2.response_gradient(index);
	  const RealVector& prev_grad  = gradient_basis_variables(c_vars,
	    sm_mi, key, cov_t1_coeffs, cov_t2_coeffs, lev-1);
	  for (v=0; v<num_v; ++v)
	    cov_t2_coeffs_lsp[v] = data_fn1_mm1 * data_grad2[v]
	      + data_fn2_mm2 * data_grad1[v] - prev_grad[v];
	}
      }
    }
  }
  else {
    // level 0 (assume this is always contained in a partial reference_key)
    index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
    cov_t1_coeffs[0][0][0] = (surrData.response_function(index) - mean_1) *
                             (s_data_2.response_function(index) - mean_2);
    ++cntr;
    // levels 1:w
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = (partial) ? reference_key[lev][1] : key[lev].size();
      cov_t1_coeffs[lev].resize(num_sets); cov_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = key[lev][set].size();
	RealVector& cov_t1_coeffs_ls = cov_t1_coeffs[lev][set];
	cov_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	  cov_t1_coeffs_ls[pt]
	    = (surrData.response_function(index) - mean_1)
	    * (s_data_2.response_function(index) - mean_2)
	    - value(surrData.continuous_variables(index), sm_mi, key,
		    cov_t1_coeffs, cov_t2_coeffs, lev-1);
	}
      }
    }
  }
}


void HierarchInterpPolyApproximation::
central_product_gradient_interpolant(
  HierarchInterpPolyApproximation* hip_approx_2, Real mean_1, Real mean_2,
  const RealVector& mean1_grad, const RealVector& mean2_grad, 
  RealMatrix2DArray& cov_t1_coeff_grads, const UShort2DArray& reference_key)
{
  // form hierarchical t1 coeff grads for (R_1 - \mu_1) (R_2 - \mu_2)
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray&           sm_mi = hsg_driver->smolyak_multi_index();
  const UShort4DArray&             key = hsg_driver->collocation_key();
  const Sizet3DArray&     colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, num_levels = key.size(), num_sets, num_tp_pts, cntr = 0,
    index, v, num_deriv_vars = expansionType1CoeffGrads[0][0].numRows();
  bool partial = !reference_key.empty();
  const SurrogateData& s_data_2 = hip_approx_2->surrData;

  // level 0 (assume this is always contained in a partial reference_key)
  cov_t1_coeff_grads.resize(num_levels); cov_t1_coeff_grads[0].resize(1);
  cov_t1_coeff_grads[0][0].shapeUninitialized(num_deriv_vars, 1);
  Real* cov_t1_coeff_grads_000 = cov_t1_coeff_grads[0][0][0];
  index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
  Real r1_mm = surrData.response_function(index) - mean_1,
       r2_mm = s_data_2.response_function(index) - mean_2;
  const RealVector& r1_grad = surrData.response_gradient(index);
  const RealVector& r2_grad = s_data_2.response_gradient(index);
  for (v=0; v<num_deriv_vars; ++v)
    cov_t1_coeff_grads_000[v] = r1_mm * (r2_grad[v] - mean2_grad[v])
                              + r2_mm * (r1_grad[v] - mean1_grad[v]);
  ++cntr;
  // levels 1:w
  for (lev=1; lev<num_levels; ++lev) {
    num_sets = (partial) ? reference_key[lev][1] : key[lev].size();
    cov_t1_coeff_grads[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = key[lev][set].size();
      RealMatrix& cov_t1_coeff_grads_ls = cov_t1_coeff_grads[lev][set];
      cov_t1_coeff_grads_ls.shapeUninitialized(num_deriv_vars, num_tp_pts);
      // type1 hierarchical interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	Real r1_mm = surrData.response_function(index) - mean_1,
             r2_mm = s_data_2.response_function(index) - mean_2;
	const RealVector& r1_grad = surrData.response_gradient(index);
	const RealVector& r2_grad = s_data_2.response_gradient(index);
	const RealVector& prev_grad
	  = gradient_nonbasis_variables(surrData.continuous_variables(index),
					sm_mi, key, cov_t1_coeff_grads, lev-1);
	Real* cov_t1_coeff_grads_lsp = cov_t1_coeff_grads_ls[pt];
	for (v=0; v<num_deriv_vars; ++v)
	  cov_t1_coeff_grads_lsp[v] = r1_mm * (r2_grad[v] - mean2_grad[v])
	    + r2_mm * (r1_grad[v] - mean1_grad[v]) - prev_grad[v];
      }
    }
  }
}


void HierarchInterpPolyApproximation::
integrate_response_moments(size_t num_moments)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::integrate_response_moments()" << std::endl;
    abort_handler(-1);
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver   = data_rep->hsg_driver();
  const UShort3DArray&      sm_mi        = hsg_driver->smolyak_multi_index();
  const UShort4DArray&      key          = hsg_driver->collocation_key();
  const Sizet3DArray&       colloc_index = hsg_driver->collocation_indices();
  size_t lev, set, pt, num_levels = key.size(), num_sets, num_tp_pts,
    cntr, index, num_v = sharedDataRep->numVars;
  int m_index, moment;

  if (numericalMoments.length() != num_moments)
    numericalMoments.sizeUninitialized(num_moments);
  Real& mean = numericalMoments[0];
  mean = expectation(expansionType1Coeffs, expansionType2Coeffs);

  // size moment coefficient arrays
  RealVector2DArray mom_t1_coeffs(num_levels);
  RealMatrix2DArray mom_t2_coeffs(num_levels);
  for (lev=0; lev<num_levels; ++lev) {
    num_sets = key[lev].size();
    mom_t1_coeffs[lev].resize(num_sets);
    mom_t2_coeffs[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = key[lev][set].size();
      mom_t1_coeffs[lev][set].sizeUninitialized(num_tp_pts);
      if (data_rep->basisConfigOptions.useDerivs)
	mom_t2_coeffs[lev][set].shapeUninitialized(num_v, num_tp_pts);
    }
  }

  for (m_index=1; m_index<num_moments; ++m_index) {
    moment = m_index+1; cntr = 0;
    if (data_rep->basisConfigOptions.useDerivs) {
      size_t v;
      // level 0
      index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
      Real data_fn_mm         = surrData.response_function(index) - mean;
      mom_t1_coeffs[0][0][0]  = std::pow(data_fn_mm, moment);
      Real* mom_t2_coeffs_000 = mom_t2_coeffs[0][0][0];
      Real deriv = moment * std::pow(data_fn_mm, m_index);
      const RealVector& data_grad = surrData.response_gradient(index);
      for (v=0; v<num_v; ++v)
	mom_t2_coeffs_000[v] = deriv * data_grad[v];
      ++cntr;
      // levels 1:w
      for (lev=1; lev<num_levels; ++lev) {
	num_sets = key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  num_tp_pts = key[lev][set].size();
	  RealVector& mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  RealMatrix& mom_t2_coeffs_ls = mom_t2_coeffs[lev][set];
	  for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	    index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	    const RealVector& c_vars = surrData.continuous_variables(index);
	    // type1 hierarchical interpolation of (R - \mu)^moment
	    data_fn_mm = surrData.response_function(index) - mean;
	    mom_t1_coeffs_ls[pt] = std::pow(data_fn_mm, moment) -
	      value(c_vars, sm_mi, key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    // type2 hierarchical interpolation of (R - \mu)^moment
	    // --> interpolated grads are moment(R-\mu)^{moment-1} R'
	    Real* mom_t2_coeffs_lsp = mom_t2_coeffs_ls[pt];
	    deriv = moment * std::pow(data_fn_mm, m_index);
	    const RealVector& data_grad = surrData.response_gradient(index);
	    const RealVector& prev_grad = gradient_basis_variables(c_vars,
	      sm_mi, key, mom_t1_coeffs, mom_t2_coeffs, lev-1);
	    for (v=0; v<num_v; ++v)
	      mom_t2_coeffs_lsp[v] = deriv * data_grad[v] - prev_grad[v];
	  }
	}
      }
    }
    else {
      // level 0
      index = (colloc_index.empty()) ? cntr : colloc_index[0][0][0];
      mom_t1_coeffs[0][0][0]
	= std::pow(surrData.response_function(index) - mean, moment);
      ++cntr;
      // levels 1:w
      for (lev=1; lev<num_levels; ++lev) {
	num_sets = key[lev].size();
	for (set=0; set<num_sets; ++set) {
	  RealVector& mom_t1_coeffs_ls = mom_t1_coeffs[lev][set];
	  num_tp_pts = key[lev][set].size();
	  // type1 hierarchical interpolation of (R - \mu)^moment
	  for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	    index = (colloc_index.empty()) ? cntr : colloc_index[lev][set][pt];
	    mom_t1_coeffs_ls[pt]
	      = std::pow(surrData.response_function(index) - mean, moment)
	      - value(surrData.continuous_variables(index), sm_mi, key,
		      mom_t1_coeffs, mom_t2_coeffs, lev-1);
	  }
	}
      }
    }

    // pass these exp coefficients into a general expectation fn
    numericalMoments[m_index] = expectation(mom_t1_coeffs, mom_t2_coeffs);
  }

  // standardize third and higher central moments, if present
  //standardize_moments(numericalMoments);

  /*
  if (numericalMoments.size() != num_moments)
    numericalMoments.size(num_moments);
  if (data_rep->basisConfigOptions.useDerivs)
    integrate_moments(expansionType1Coeffs, expansionType2Coeffs,
		      hsg_driver->type1_weight_set_arrays(),
		      hsg_driver->type2_weight_set_arrays(), numericalMoments);
  else
    integrate_moments(expansionType1Coeffs,
		      hsg_driver->type1_weight_set_arrays(), numericalMoments);
  */
}


void HierarchInterpPolyApproximation::
integrate_expansion_moments(size_t num_moments)
{
  // for now: nested interpolation is exact
  expansionMoments = numericalMoments;

  // a couple different ways to go with this in the future:
  // (1) evaluate hierarchical value(lev) - value(lev-1) with HSGDriver wts
  // (2) evaluate value() with CSGDriver wts
  //  > promote Nodal implementation of this function to base class
  //  > redefine HierarchSparseGridDriver::type1_weight_sets() to generate
  //    from 1D weights array in CSG-style approach (not simple concatenation)
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
  variance = expectation(prod_member_t1_coeffs, member_t1_wts,
			 prod_member_t2_coeffs, member_t2_wts);

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
    complement_variance = expectation(cprod_member_t1_coeffs, member_t1_wts,
				      cprod_member_t2_coeffs, member_t2_wts);

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

  // Perform inner integral over complementary set u' (non-member vars) to
  // form new weighted expansion coefficients h (stored as member_t1_coeffs)
  size_t v, num_v = sharedDataRep->numVars, lev, set, pt,
    num_levels = expansionType1Coeffs.size(), num_sets, num_tp_pts,
    num_member_coeffs, v_cntr, p_cntr = 0, m_index;
  Real member_wt, nonmember_wt;
  SizetArray indexing_factor; // for indexing member coeffs,wts
  UShortArray delta_sizes; UShort2DArray delta_keys;
  member_t1_coeffs.resize(num_levels);  member_t1_wts.resize(num_levels);
  member_t2_coeffs.resize(num_levels);  member_t2_wts.resize(num_levels);
  member_colloc_key.resize(num_levels); member_colloc_index.resize(num_levels);
  for (lev=0; lev<num_levels; ++lev) {
    const RealVectorArray& t1_coeffs_l = expansionType1Coeffs[lev];
    num_sets = t1_coeffs_l.size();
    member_t1_coeffs[lev].resize(num_sets); member_t1_wts[lev].resize(num_sets);
    member_t2_coeffs[lev].resize(num_sets); member_t2_wts[lev].resize(num_sets);
    member_colloc_key[lev].resize(num_sets);
    member_colloc_index[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      const RealVector& t1_coeffs_ls = t1_coeffs_l[set];
      const RealMatrix& t2_coeffs_ls = expansionType2Coeffs[lev][set];
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
	if (data_rep->basisConfigOptions.useDerivs) {
	  m_t2_coeffs_ls.shape(num_v, num_member_coeffs);
	  m_t2_wts_ls.shape(num_v, num_member_coeffs);
	}
	m_key_ls.resize(num_member_coeffs);
	m_index_ls.resize(num_member_coeffs);
      }

      for (pt=0; pt<num_tp_pts; ++pt, ++p_cntr) {
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
	m_index_ls[m_index]  = (colloc_index.empty()) ? p_cntr :
	  colloc_index[lev][set][pt];   // links back to surrData c_vars
	m_key_ls[m_index]    = key_lsp; // links back to interp polynomials

	// now do the same for the type2 coeffs and weights
	if (data_rep->basisConfigOptions.useDerivs) {
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
	    << "\nmember_t1_coeffs[" << lev << "][" << set << "]:\n";
      write_data(PCout, member_t1_coeffs[lev][set]);
      PCout << "member_t1_wts[" << lev << "][" << set << "]:\n";
      write_data(PCout, member_t1_wts[lev][set]);
      if (data_rep->basisConfigOptions.useDerivs) {
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
  const UShort4DArray& m_colloc_key,    const Sizet3DArray&  m_colloc_index,
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

  size_t v, num_v = sharedDataRep->numVars, lev, set, pt,
    num_levels = m_t1_coeffs.size(), num_sets, num_tp_pts, index,
    h_level = num_levels - 1;
  SizetList member_indices;
  for (v=0; v<num_v; ++v)
    if (m_bits[v])
      member_indices.push_back(v);

  // form hierarchical t1/t2 coeffs for (h - mean)^2
  cprod_m_t1_coeffs.resize(num_levels); cprod_m_t1_coeffs[0].resize(1);
  cprod_m_t2_coeffs.resize(num_levels); cprod_m_t2_coeffs[0].resize(1);
  // level 0 type1
  cprod_m_t1_coeffs[0][0].sizeUninitialized(1);
  index = m_colloc_index[0][0][0];
  Real h_val_mm =
    value(surrData.continuous_variables(index), sm_mi, m_colloc_key,
	  m_t1_coeffs, m_t2_coeffs, h_level, member_indices) - mean;
  cprod_m_t1_coeffs[0][0][0] = h_val_mm * h_val_mm;
  if (data_rep->basisConfigOptions.useDerivs) {
    // level 0 type2
    const RealVector& h_grad_000 =
      gradient_basis_variables(surrData.continuous_variables(index), sm_mi,
      m_colloc_key, m_t1_coeffs, m_t2_coeffs, h_level, member_indices);
    cprod_m_t2_coeffs[0][0].shapeUninitialized(num_v, 1);
    Real *cprod_m_t2_coeffs_000 = cprod_m_t2_coeffs[0][0][0];
    for (v=0; v<num_v; ++v)
      cprod_m_t2_coeffs_000[v] = 2. * h_val_mm * h_grad_000[v];
    // levels 1:w type1 and type2
    for (lev=1; lev<num_levels; ++lev) {
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
	  index = m_colloc_index[lev][set][pt];
	  const RealVector& c_vars = surrData.continuous_variables(index);
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
	PCout << "cprod_m_t1_coeffs[" << lev << "][" << set << "]:\n";
	write_data(PCout, cprod_m_t1_coeffs[lev][set]);
	PCout << "cprod_m_t2_coeffs[" << lev << "][" << set << "]:\n";
	write_data(PCout, cprod_m_t2_coeffs[lev][set], false, true, true);
	PCout << std::endl;
#endif // VBD_DEBUG
      }
    }
  }
  else {
    // levels 1:w type1
    for (lev=1; lev<num_levels; ++lev) {
      num_sets = m_colloc_key[lev].size();
      cprod_m_t1_coeffs[lev].resize(num_sets);
      cprod_m_t2_coeffs[lev].resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = m_colloc_key[lev][set].size();
	RealVector& cprod_m_t1_coeffs_ls = cprod_m_t1_coeffs[lev][set];
	cprod_m_t1_coeffs_ls.sizeUninitialized(num_tp_pts);
	// type1 hierarchical interpolation of (h - mean)^2
	for (pt=0; pt<num_tp_pts; ++pt) {
	  index = m_colloc_index[lev][set][pt];
	  const RealVector& c_vars = surrData.continuous_variables(index);
	  h_val_mm = value(c_vars, sm_mi, m_colloc_key, m_t1_coeffs,
			   m_t2_coeffs, h_level, member_indices) - mean;
	  cprod_m_t1_coeffs_ls[pt] = h_val_mm * h_val_mm -
	    value(c_vars, sm_mi, m_colloc_key, cprod_m_t1_coeffs,
		  cprod_m_t2_coeffs, lev-1, member_indices);
	}
#ifdef VBD_DEBUG
	PCout << "cprod_m_t1_coeffs[" << lev << "][" << set << "]:\n";
	write_data(PCout, cprod_m_t1_coeffs[lev][set]);
	PCout << std::endl;
#endif // VBD_DEBUG
      }
    }
  }
}

}
