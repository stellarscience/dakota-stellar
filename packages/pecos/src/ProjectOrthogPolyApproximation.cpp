/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        ProjectOrthogPolyApproximation
//- Description:  Implementation code for ProjectOrthogPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "ProjectOrthogPolyApproximation.hpp"
#include "TensorProductDriver.hpp"
#include "IncrementalSparseGridDriver.hpp"
#include "CubatureDriver.hpp"
#include "pecos_global_defs.hpp"

//#define DEBUG

namespace Pecos {


int ProjectOrthogPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface in multiple dimensions
  return (expansionCoeffFlag || expansionCoeffGradFlag) ? 1 : 0;
}


void ProjectOrthogPolyApproximation::allocate_arrays()
{
  // SharedProjectOrthogPolyApproxData::allocate_data() has already executed

  OrthogPolyApproximation::allocate_arrays();

  // integration-specific allocations:
  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case INCREMENTAL_SPARSE_GRID:
    //if (data_rep->expConfigOptions.refineControl) {
    std::shared_ptr<IncrementalSparseGridDriver> isg_driver =
      std::static_pointer_cast<IncrementalSparseGridDriver>(data_rep->driver());
    size_t num_smolyak_indices = isg_driver->smolyak_multi_index().size();
    const ActiveKey& key = data_rep->activeKey;
    tpExpansionCoeffs[key].resize(num_smolyak_indices);
    tpExpansionCoeffGrads[key].resize(num_smolyak_indices);
    //}
    break;
  }
}


void ProjectOrthogPolyApproximation::integration_checks()
{
  if (surrData.anchor()) {
    PCerr << "Error: anchor point not supported for numerical integration in "
	  << "ProjectOrthogPolyApproximation." << std::endl;
    abort_handler(-1);
  }
  std::shared_ptr<SharedPolyApproxData> data_rep =
    std::static_pointer_cast<SharedPolyApproxData>(sharedDataRep);
  std::shared_ptr<IntegrationDriver> driver_rep = data_rep->driverRep;

  if (!driver_rep) {
    PCerr << "Error: pointer to integration driver required in "
	  << "ProjectOrthogPolyApproximation." << std::endl;
    abort_handler(-1);
  }
  size_t num_data_pts = surrData.points(),
         num_grid_pts = driver_rep->grid_size();
  if (num_data_pts != num_grid_pts) {
    PCerr << "Error: number of current points (" << num_data_pts << ") is "
	  << "not consistent with\n       number of points/weights ("
	  << num_grid_pts << ") from integration driver in\n       "
	  << "ProjectOrthogPolyApproximation." << std::endl;
    abort_handler(-1);
  }
}


void ProjectOrthogPolyApproximation::compute_coefficients()
{
  PolynomialApproximation::compute_coefficients();
  if (!expansionCoeffFlag && !expansionCoeffGradFlag)
    return;

  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
#ifdef DEBUG
  data_rep->gradient_check();
#endif // DEBUG

  // Size coefficient and Sobol arrays
  allocate_arrays();

  // calculate expansion coefficients
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case CUBATURE: // single expansion integration
    integration_checks();
    integrate_expansion(data_rep->multi_index(), surrData.variables_data(),
			surrData.response_data(),
			data_rep->driver()->type1_weight_sets(),
			expCoeffsIter->second, expCoeffGradsIter->second);
    break;
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    integration_checks();
    RealVector& exp_coeffs      =     expCoeffsIter->second;
    RealMatrix& exp_coeff_grads = expCoeffGradsIter->second;
    // multiple tensor expansion integrations
    if (expansionCoeffFlag)     exp_coeffs      = 0.;
    if (expansionCoeffGradFlag) exp_coeff_grads = 0.;
    std::shared_ptr<CombinedSparseGridDriver> csg_driver =
      std::static_pointer_cast<CombinedSparseGridDriver>(data_rep->driver());
    const IntArray&      sm_coeffs = csg_driver->smolyak_coefficients();
    const ActiveKey&   key       = data_rep->activeKey;
    const UShort3DArray& tp_mi     = data_rep->tpMultiIndex[key];
    const Sizet2DArray&  tp_mi_map = data_rep->tpMultiIndexMap[key];
    RealVectorArray& tp_exp_coeffs      = tpExpansionCoeffs[key];
    RealMatrixArray& tp_exp_coeff_grads = tpExpansionCoeffGrads[key];
    size_t i, num_tensor_grids = tp_mi.size(); int coeff;
    SDVArray tp_data_vars; SDRArray tp_data_resp;
    RealVector tp_wts, tp_coeffs; RealMatrix tp_coeff_grads;
    bool store_tp = (data_rep->expConfigOptions.refineControl);
    // loop over tensor-products, forming sub-expansions, and sum them up
    // Note: SharedOrthogPolyApproxData::allocate_data() uses
    // sparse_grid_multi_index() to build multiIndex with append_multi_index()
    for (i=0; i<num_tensor_grids; ++i) {
      // form tp_data_vars, tp_data_resp, tp_wts using collocKey et al.
      integration_data(i, tp_data_vars, tp_data_resp, tp_wts);

      // form tp_multi_index from tpMultiIndexMap
      //map_tensor_product_multi_index(tp_multi_index, i);

      // form tp expansion coeffs
      RealVector& tp_coeffs_i = (store_tp) ? tp_exp_coeffs[i] : tp_coeffs;
      RealMatrix& tp_grads_i  = (store_tp) ?
	tp_exp_coeff_grads[i] : tp_coeff_grads;
      integrate_expansion(tp_mi[i], tp_data_vars, tp_data_resp, tp_wts,
			  tp_coeffs_i, tp_grads_i);

      // sum tensor product coeffs/grads into expansion coeffs/grads
      coeff = sm_coeffs[i];
      if (coeff)
	overlay_expansion(tp_mi_map[i], tp_coeffs_i, tp_grads_i, coeff,
			  exp_coeffs, exp_coeff_grads);
    }
    break;
  }
  case SAMPLING:
    surrData.data_checks(); // defines failed resp map
    expectation();
    break;
  default:
    PCerr << "Error: unsupported expCoeffsSolnApproach in ProjectOrthogPoly"
	  << "Approximation::compute_coefficients()" << std::endl;
    abort_handler(-1);
    break;
  }

  clear_computed_bits();
}


void ProjectOrthogPolyApproximation::increment_coefficients()
{
  // TO DO: partial sync of new TP data set, e.g. update_surrogate_data() ?
  synchronize_surrogate_data();

  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
  const ActiveKey& key = data_rep->activeKey;

  // synchronize expansionCoeff{s,Grads} and approxData
  update_active_iterators(key);
  // resize component Sobol' array sizes to pick up new interaction terms
  // (based on size of data_rep->sobolIndexMap)
  allocate_component_sobol();

  // for use in pop_coefficients()
  prevExpCoeffs     = expCoeffsIter->second;     // copy
  prevExpCoeffGrads = expCoeffGradsIter->second; // copy

  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case INCREMENTAL_SPARSE_GRID: {
    // tpMultiIndex{,Map,MapRef} already updated in
    // SharedProjectOrthogPolyApproxData::increment_data()
    const UShort3DArray& tp_mi = data_rep->tpMultiIndex[key];
    RealVectorArray& tp_exp_coeffs      = tpExpansionCoeffs[key];
    RealMatrixArray& tp_exp_coeff_grads = tpExpansionCoeffGrads[key];
    size_t start_append = tp_exp_coeffs.size();
    SDVArray tp_data_vars; SDRArray tp_data_resp; RealVector tp_wts;

    switch (data_rep->expConfigOptions.refineControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
      RealVector rv; tp_exp_coeffs.push_back(rv);
      RealMatrix rm; tp_exp_coeff_grads.push_back(rm);

      // form tp_data_pts, tp_wts using collocKey et al.
      integration_data(start_append, tp_data_vars, tp_data_resp, tp_wts);
      // form trial expansion coeffs/grads
      integrate_expansion(tp_mi[start_append], tp_data_vars, tp_data_resp,
			  tp_wts, tp_exp_coeffs[start_append],
			  tp_exp_coeff_grads[start_append]);
      break;
    }
    default: { // multiple index sets from iso/aniso sparse grid refinement
      size_t i, num_tp_mi = tp_mi.size();
      tp_exp_coeffs.resize(num_tp_mi);  tp_exp_coeff_grads.resize(num_tp_mi);
      // loop over tensor-products, forming sub-expansions, and sum them up
      for (i=start_append; i<num_tp_mi; ++i) {
	// form tp_data_vars, tp_data_resp, tp_wts using collocKey et al.
	integration_data(i, tp_data_vars, tp_data_resp, tp_wts);
	// form tp expansion coeffs
	integrate_expansion(tp_mi[i], tp_data_vars, tp_data_resp, tp_wts,
			    tp_exp_coeffs[i], tp_exp_coeff_grads[i]);
      }
      break;
    }
    }

    // sum trial expansion into expansionCoeffs/expansionCoeffGrads
    append_tensor_expansions(start_append);
    break;
  }
  case QUADRATURE: case CUBATURE:
    integration_checks();
    integrate_expansion(data_rep->multi_index(), surrData.variables_data(),
			surrData.response_data(),
			data_rep->driver()->type1_weight_sets(),
			expCoeffsIter->second, expCoeffGradsIter->second);
    break;
  }

  clear_computed_bits();
}


void ProjectOrthogPolyApproximation::pop_coefficients(bool save_data)
{
  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
  const ActiveKey& key = data_rep->activeKey;

  // likely overkill, but multilevel roll up after increment modifies and
  // then restores active key
  update_active_iterators(key);

  if (save_data)
    switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
    case INCREMENTAL_SPARSE_GRID: {
      RealVectorArray& tp_exp_coeffs      = tpExpansionCoeffs[key];
      RealMatrixArray& tp_exp_coeff_grads = tpExpansionCoeffGrads[key];

      switch (data_rep->expConfigOptions.refineControl) {
      case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED:
	// reset tensor-product bookkeeping and save restorable data
	poppedExpCoeffs[key].push_back(tp_exp_coeffs.back());
	poppedExpCoeffGrads[key].push_back(tp_exp_coeff_grads.back());
	tp_exp_coeffs.pop_back();  tp_exp_coeff_grads.pop_back();
	break;
      default: { // multiple index sets from iso/aniso sparse grid refinement
	const UShort3DArray& tp_mi = data_rep->tpMultiIndex[key];
	RealVectorDeque& pop_exp_coeffs      = poppedExpCoeffs[key];
	RealMatrixDeque& pop_exp_coeff_grads = poppedExpCoeffGrads[key];
	size_t i, num_tp_mi = tp_mi.size(), num_tp_exp = tp_exp_coeffs.size();
	RealVectorArray::iterator tp_ec_it  = tp_exp_coeffs.begin();
	RealMatrixArray::iterator tp_ecg_it = tp_exp_coeff_grads.begin();
	std::advance(tp_ec_it, num_tp_mi);  std::advance(tp_ecg_it, num_tp_mi);
	pop_exp_coeffs.insert(pop_exp_coeffs.end(), tp_ec_it,
			      tp_exp_coeffs.end());
	pop_exp_coeff_grads.insert(pop_exp_coeff_grads.end(), tp_ecg_it,
				   tp_exp_coeff_grads.end());
	tp_exp_coeffs.resize(num_tp_mi);  tp_exp_coeff_grads.resize(num_tp_mi);
	break;
      }
      }
      break;
    }
    case QUADRATURE: case CUBATURE:
      poppedExpCoeffs[key].push_back(expCoeffsIter->second);
      poppedExpCoeffGrads[key].push_back(expCoeffGradsIter->second);
      break;
    }

  // reset expansion{Coeffs,CoeffGrads}
  expCoeffsIter->second     = prevExpCoeffs;
  expCoeffGradsIter->second = prevExpCoeffGrads;
  // don't update Sobol' array sizes for decrement, push, or finalize

  clear_computed_bits(); // previous moments are not stored
}


void ProjectOrthogPolyApproximation::push_coefficients()
{
  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
  const ActiveKey& key = data_rep->activeKey;

  // synchronize expansionCoeff{s,Grads} and approxData
  update_active_iterators(key);

  // for use in decrement; both pushes and new increments can be popped
  prevExpCoeffs     = expCoeffsIter->second;     // copy
  prevExpCoeffGrads = expCoeffGradsIter->second; // copy

  RealVectorDeque& pop_exp_coeffs      = poppedExpCoeffs[key];
  RealMatrixDeque& pop_exp_coeff_grads = poppedExpCoeffGrads[key];
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case INCREMENTAL_SPARSE_GRID: {
    RealVectorArray& tp_exp_coeffs      = tpExpansionCoeffs[key];
    RealMatrixArray& tp_exp_coeff_grads = tpExpansionCoeffGrads[key];
    size_t start_append = tp_exp_coeffs.size(); // before push_back

    switch (data_rep->expConfigOptions.refineControl) {
    case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
      // move previous expansion data to current expansion
      size_t index_star = data_rep->push_index();

      RealVectorDeque::iterator cit = pop_exp_coeffs.begin();
      RealMatrixDeque::iterator git = pop_exp_coeff_grads.begin();
      std::advance(cit, index_star); std::advance(git, index_star);

      tp_exp_coeffs.push_back(*cit);
      tp_exp_coeff_grads.push_back(*git);
      pop_exp_coeffs.erase(cit); pop_exp_coeff_grads.erase(git);
      break;
    }
    default: // multiple index sets from one sparse grid refinement candidate
      tp_exp_coeffs.insert(tp_exp_coeffs.end(), pop_exp_coeffs.begin(),
			   pop_exp_coeffs.end());
      tp_exp_coeff_grads.insert(tp_exp_coeff_grads.end(),
				pop_exp_coeff_grads.begin(),
				pop_exp_coeff_grads.end());
      pop_exp_coeffs.clear();  pop_exp_coeff_grads.clear();
      break;
    }

    // sum trial expansion into expansionCoeffs/expansionCoeffGrads
    append_tensor_expansions(start_append);
    break;
  }
  case QUADRATURE: case CUBATURE:
    expCoeffsIter->second     = pop_exp_coeffs.back();
    expCoeffGradsIter->second = pop_exp_coeff_grads.back();
    pop_exp_coeffs.pop_back();  pop_exp_coeff_grads.pop_back();
    break;    
  }

  // don't update Sobol' array sizes for decrement, push, or finalize

  clear_computed_bits();
}


void ProjectOrthogPolyApproximation::finalize_coefficients()
{
  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
  const ActiveKey& key = data_rep->activeKey;

  // synchronize expansionCoeff{s,Grads} and approxData
  update_active_iterators(key);

  RealVectorDeque& pop_exp_coeffs      = poppedExpCoeffs[key];
  RealMatrixDeque& pop_exp_coeff_grads = poppedExpCoeffGrads[key];
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case INCREMENTAL_SPARSE_GRID: {

    // Note: finalization fns only used for generalized sparse grids, but
    // would be the same for a single iso/aniso refinement candidate with
    // multiple index sets

    RealVectorArray& tp_exp_coeffs      = tpExpansionCoeffs[key];
    RealMatrixArray& tp_exp_coeff_grads = tpExpansionCoeffGrads[key];
    // don't update Sobol' array sizes for decrement, push, or finalize
    size_t start_append = tp_exp_coeffs.size(); // before insertion
    // move previous expansion data to current expansion

    // Note: popped sets are not explicitly added in computed_trial_sets()
    //       order as in IncrementalSparseGridDriver::finalize_sets().
    //       However, poppedLevMultiIndex et al. become ordered due to
    //       enumeration of ordered active_multi_index().  Rather than
    //       incurring additional overhead by mapping indices, a compile-time
    //       verification block is defined in shared pre_finalize_data().
    tp_exp_coeffs.insert(tp_exp_coeffs.end(), pop_exp_coeffs.begin(),
			 pop_exp_coeffs.end());
    tp_exp_coeff_grads.insert(tp_exp_coeff_grads.end(),
			      pop_exp_coeff_grads.begin(),
			      pop_exp_coeff_grads.end());
    // sum remaining trial expansions into expansionCoeff{s,Grads}.  For
    // finalize, don't need to cache prevExpCoeff{s,Grads} prior to append.
    append_tensor_expansions(start_append);
    break;
  }
  case QUADRATURE: case CUBATURE: // for completeness (not used)
    if (!pop_exp_coeffs.empty()) expCoeffsIter->second = pop_exp_coeffs.back();
    if (!pop_exp_coeff_grads.empty())
      expCoeffGradsIter->second = pop_exp_coeff_grads.back();
    break;
  }

  pop_exp_coeffs.clear();  pop_exp_coeff_grads.clear();
  clear_computed_bits();
}


void ProjectOrthogPolyApproximation::
append_tensor_expansions(size_t start_tp_index)
{
  // synchonize expansionCoeff{s,Grads} size with updated multiIndex
  // (following any caching of previous states)
  resize_expansion();

  // update expansion{Coeffs,CoeffGrads} using an incremental update
  // rather than building from scratch
  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
  std::shared_ptr<IncrementalSparseGridDriver> isg_driver =
    std::static_pointer_cast<IncrementalSparseGridDriver>(data_rep->driver());
  const IntArray&     sm_coeffs = isg_driver->smolyak_coefficients();
  const IntArray& sm_coeffs_ref = isg_driver->smolyak_coefficients_reference();
#ifdef DEBUG
  PCout << "In ProjectOrthogPolyApproximation::append_tensor_expansions() with "
	<< "start index " << start_tp_index << "\nsm_coeffs:\n" << sm_coeffs
	<< "sm_coeffs_ref:\n" << sm_coeffs_ref << std::endl;
#endif // DEBUG
  const ActiveKey&  key       = data_rep->activeKey;
  const Sizet2DArray& tp_mi_map = data_rep->tpMultiIndexMap[key];
  RealVectorArray& tp_exp_coeffs      = tpExpansionCoeffs[key];
  RealMatrixArray& tp_exp_coeff_grads = tpExpansionCoeffGrads[key];

  // add trial expansions
  size_t index, num_tensor_grids = sm_coeffs.size();
  int coeff, delta_coeff;
  RealVector& exp_coeffs      =     expCoeffsIter->second;
  RealMatrix& exp_coeff_grads = expCoeffGradsIter->second;
  for (index=start_tp_index; index<num_tensor_grids; ++index) {
    coeff = sm_coeffs[index];
    if (coeff)
      overlay_expansion(tp_mi_map[index], tp_exp_coeffs[index],
			tp_exp_coeff_grads[index], coeff, exp_coeffs,
			exp_coeff_grads);
#ifdef DEBUG
    PCout << "Trial set sm_coeff = "  << coeff
	  << "\ntpExpansionCoeffs:\n" << tp_exp_coeffs[index]
	  << "\ntpMultiIndexMap:\n"   << tp_mi_map[index] << '\n';
#endif // DEBUG
  }
  // update other expansion contributions with a changed smolyak coefficient
  for (index=0; index<start_tp_index; ++index) {
    // add new, subtract previous
    delta_coeff = sm_coeffs[index] - sm_coeffs_ref[index];
#ifdef DEBUG
    PCout << "Old set delta_coeff = " << delta_coeff
	  << "\ntpExpansionCoeffs:\n" << tp_exp_coeffs[index]
	  << "\ntpMultiIndexMap:\n"   << tp_mi_map[index] << '\n';
#endif // DEBUG
    if (delta_coeff)
      overlay_expansion(tp_mi_map[index], tp_exp_coeffs[index],
			tp_exp_coeff_grads[index], delta_coeff, exp_coeffs,
			exp_coeff_grads);
  }
}


void ProjectOrthogPolyApproximation::
integration_data(size_t tp_index, SDVArray& tp_data_vars,
		 SDRArray& tp_data_resp, RealVector& tp_weights)
{
  // extract tensor vars/resp from surrData and tensor wts from
  // type1CollocWts1D
  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
  std::shared_ptr<CombinedSparseGridDriver> csg_driver =
    std::static_pointer_cast<CombinedSparseGridDriver>(data_rep->driver());
  const UShortArray&    sm_index = csg_driver->smolyak_multi_index()[tp_index];
  const UShort2DArray&       key = csg_driver->collocation_key()[tp_index];
  const SizetArray&  colloc_index = csg_driver->collocation_indices()[tp_index];
  const Real3DArray& colloc_wts_1d = csg_driver->type1_collocation_weights_1d();
  const SDVArray& data_vars = surrData.variables_data();
  const SDRArray& data_resp = surrData.response_data();
  size_t i, j, index, num_tp_pts = colloc_index.size(),
    num_v = sharedDataRep->numVars;
  tp_data_vars.resize(num_tp_pts); tp_data_resp.resize(num_tp_pts);
  tp_weights.resize(num_tp_pts);
  for (i=0; i<num_tp_pts; ++i) {
    // tensor-product vars/resp
    index = colloc_index[i];
    tp_data_vars[i] = data_vars[index];
    tp_data_resp[i] = data_resp[index];
    // tensor-product weight
    Real& tp_wts_i = tp_weights[i]; tp_wts_i = 1.;
    const UShortArray& key_i = key[i];
    for (j=0; j<num_v; ++j)
      tp_wts_i *= colloc_wts_1d[sm_index[j]][j][key_i[j]];
  }
}


/** The coefficients of the PCE for the response are calculated using a
    spectral projection of the response against each multivariate orthogonal
    polynomial basis fn using the inner product ratio <f,Psi>/<Psi^2>, where
    inner product <a,b> is the n-dimensional integral of a*b*weighting over
    the support range of the n-dimensional (composite) weighting function.
    1-D quadrature rules are defined for specific 1-D weighting functions
    and support ranges and approximate the integral of f*weighting as the
    Sum_i of w_i f_i.  To extend this to n-dimensions, a tensor product
    quadrature rule, cubature, or Smolyak sparse grid rule is applied.  
    It is not necessary to approximate the integral for the denominator
    numerically, since this is available analytically. */
void ProjectOrthogPolyApproximation::
integrate_expansion(const UShort2DArray& multi_index,
		    const SDVArray& data_vars, const SDRArray& data_resp,
		    const RealVector& wt_sets, RealVector& exp_coeffs,
		    RealMatrix& exp_coeff_grads)
{
  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);

  // Perform numerical integration via tensor-product quadrature/cubature/
  // Smolyak sparse grids.  Quadrature/cubature use a single application of
  // point and weight sets computed by TensorProductDriver/CubatureDriver, and
  // sparse grids could do this as well, but it is better to integrate the
  // sparse grid on a per-tensor-product basis folowed by summing the
  // corresponding PC expansions.
  if (data_resp[0].is_null()) {
    PCerr << "Error: null SDR in ProjectOrthogPolyApproximation::"
	  << "integrate_expansion()" << std::endl;
    abort_handler(-1);
  }
  size_t i, j, k, num_exp_terms = multi_index.size(),
    num_pts = std::min(data_vars.size(), data_resp.size()),
    num_deriv_vars = data_resp[0].response_gradient().length();
  Real wt_resp_fn_i, Psi_ij; Real* exp_grad;
  RealVector wt_resp_grad_i;
  if (expansionCoeffFlag) { // shape if needed and zero out
    if (exp_coeffs.length() != num_exp_terms)
      exp_coeffs.size(num_exp_terms); // init to 0
    else
      exp_coeffs = 0.;
  }
  if (expansionCoeffGradFlag) {
    if (exp_coeff_grads.numRows() != num_deriv_vars ||
	exp_coeff_grads.numCols() != num_exp_terms)
      exp_coeff_grads.shape(num_deriv_vars, num_exp_terms); // init to 0
    else
      exp_coeff_grads = 0.;
    wt_resp_grad_i.sizeUninitialized(num_deriv_vars);
  }
  for (i=0; i<num_pts; ++i) {
    if (expansionCoeffFlag)
      wt_resp_fn_i = wt_sets[i] * data_resp[i].response_function();
    if (expansionCoeffGradFlag) {
      wt_resp_grad_i = data_resp[i].response_gradient(); // copy
      wt_resp_grad_i.scale(wt_sets[i]);
    }
#ifdef DEBUG
    PCout << "wt = " << wt_sets[i] << " resp = "
	  << data_resp[i].response_function() << std::endl;
#endif //DEBUG
    const RealVector& c_vars_i = data_vars[i].continuous_variables();
    for (j=0; j<num_exp_terms; ++j) {
      Psi_ij = data_rep->multivariate_polynomial(c_vars_i, multi_index[j]);
      if (expansionCoeffFlag) {
	exp_coeffs[j] += Psi_ij * wt_resp_fn_i;
#ifdef DEBUG
	PCout << "Psi[" << i << "][" << j << "] = " << Psi_ij
	      << " exp_coeffs[" << j << "] = " << exp_coeffs[j] << std::endl;
#endif //DEBUG
      }
      if (expansionCoeffGradFlag) {
	exp_grad = exp_coeff_grads[j];
	for (k=0; k<num_deriv_vars; ++k)
	  exp_grad[k] += Psi_ij * wt_resp_grad_i[k];
      }
    }
  }

  for (i=0; i<num_exp_terms; ++i) {
    Real norm_sq = data_rep->norm_squared(multi_index[i]);
    if (expansionCoeffFlag)
      exp_coeffs[i] /= norm_sq;
    if (expansionCoeffGradFlag) {
      exp_grad = exp_coeff_grads[i];
      for (k=0; k<num_deriv_vars; ++k)
	exp_grad[k] /= norm_sq;
    }
  }
#ifdef DEBUG
  PCout << "expansion_coeffs:\n" << exp_coeffs;
  if (exp_coeff_grads.numRows())
    PCout << "expansion_coeff_grads:\n" << exp_coeff_grads;
  PCout << "\n\n";
#endif // DEBUG
}


/** The coefficients of the PCE for the response are calculated using a
    spectral projection of the response against each multivariate orthogonal
    polynomial basis fn using the inner product ratio <f,Psi>/<Psi^2>,
    where inner product <a,b> is the n-dimensional integral of a*b*weighting
    over the support range of the n-dimensional (composite) weighting
    function.  When interpreting the weighting function as a probability
    density function, <a,b> = expected value of a*b, which can be evaluated
    by sampling from the probability density function and computing the mean
    statistic.  It is not necessary to compute the mean statistic for the
    denominator, since this is available analytically. */
void ProjectOrthogPolyApproximation::expectation()
{
  RealVector& exp_coeffs      = expCoeffsIter->second;
  RealMatrix& exp_coeff_grads = expCoeffGradsIter->second;

  const SDVArray& sdv_array   = surrData.variables_data();
  const SDRArray& sdr_array   = surrData.response_data();

  // "lhs" or "random", no weights needed
  size_t i, j, k, num_deriv_vars = exp_coeff_grads.numRows(),
    num_surr_data_pts = surrData.points(), num_failed_surr_fn = 0,
    num_failed_surr_grad = 0;
  SizetShortMap::const_iterator fit;
  const SizetShortMap& failed_resp_data = surrData.failed_response_data();
  for (fit=failed_resp_data.begin(); fit!=failed_resp_data.end(); ++fit) {
    if (fit->second & 1) ++num_failed_surr_fn;
    if (fit->second & 2) ++num_failed_surr_grad;
  }
  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
  const UShort2DArray& mi = data_rep->multi_index();
  size_t num_data_pts_fn = num_surr_data_pts - num_failed_surr_fn,
    num_data_pts_grad    = num_surr_data_pts - num_failed_surr_grad,
    num_exp_terms = mi.size();
  if (expansionCoeffFlag)
    PCout << "Expectations of " << num_exp_terms << " chaos coefficients "
	  << "using " << num_data_pts_fn << " observations.\n";
  if (expansionCoeffGradFlag)
    PCout << "Expectations of gradients of " << num_exp_terms << " chaos "
	  << "coefficients using " << num_data_pts_grad << " observations.\n";

  /*
  // The following implementation evaluates all PCE coefficients
  // using a consistent expectation formulation
  for (i=0; i<num_exp_terms; ++i) {
    Real& exp_coeff_i = exp_coeffs[i]; exp_coeff_i = 0.;
    const UShortArray& mi_i = mi[i];
    for (j=0; j<num_data_pts; ++j)
      exp_coeff_i += sdr_array[j].response_function() * data_rep->
        multivariate_polynomial(sdv_array[j].continuous_variables(), mi_i);
    exp_coeff_i /= num_data_pts * data_rep->norm_squared(mi_i);
#ifdef DEBUG
    PCout << "coeff[" << i << "] = " << exp_coeff_i << " norm squared[" << i
          << "] = " << data_rep->norm_squared(mi_i) << '\n';
#endif // DEBUG
  }
  */

  // This alternate implementation evaluates the first PCE coefficient (the
  // response mean) as an expectation and then removes the mean from the
  // expectation evaluation of all subsequent coefficients.  This approach
  // has been observed to result in better results for small sample sizes.
  Real empty_r;
  Real& mean      = (expansionCoeffFlag)     ? exp_coeffs[0] : empty_r;
  Real* mean_grad = (expansionCoeffGradFlag) ? exp_coeff_grads[0] : NULL;
  if (expansionCoeffFlag)     exp_coeffs      = 0.;
  if (expansionCoeffGradFlag) exp_coeff_grads = 0.;
  for (k=0, fit=failed_resp_data.begin(); k<num_surr_data_pts; ++k) {
    bool add_val = expansionCoeffFlag, add_grad = expansionCoeffGradFlag;
    fail_booleans(fit, k, add_val, add_grad);
    const SurrogateDataResp& sdr = sdr_array[k];
    if (add_val)
      mean += sdr.response_function();
    if (add_grad) {
      const RealVector& curr_pt_grad = sdr.response_gradient();
      for (j=0; j<num_deriv_vars; ++j)
	mean_grad[j] += curr_pt_grad[j];
    }
  }
  if (expansionCoeffFlag)
    mean /= num_data_pts_fn;
  if (expansionCoeffGradFlag)
    for (j=0; j<num_deriv_vars; ++j)
      mean_grad[j] /= num_data_pts_grad;

  Real chaos_sample, resp_fn_minus_mean, norm_sq; Real* exp_grad_i;
  RealVector resp_grad_minus_mean;
  if (expansionCoeffGradFlag)
    resp_grad_minus_mean.sizeUninitialized(num_deriv_vars);
  for (k=0, fit=failed_resp_data.begin(); k<num_surr_data_pts; ++k) {
    bool add_val = expansionCoeffFlag, add_grad = expansionCoeffGradFlag;
    fail_booleans(fit, k, add_val, add_grad);
    const SurrogateDataResp& sdr = sdr_array[k];
    if (add_val)
      resp_fn_minus_mean = sdr.response_function() - mean;
    if (add_grad) {
      const RealVector& resp_grad = sdr.response_gradient();
      for (j=0; j<num_deriv_vars; ++j)
	resp_grad_minus_mean[j] = resp_grad[j] - mean_grad[j];
    }
    const RealVector& c_vars = sdv_array[k].continuous_variables();
    for (i=1; i<num_exp_terms; ++i) {
      chaos_sample = data_rep->multivariate_polynomial(c_vars, mi[i]);
      if (add_val)
	exp_coeffs[i] += resp_fn_minus_mean * chaos_sample;
      if (add_grad) {
	exp_grad_i = exp_coeff_grads[i];
	for (j=0; j<num_deriv_vars; ++j)
	  exp_grad_i[j] += resp_grad_minus_mean[j] * chaos_sample;
      }
    }
  }
  for (i=1; i<num_exp_terms; ++i) {
    norm_sq = data_rep->norm_squared(mi[i]);
    if (expansionCoeffFlag)
      exp_coeffs[i] /= norm_sq * num_data_pts_fn;
    if (expansionCoeffGradFlag) {
      exp_grad_i = exp_coeff_grads[i];
      for (j=0; j<num_deriv_vars; ++j)
	exp_grad_i[j] /= norm_sq * num_data_pts_grad;
    }
#ifdef DEBUG
    PCout << "coeff[" << i << "] = " << exp_coeffs[i]
        //<< "coeff_grad[" << i <<"] = " << exp_grad_i
	  << " norm squared[" << i <<"] = " << norm_sq << '\n';
#endif // DEBUG
  }
}


void ProjectOrthogPolyApproximation::
synthetic_surrogate_data(SurrogateData& surr_data)
{
  // Update the active key of surr_data with synthetic data based on the
  // active grid from csg_driver

  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);

  // CombinedSparseGridDriver::combined_to_active() transfers all data except
  // collocation indices, which are invalidated by the combination.  In support
  // of the synthetic data to be created, new colloc indices are defined at the
  // end of CSGDriver::combined_to_active(), and ordering is preserved below. 
  const RealMatrix&   var_sets = data_rep->driver()->variable_sets();
  const RealVector& exp_coeffs = expCoeffsIter->second;

  // shallow copy vars array data using CombinedSparseGridDriver::variableSets
  surr_data.clear_all_active();
  size_t num_v = var_sets.numRows(), num_pts = var_sets.numCols();
  surr_data.resize(num_pts, 1, num_v); // no deriv data
  
  // use interpolant to produce data values that are being interpolated
  // Note: SurrogateData reconstruction follows order of unique variable sets
  SDVArray& sdv_array = surr_data.variables_data(); 
  SDRArray& sdr_array = surr_data.response_data();
  for (size_t pt=0; pt<num_pts; ++pt) {
    RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[pt]),
		      (int)num_v);
    sdv_array[pt].continuous_variables(c_vars);
    sdr_array[pt].response_function(value(c_vars));
  }
}


void ProjectOrthogPolyApproximation::
integrate_response_moments(size_t num_moments)//, bool combined_stats)
{
  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);

  // use surrData: original single-level or synthetic combined
  const SDRArray& sdr_array = surrData.response_data();
  size_t i, num_pts = sdr_array.size();
  RealVector data_fns(num_pts);
  for (i=0; i<num_pts; ++i)
    data_fns[i] = sdr_array[i].response_function();

  /*
  // For MLMF, augment surrData using evaluations from stored expansions.
  // This fn is used in final expansion post-processing so combined_to_active()
  // has been applied.
  //
  // Note: this is a mixed implementation --> mixes direct response integration
  // (using surrData for maximal grid / active expansion, ignoring combined
  // coefficients) with surrogate integration (prior to clear_inactive(),
  // interpolate uncombined/non-active expansions onto maximal grid since
  // these response data aren't available for non-maximal surrData keys).
  // > similar to NodalInterpPolyApproximation::combine_coefficients(), which
  //   is then used in NIPA::integrate_response_moments()
  // > consider deactivating response moment integration for ML expansions
  //   (the default at OPA used by ROPA) since some of the intent (a numerical
  //   integration sanity check) is lost.  However, a useful benefit is it
  //   supports higher moment estimates.
  // > consider just evaluating the active expansion (promoted from combined)
  //   on the maximal grid (combined_stats is then easily supported, swapping
  //   expansionCoeffs and combinedExpCoeffs).  But numerical integration of
  //   the ML expansion doesn't provide much value beyond analytic moments of
  //   the same expansion --> if preserved, the current implementation seems
  //   the best compromise, short of the next item:
  // > the only rigorous way to do this would be to roll up the numerical
  //   quadrature contributions at each level --> perhaps NIPA could be
  //   (improved and) leveraged for this in the future.
  const std::map<ActiveKey, UShort2DArray>& mi = data_rep->multiIndex;
  if (mi.size() > 1) { //if (combined_stats) {
    std::map<ActiveKey, UShort2DArray>::const_iterator mi_cit = mi.begin();
    std::map<ActiveKey, RealVector>::const_iterator ec_cit;
    short combine_type = data_rep->expConfigOptions.combineType;
    const SDVArray& sdv_array = surrData.variables_data();
    for (ec_cit = expansionCoeffs.begin(); ec_cit != expansionCoeffs.end();
	 ++ec_cit, ++mi_cit)
      if (ec_cit != expCoeffsIter) {
	const UShort2DArray& mi_i = mi_cit->second;
	const RealVector&    ec_i = ec_cit->second;
	switch (combine_type) {
	case MULT_COMBINE:
	  for (i=0; i<num_pts; ++i)
	    data_fns[i] *= OrthogPolyApproximation::
	      value(sdv_array[i].continuous_variables(), mi_i, ec_i);
	  break;
	default: //case ADD_COMBINE: (correction specification not required)
	  for (i=0; i<num_pts; ++i)
	    data_fns[i] += OrthogPolyApproximation::
	      value(sdv_array[i].continuous_variables(), mi_i, ec_i);
	  break;
	}
      }
  }
  */

  // update numerical moments (secondaryMoments) based on data_fns
  if (secondaryMoments.length() != num_moments)
    secondaryMoments.sizeUninitialized(num_moments);
  integrate_moments(data_fns, data_rep->driverRep->type1_weight_sets(),
		    secondaryMoments);
}


Real ProjectOrthogPolyApproximation::value(const RealVector& x)
{
  // sum expansion to get response value prediction

  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    if (data_rep->expConfigOptions.combineType) // not guaranteed to use
      return OrthogPolyApproximation::value(x); // tensor indexing
    else { // Horner's rule approach applicable for tensor indexing
      if (!expansionCoeffFlag) { // check for required data
	PCerr << "Error: expansion coefficients not defined in "
	      << "ProjectOrthogPolyApproximation::value()" << std::endl;
	abort_handler(-1);
      }
      RealVector accumulator(sharedDataRep->numVars); // init to 0.
      return data_rep->
	tensor_product_value(x, expansionCoeffs[data_rep->activeKey],
			     data_rep->expansion_order(),
			     data_rep->multi_index(), accumulator);
    }
    break;
  /*
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    // Horner's rule approach requires storage of tpExpansionCoeffs in
    // compute_coefficients().  For now, leave store_tp as is and use
    // default approach if tpExpansionCoeffs is empty.  In addition,
    // tp arrays are not currently updated for expansion combinations.
    if (tpExpansionCoeffs.empty() || data_rep->expConfigOptions.combineType)
      // most cases
      return OrthogPolyApproximation::value(x);
    else { // generalized sparse grid case
      // Error check for required data
      if (!expansionCoeffFlag) {
	PCerr << "Error: expansion coefficients not defined in "
	      << "ProjectOrthogPolyApproximation::value()" << std::endl;
	abort_handler(-1);
      }
      CombinedSparseGridDriver* csg_driver
        = (CombinedSparseGridDriver*)data_rep->driver();
      const UShort2DArray& sm_mi     = csg_driver->smolyak_multi_index();
      const IntArray&      sm_coeffs = csg_driver->smolyak_coefficients();
      RealVector accumulator(sharedDataRep->numVars); // init to 0.
      Real approx_val = 0.;
      size_t i, num_sm_mi = sm_mi.size(); int sm_coeff;
      for (i=0; i<num_sm_mi; ++i) {
	sm_coeff = sm_coeffs[i];
	if (sm_coeff)
	  approx_val += sm_coeff * data_rep->
	    tensor_product_value(x, tp_exp_coeffs[i], tpApproxOrders[i],// TO DO
				 tp_mi[i], accumulator);
      }
      return approx_val;
    }
    break;
  }
  */
  default: // other cases are total-order expansions
    return OrthogPolyApproximation::value(x);
    break;
  }
}


Real ProjectOrthogPolyApproximation::
stored_value(const RealVector& x, const ActiveKey& key)
{
  // sum expansion to get response value prediction

  std::shared_ptr<SharedProjectOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedProjectOrthogPolyApproxData>(sharedDataRep);
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: { // Horner's rule approach
    // Note: requires tensor indexing in each multiIndex
    RealVector accumulator(sharedDataRep->numVars); // init to 0.
    return data_rep->
      tensor_product_value(x, expansionCoeffs[key],
			   data_rep->keyed_expansion_order(key),
			   data_rep->multi_index(key), accumulator);
    break;
  }
  // Horner's rule approach would require storage of tensor product components
  //case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID:
    //break;
  default: // other cases are total-order expansions
    return OrthogPolyApproximation::stored_value(x, key);
    break;
  }
}

} // namespace Pecos
