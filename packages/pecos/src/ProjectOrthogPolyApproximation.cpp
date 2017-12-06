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
#include "SharedProjectOrthogPolyApproxData.hpp"
#include "TensorProductDriver.hpp"
#include "CombinedSparseGridDriver.hpp"
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
  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case COMBINED_SPARSE_GRID:
    if (data_rep->expConfigOptions.refinementControl ==
	DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
      CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
      size_t num_smolyak_indices = csg_driver->smolyak_multi_index().size();
      tpExpansionCoeffs.resize(num_smolyak_indices);
      tpExpansionCoeffGrads.resize(num_smolyak_indices);
    }
    break;
  }
}


void ProjectOrthogPolyApproximation::compute_coefficients()
{
  if (!expansionCoeffFlag && !expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "ProjectOrthogPolyApproximation::compute_coefficients().\n"
	  << "         Bypassing approximation construction." << std::endl;
    return;
  }

  // For testing of anchor point logic:
  //size_t index = surrData.points() - 1;
  //surrData.anchor_point(surrData.variables_data()[index],
  //                      surrData.response_data()[index]);
  //surrData.pop(1);

  // anchor point, if present, is handled differently for different
  // expCoeffsSolnApproach settings:
  //   SAMPLING:   treat it as another data point
  //   QUADRATURE/CUBATURE/COMBINED_SPARSE_GRID: error
  //   LEAST_SQ_REGRESSION: use equality-constrained least squares
  size_t i, j, num_total_pts = surrData.points(),
    num_v = sharedDataRep->numVars;
  if (surrData.anchor())
    ++num_total_pts;
  if (!num_total_pts) {
    PCerr << "Error: nonzero number of sample points required in ProjectOrthog"
	  << "PolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }

  // Array sizing can be divided into two parts:
  // > data used in all cases (size in allocate_arrays())
  // > data not used in expansion import case (size here)
  allocate_arrays();
#ifdef DEBUG
  gradient_check();
#endif // DEBUG

  // calculate polynomial chaos coefficients
  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    // verify quad_order stencil matches num_total_pts
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    const UShortArray&   quad_order = tpq_driver->quadrature_order();
    if (quad_order.size() != num_v) {
      PCerr << "Error: quadrature order array is not consistent with number of "
	    << "variables (" << num_v << ")\n       in ProjectOrthogPoly"
	    << "Approximation::compute_coefficients()." << std::endl;
      abort_handler(-1);
    }
    size_t num_colloc_pts = 1;
    for (i=0; i<num_v; ++i)
      num_colloc_pts *= quad_order[i];
    if (num_total_pts != num_colloc_pts) {
      PCerr << "Error: number of current points (" << num_total_pts
	    << ") is not consistent with\n       quadrature data in Project"
	    << "OrthogPolyApproximation::compute_coefficients()." << std::endl;
      abort_handler(-1);
    }

    // single expansion integration
    integration_checks();
    integrate_expansion(data_rep->multiIndex, surrData.variables_data(),
			surrData.response_data(),
			tpq_driver->type1_weight_sets(),
			expansionCoeffs, expansionCoeffGrads);
    break;
  }
  case CUBATURE: {
    // single expansion integration
    integration_checks();
    CubatureDriver* cub_driver = data_rep->cub_driver();
    integrate_expansion(data_rep->multiIndex, surrData.variables_data(),
			surrData.response_data(),
			cub_driver->type1_weight_sets(),
			expansionCoeffs, expansionCoeffGrads);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    // multiple tensor expansion integrations
    if (expansionCoeffFlag)     expansionCoeffs = 0.;
    if (expansionCoeffGradFlag) expansionCoeffGrads = 0.;
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const IntArray& sm_coeffs = csg_driver->smolyak_coefficients();
    size_t i, num_tensor_grids = data_rep->tpMultiIndex.size(); int coeff;
    SDVArray tp_data_vars; SDRArray tp_data_resp;
    RealVector tp_wts, tp_coeffs; RealMatrix tp_coeff_grads;
    bool store_tp = (data_rep->expConfigOptions.refinementControl ==
		     DIMENSION_ADAPTIVE_CONTROL_GENERALIZED);
    // loop over tensor-products, forming sub-expansions, and sum them up
    // Note: SharedOrthogPolyApproxData::allocate_data() uses
    // sparse_grid_multi_index() to build multiIndex with append_multi_index()
    for (i=0; i<num_tensor_grids; ++i) {
      // form tp_data_vars, tp_data_resp, tp_wts using collocKey et al.
      integration_data(i, tp_data_vars, tp_data_resp, tp_wts);

      // form tp_multi_index from tpMultiIndexMap
      //map_tensor_product_multi_index(tp_multi_index, i);

      // form tp expansion coeffs
      RealVector& tp_coeffs_i = (store_tp) ? tpExpansionCoeffs[i] : tp_coeffs;
      RealMatrix& tp_grads_i
	= (store_tp) ? tpExpansionCoeffGrads[i] : tp_coeff_grads;
      integrate_expansion(data_rep->tpMultiIndex[i], tp_data_vars, tp_data_resp,
			  tp_wts, tp_coeffs_i, tp_grads_i);

      // sum tensor product coeffs/grads into expansion coeffs/grads
      coeff = sm_coeffs[i];
      if (coeff)
	overlay_expansion(data_rep->tpMultiIndexMap[i],
			  tp_coeffs_i, tp_grads_i, coeff);
    }
    break;
  }
  case SAMPLING:
    surrData.data_checks();
    expectation();
    break;
  default:
    PCerr << "Error: unsupported expCoeffsSolnApproach in ProjectOrthogPoly"
	  << "Approximation::compute_coefficients()" << std::endl;
    abort_handler(-1);
    break;
  }

  computedMean = computedVariance = 0;
}


void ProjectOrthogPolyApproximation::increment_coefficients()
{
  // tpMultiIndex{,Map,MapRef} already updated in
  // SharedProjectOrthogPolyApproxData::increment_data()
  size_t last_index = tpExpansionCoeffs.size();
  RealVector rv; tpExpansionCoeffs.push_back(rv);
  RealMatrix rm; tpExpansionCoeffGrads.push_back(rm);

  // resize component Sobol' array sizes to pick up new interaction terms
  allocate_component_sobol();

  // resize the PCE
  resize_expansion();

  // form tp_data_pts, tp_wts using collocKey et al.
  SDVArray tp_data_vars; SDRArray tp_data_resp; RealVector tp_wts;
  integration_data(last_index, tp_data_vars, tp_data_resp, tp_wts);
  // form trial expansion coeffs/grads
  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;
  integrate_expansion(data_rep->tpMultiIndex[last_index], tp_data_vars,
		      tp_data_resp, tp_wts, tpExpansionCoeffs[last_index],
		      tpExpansionCoeffGrads[last_index]);
  // sum trial expansion into expansionCoeffs/expansionCoeffGrads
  append_tensor_expansions(last_index);

  computedMean = computedVariance = 0;
}


void ProjectOrthogPolyApproximation::decrement_coefficients()
{
  // reset expansion{Coeffs,CoeffGrads}: (set in append_tensor_expansions())
  expansionCoeffs     = prevExpCoeffs;
  expansionCoeffGrads = prevExpCoeffGrads;

  // don't update Sobol' array sizes for decrement, push, or finalize

  // expansion resize not necessary since (1) already updated from prevExp
  // and (2) not updating expansion on decrement (next increment updates).
  //resize_expansion();

  // reset tensor-product bookkeeping and save restorable data
  poppedTPExpCoeffs.push_back(tpExpansionCoeffs.back());
  poppedTPExpCoeffGrads.push_back(tpExpansionCoeffGrads.back());
  tpExpansionCoeffs.pop_back();  tpExpansionCoeffGrads.pop_back();

  computedMean = computedVariance = 0;
}


void ProjectOrthogPolyApproximation::push_coefficients()
{
  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;

  // move previous expansion data to current expansion
  size_t last_index = tpExpansionCoeffs.size();
  size_t index_star = data_rep->pushIndex;

  std::deque<RealVector>::iterator cit = poppedTPExpCoeffs.begin();
  std::deque<RealMatrix>::iterator git = poppedTPExpCoeffGrads.begin();
  std::advance(cit, index_star); std::advance(git, index_star);

  tpExpansionCoeffs.push_back(*cit);     poppedTPExpCoeffs.erase(cit);
  tpExpansionCoeffGrads.push_back(*git); poppedTPExpCoeffGrads.erase(git);

  // don't update Sobol' array sizes for decrement, push, or finalize

  // resize the PCE
  resize_expansion();

  // sum trial expansion into expansionCoeffs/expansionCoeffGrads
  append_tensor_expansions(last_index);

  computedMean = computedVariance = 0;
}


void ProjectOrthogPolyApproximation::finalize_coefficients()
{
  size_t start_index = tpExpansionCoeffs.size();

  // don't update Sobol' array sizes for decrement, push, or finalize
  resize_expansion();
  // move previous expansion data to current expansion
  tpExpansionCoeffs.insert(tpExpansionCoeffs.end(), poppedTPExpCoeffs.begin(),
    poppedTPExpCoeffs.end());
  tpExpansionCoeffGrads.insert(tpExpansionCoeffGrads.end(),
    poppedTPExpCoeffGrads.begin(), poppedTPExpCoeffGrads.end());

  poppedTPExpCoeffs.clear();       poppedTPExpCoeffGrads.clear();
  // sum remaining trial expansions into expansionCoeffs/expansionCoeffGrads
  append_tensor_expansions(start_index);

  computedMean = computedVariance = 0;
}


void ProjectOrthogPolyApproximation::
append_tensor_expansions(size_t start_index)
{
  // for use in decrement_coefficients()
  prevExpCoeffs = expansionCoeffs; prevExpCoeffGrads = expansionCoeffGrads;

  // update expansion{Coeffs,CoeffGrads} using a hierarchical update
  // rather than building from scratch
  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;
  CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
  const IntArray&     sm_coeffs = csg_driver->smolyak_coefficients();
  const IntArray& sm_coeffs_ref = csg_driver->smolyak_coefficients_reference();
#ifdef DEBUG
  PCout << "In ProjectOrthogPolyApproximation::append_tensor_expansions() with "
	<< "start index " << start_index << "\nsm_coeffs:\n" << sm_coeffs
	<< "sm_coeffs_ref:\n" << sm_coeffs_ref << std::endl;
#endif // DEBUG

  // add trial expansions
  size_t index, num_tensor_grids = sm_coeffs.size();
  int coeff, delta_coeff;
  for (index=start_index; index<num_tensor_grids; ++index) {
    coeff = sm_coeffs[index];
    if (coeff)
      overlay_expansion(data_rep->tpMultiIndexMap[index],
			tpExpansionCoeffs[index],
			tpExpansionCoeffGrads[index], coeff);
#ifdef DEBUG
    PCout << "Trial set sm_coeff = " << coeff << "\ntpExpansionCoeffs:\n";
    write_data(PCout, tpExpansionCoeffs[index]);
    PCout << "\ntpMultiIndexMap:\n" << data_rep->tpMultiIndexMap[index] << '\n';
#endif // DEBUG
  }
  // update other expansion contributions with a changed smolyak coefficient
  for (index=0; index<start_index; ++index) {
    // add new, subtract previous
    delta_coeff = sm_coeffs[index] - sm_coeffs_ref[index];
#ifdef DEBUG
    PCout << "Old set delta_coeff = " << delta_coeff
	  << "\ntpExpansionCoeffs:\n";
    write_data(PCout, tpExpansionCoeffs[index]);
    PCout << "\ntpMultiIndexMap:\n" << data_rep->tpMultiIndexMap[index] << '\n';
#endif // DEBUG
    if (delta_coeff)
      overlay_expansion(data_rep->tpMultiIndexMap[index],
			tpExpansionCoeffs[index],
			tpExpansionCoeffGrads[index], delta_coeff);
  }
}


void ProjectOrthogPolyApproximation::
integration_data(size_t tp_index, SDVArray& tp_data_vars,
		 SDRArray& tp_data_resp, RealVector& tp_weights)
{
  // extract tensor vars/resp from surrData and tensor wts from type1CollocWts1D
  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;
  CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
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
  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;

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
  PCout << "expansion_coeffs:\n"; write_data(PCout, exp_coeffs);
  if (exp_coeff_grads.numRows()) {
    PCout << "expansion_coeff_grads:\n";
    write_data(PCout, exp_coeff_grads, true, true, true);
  }
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
  // "lhs" or "random", no weights needed
  size_t i, j, k, num_surr_data_pts = surrData.points(), num_failed_surr_fn = 0,
    num_failed_surr_grad = 0, num_deriv_vars = expansionCoeffGrads.numRows();
  SizetShortMap::const_iterator fit;
  short                failed_anchor_data = surrData.failed_anchor_data();
  const SizetShortMap& failed_resp_data   = surrData.failed_response_data();
  for (fit=failed_resp_data.begin(); fit!=failed_resp_data.end(); ++fit) {
    if (fit->second & 1) ++num_failed_surr_fn;
    if (fit->second & 2) ++num_failed_surr_grad;
  }
  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;
  size_t num_data_pts_fn = num_surr_data_pts - num_failed_surr_fn,
    num_data_pts_grad    = num_surr_data_pts - num_failed_surr_grad,
    num_total_pts_fn = num_data_pts_fn, num_total_pts_grad = num_data_pts_grad,
    num_exp_terms = data_rep->multiIndex.size();
  bool anchor_fn = false, anchor_grad = false;
  if (surrData.anchor()) {
    if (expansionCoeffFlag     && !(failed_anchor_data & 1))
      { anchor_fn   = true; ++num_total_pts_fn; }
    if (expansionCoeffGradFlag && !(failed_anchor_data & 2))
      { anchor_grad = true; ++num_total_pts_grad; }
  }
  if (expansionCoeffFlag)
    PCout << "Expectations of " << num_exp_terms << " chaos coefficients "
	  << "using " << num_total_pts_fn << " observations.\n";
  if (expansionCoeffGradFlag)
    PCout << "Expectations of gradients of " << num_exp_terms << " chaos "
	  << "coefficients using " << num_total_pts_grad << " observations.\n";

  /*
  // The following implementation evaluates all PCE coefficients
  // using a consistent expectation formulation
  for (i=0; i<num_exp_terms; ++i) {
    Real& exp_coeff_i = expansionCoeffs[i];
    const UShortArray& mi_i = data_rep->multiIndex[i];
    exp_coeff_i = (anchor_fn) ?
      surrData.anchor_function() * data_rep->multivariate_polynomial(
        surrData.anchor_continuous_variables(), mi_i) : 0.;
    for (j=0; j<num_data_pts; ++j)
      exp_coeff_i += surrData.response_function(j) * data_rep->
        multivariate_polynomial(surrData.continuous_variables(j), mi_i);
    exp_coeff_i /= num_total_pts * data_rep->norm_squared(mi_i);
#ifdef DEBUG
    PCout << "coeff[" << i << "] = " << exp_coeff_i
	  << " norm squared[" << i <<"] = " << data_rep->norm_squared(mi_i)
	  << '\n';
#endif // DEBUG
  }
  */

  // This alternate implementation evaluates the first PCE coefficient (the
  // response mean) as an expectation and then removes the mean from the
  // expectation evaluation of all subsequent coefficients.  This approach
  // has been observed to result in better results for small sample sizes.
  Real empty_r;
  Real& mean      = (expansionCoeffFlag)     ? expansionCoeffs[0] : empty_r;
  Real* mean_grad = (expansionCoeffGradFlag) ? expansionCoeffGrads[0] : NULL;
  if (expansionCoeffFlag) {
    if (anchor_fn)   mean = surrData.anchor_function();
    else             expansionCoeffs = 0.;
  }
  if (expansionCoeffGradFlag) {
    if (anchor_grad) copy_data(surrData.anchor_gradient().values(),
			       num_deriv_vars, mean_grad);
    else             expansionCoeffGrads = 0.;
  }
  for (k=0, fit=failed_resp_data.begin(); k<num_surr_data_pts; ++k) {
    bool add_val = expansionCoeffFlag, add_grad = expansionCoeffGradFlag;
    fail_booleans(fit, k, add_val, add_grad);
    if (add_val)
      mean += surrData.response_function(k);
    if (add_grad) {
      const RealVector& curr_pt_grad = surrData.response_gradient(k);
      for (j=0; j<num_deriv_vars; ++j)
	mean_grad[j] += curr_pt_grad[j];
    }
  }
  if (expansionCoeffFlag)
    mean /= num_total_pts_fn;
  if (expansionCoeffGradFlag)
    for (j=0; j<num_deriv_vars; ++j)
      mean_grad[j] /= num_total_pts_grad;

  Real chaos_sample, resp_fn_minus_mean, norm_sq; Real* exp_grad_i;
  RealVector resp_grad_minus_mean;
  if (expansionCoeffGradFlag)
    resp_grad_minus_mean.sizeUninitialized(num_deriv_vars);
  if (anchor_fn || anchor_grad) {
    if (anchor_fn)
      resp_fn_minus_mean = surrData.anchor_function() - mean;
    if (anchor_grad) {
      const RealVector& anch_grad = surrData.anchor_gradient();
      for (j=0; j<num_deriv_vars; ++j)
	resp_grad_minus_mean[j] = anch_grad[j] - mean_grad[j];
    }
    const RealVector& c_vars = surrData.anchor_continuous_variables();
    for (i=1; i<num_exp_terms; ++i) {
      chaos_sample
	= data_rep->multivariate_polynomial(c_vars, data_rep->multiIndex[i]);
      if (anchor_fn)
	expansionCoeffs[i] = resp_fn_minus_mean * chaos_sample;
      if (anchor_grad) {
	exp_grad_i = expansionCoeffGrads[i];
	for (j=0; j<num_deriv_vars; ++j)
	  exp_grad_i[j] = resp_grad_minus_mean[j] * chaos_sample;
      }
    }
  }
  for (k=0, fit=failed_resp_data.begin(); k<num_surr_data_pts; ++k) {
    bool add_val = expansionCoeffFlag, add_grad = expansionCoeffGradFlag;
    fail_booleans(fit, k, add_val, add_grad);
    if (add_val)
      resp_fn_minus_mean = surrData.response_function(k) - mean;
    if (add_grad) {
      const RealVector& resp_grad = surrData.response_gradient(k);
      for (j=0; j<num_deriv_vars; ++j)
	resp_grad_minus_mean[j] = resp_grad[j] - mean_grad[j];
    }
    const RealVector& c_vars = surrData.continuous_variables(k);
    for (i=1; i<num_exp_terms; ++i) {
      chaos_sample
	= data_rep->multivariate_polynomial(c_vars, data_rep->multiIndex[i]);
      if (add_val)
	expansionCoeffs[i] += resp_fn_minus_mean * chaos_sample;
      if (add_grad) {
	exp_grad_i = expansionCoeffGrads[i];
	for (j=0; j<num_deriv_vars; ++j)
	  exp_grad_i[j] += resp_grad_minus_mean[j] * chaos_sample;
      }
    }
  }
  for (i=1; i<num_exp_terms; ++i) {
    norm_sq = data_rep->norm_squared(data_rep->multiIndex[i]);
    if (expansionCoeffFlag)
      expansionCoeffs[i] /= norm_sq * num_total_pts_fn;
    if (expansionCoeffGradFlag) {
      exp_grad_i = expansionCoeffGrads[i];
      for (j=0; j<num_deriv_vars; ++j)
	exp_grad_i[j] /= norm_sq * num_total_pts_grad;
    }
#ifdef DEBUG
    PCout << "coeff[" << i << "] = " << expansionCoeffs[i]
        //<< "coeff_grad[" << i <<"] = " << exp_grad_i
	  << " norm squared[" << i <<"] = " << norm_sq << '\n';
#endif // DEBUG
  }
}


void ProjectOrthogPolyApproximation::
integrate_response_moments(size_t num_moments)
{
  size_t i, s, num_pts = surrData.points(), num_stored = storedExpCoeffs.size();
  bool anchor_pt = surrData.anchor();
  if (anchor_pt) ++num_pts;

  // define data_coeffs
  RealVector data_coeffs(num_pts);
  if (anchor_pt) {
    data_coeffs[0] = surrData.anchor_function();
    for (i=1; i<num_pts; ++i)
      data_coeffs[i] = surrData.response_function(i-1);
  }
  else
    for (i=0; i<num_pts; ++i)
      data_coeffs[i] = surrData.response_function(i);

  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;
  if (data_rep->storedExpCombineType && num_stored) {
    // update data_coeffs using evaluations from stored expansions
    switch (data_rep->storedExpCombineType) {
    case ADD_COMBINE:
      if (anchor_pt) {
	const RealVector& a_c_vars = surrData.anchor_continuous_variables();
	for (s=0; s<num_stored; ++s)
	  data_coeffs[0] += stored_value(a_c_vars, s);
	for (i=1; i<num_pts; ++i) {
	  const RealVector& c_vars = surrData.continuous_variables(i-1);
	  for (s=0; s<num_stored; ++s)
	    data_coeffs[i] += stored_value(c_vars, s);
	}
      }
      else
	for (i=0; i<num_pts; ++i) {
	  const RealVector& c_vars = surrData.continuous_variables(i);
	  for (s=0; s<num_stored; ++s)
	    data_coeffs[i] += stored_value(c_vars, s);
	}
      break;
    case MULT_COMBINE:
      if (anchor_pt) {
	const RealVector& a_c_vars = surrData.anchor_continuous_variables();
	for (s=0; s<num_stored; ++s)
	  data_coeffs[0] *= stored_value(a_c_vars, s);
	for (i=1; i<num_pts; ++i) {
	  const RealVector& c_vars = surrData.continuous_variables(i-1);
	  for (s=0; s<num_stored; ++s)
	    data_coeffs[i] *= stored_value(c_vars, s);
	}
      }
      else
	for (i=0; i<num_pts; ++i) {
	  const RealVector& c_vars = surrData.continuous_variables(i);
	  for (s=0; s<num_stored; ++s)
	    data_coeffs[i] *= stored_value(c_vars, s);
	}
      break;
    }
    // stored data may now be cleared
    if (expansionCoeffFlag)     storedExpCoeffs.clear();
    if (expansionCoeffGradFlag) storedExpCoeffGrads.clear();
  }

  // update numericalMoments based on data_coeffs
  if (numericalMoments.length() != num_moments)
    numericalMoments.sizeUninitialized(num_moments);
  integrate_moments(data_coeffs, data_rep->driverRep->type1_weight_sets(),
		    numericalMoments);
}


Real ProjectOrthogPolyApproximation::value(const RealVector& x)
{
  // sum expansion to get response value prediction

  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE:
    if (data_rep->storedExpCombineType) // not guaranteed to use tensor indexing
      return OrthogPolyApproximation::value(x);
    else { // Horner's rule approach applicable for tensor indexing
      // Error check for required data
      if (!expansionCoeffFlag) {
	PCerr << "Error: expansion coefficients not defined in "
	      << "ProjectOrthogPolyApproximation::value()" << std::endl;
	abort_handler(-1);
      }
      TensorProductDriver* tpq_driver = data_rep->tpq_driver();
      RealVector accumulator(sharedDataRep->numVars); // init to 0.
      return data_rep->tensor_product_value(x, expansionCoeffs,
	data_rep->approxOrder, data_rep->multiIndex, accumulator);
    }
    break;
  /*
  case COMBINED_SPARSE_GRID: {
    // Horner's rule approach requires storage of tpExpansionCoeffs in
    // compute_coefficients().  For now, leave store_tp as is and use
    // default approach if tpExpansionCoeffs is empty.  In addition,
    // tp arrays are not currently updated for expansion combinations.
    if (tpExpansionCoeffs.empty() || data_rep->storedExpCombineType)//most cases
      return OrthogPolyApproximation::value(x);
    else { // generalized sparse grid case
      // Error check for required data
      if (!expansionCoeffFlag) {
	PCerr << "Error: expansion coefficients not defined in "
	      << "ProjectOrthogPolyApproximation::value()" << std::endl;
	abort_handler(-1);
      }
      CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
      const UShort2DArray& sm_mi     = csg_driver->smolyak_multi_index();
      const IntArray&      sm_coeffs = csg_driver->smolyak_coefficients();
      RealVector accumulator(sharedDataRep->numVars); // init to 0.
      Real approx_val = 0.;
      size_t i, num_sm_mi = sm_mi.size(); int sm_coeff;
      for (i=0; i<num_sm_mi; ++i) {
	sm_coeff = sm_coeffs[i];
	if (sm_coeff)
	  approx_val += sm_coeff * data_rep->
	    tensor_product_value(x, tpExpansionCoeffs[i],
				 tpApproxOrders[i], // TO DO
				 tpMultiIndex[i], accumulator);
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
stored_value(const RealVector& x, size_t index)
{
  // sum expansion to get response value prediction

  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: { // Horner's rule approach
    // Error check for required data
    size_t i, num_stored_terms = data_rep->storedMultiIndex[index].size();
    if (!num_stored_terms ||
	storedExpCoeffs[index].length() != num_stored_terms) {
      PCerr << "Error: stored expansion coefficients not available in "
	    << "ProjectOrthogPolyApproximation::stored_value()" << std::endl;
      abort_handler(-1);
    }
    // Note: requires tensor indexing in storedMultiIndex (see OPA::value(x)),
    // which is safe to assume prior to support of >2 levels of fidelity.
    RealVector accumulator(sharedDataRep->numVars); // init to 0.
    return data_rep->
      tensor_product_value(x, storedExpCoeffs[index],
			   data_rep->storedApproxOrder[index],
			   data_rep->storedMultiIndex[index], accumulator);
    break;
  }
  // Horner's rule approach would require storage of tensor product components
  //case COMBINED_SPARSE_GRID:
    //break;
  default: // other cases are total-order expansions
    return OrthogPolyApproximation::stored_value(x, index);
    break;
  }
}

} // namespace Pecos
