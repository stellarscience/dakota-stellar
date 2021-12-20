/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        PolynomialApproximation
//- Description:  Implementation code for PolynomialApproximation class
//-               
//- Owner:        Mike Eldred

#include "PolynomialApproximation.hpp"
#include "BasisPolynomial.hpp"
#include "SparseGridDriver.hpp"
#include "NumericGenOrthogPolynomial.hpp"
#include "DiscrepancyCalculator.hpp"
#include "ActiveKey.hpp"

//#define DEBUG


namespace Pecos {

void PolynomialApproximation::compute_coefficients()
{
  if (!expansionCoeffFlag && !expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in Polynomial"
	  << "Approximation::compute_coefficients().\n         Bypassing "
	  << "approximation construction." << std::endl;
    return;
  }

  // update surrData (if active aggregation key)
  synchronize_surrogate_data();

  // For testing of anchor point logic:
  //surrData.anchor_index(0); // treat 1st SDV,SDR as anchor point

  // anchor point, if present, is handled differently for different
  // expCoeffsSolnApproach settings:
  //   SAMPLING:   treat it as another data point
  //   QUADRATURE/CUBATURE/*_SPARSE_GRID: error
  //   LEAST_SQ_REGRESSION: use equality-constrained least squares
  if (!surrData.points()) {
    PCerr << "Error: nonzero number of sample points required in Polynomial"
	  << "Approximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }
}


void PolynomialApproximation::synchronize_surrogate_data()
{
  std::shared_ptr<SharedPolyApproxData> data_rep =
    std::static_pointer_cast<SharedPolyApproxData>(sharedDataRep);
  const ActiveKey& active_key = data_rep->activeKey;
  if (active_key != surrData.active_key()) {
    PCerr << "Error: active key mismatch in PolynomialApproximation::"
	  << "synchronize_surrogate_data()." << std::endl;
    abort_handler(-1);
  }

  // level 0: surrData non-aggregated key stores raw data
  short discrepancy  = data_rep->expConfigOptions.discrepReduction,
        combine_type = data_rep->expConfigOptions.combineType;
  if (!discrepancy || !active_key.aggregated() ||
      !active_key.raw_with_reduction_data())
    return;

  switch (discrepancy) {
  case RECURSIVE_DISCREPANCY:
    // When using a recursive discrepancy with additive/multiplicative corr,
    // we will subtract/divide the current polynomial approx prediction from
    // the new surrData so that we form an expansion on the surplus.  Prior
    // to using compute() to form the surplus, LF-hat must be generated and
    // will be stored within surrData in a format that compute() can utilize.
    generate_synthetic_data(surrData, active_key, combine_type);
    break;
  //case DISTINCT_DISCREPANCY:
    // When using a distinct discrepancy with additive/multiplicative corr,
    // we will subtract/divide the HF,LF pairs.  In this case, the data is
    // already provided within surrData and specific pairings are identified
    // by data groups.
  }
  // now compute the discrepancy between {HF,LF} or {HF,LF-hat} datasets
  DiscrepancyCalculator::compute(surrData, active_key, combine_type);
}


/** Compute the combined expansion prediction that corresponds to new surrData,
    prior to forming an expansion on the difference (surplus). */
void PolynomialApproximation::
generate_synthetic_data(SurrogateData& surr_data, const ActiveKey& active_key,
			short combine_type)
{
  // generate synthetic low-fidelity data at the high-fidelity points
  // and store within surr_data[lf_hat_key] where lf_hat_key is the trailing
  // portion of active_key.  This synthetic data then enables the computation
  // and emulation of a recursive discrepancy from hf - lf_hat differences
  // (surpluses) at the high-fidelity points
  ActiveKey hf_key, lf_hat_key; // LF-hat in surplus case
  active_key.extract_keys(hf_key, lf_hat_key);
  ActiveKey lf0_key = surr_data.filtered_key(SINGLETON_FILTER, 0); // *** Note: ActiveKey first sorts on group id

  // initialize surr_data[lf_hat_key]
  surr_data.active_key(lf_hat_key); // active key restored at fn end
  surr_data.variables_data(surr_data.variables_data(hf_key)); // shallow copies
  surr_data.anchor_index(surr_data.anchor_index(hf_key));
  surr_data.pop_count_stack(surr_data.pop_count_stack(hf_key));

  const SDRArray& hf_sdr_array = surr_data.response_data(hf_key);
  surr_data.size_active_sdr(hf_sdr_array); // size lf_hat_sdr_array
  const SDVArray&  sdv_array = surr_data.variables_data();
  SDRArray& lf_hat_sdr_array = surr_data.response_data();

  // extract all discrepancy data sets (which have expansions supporting
  // stored_{value,gradient} evaluations)
  const std::map<ActiveKey, SDRArray>& discrep_resp_map
    = surr_data.filtered_response_data_map(RAW_WITH_REDUCTION_DATA_FILTER);
  std::map<ActiveKey, SDRArray>::const_iterator cit;
  size_t i, num_pts = hf_sdr_array.size();
  switch (combine_type) {
  case MULT_COMBINE: {
    Real stored_val, fn_val_j, fn_val_jm1;
    RealVector fn_grad_j, fn_grad_jm1;
    size_t j, k, num_deriv_vars = surr_data.num_derivative_variables();
    for (i=0; i<num_pts; ++i) {
      const RealVector&       c_vars = sdv_array[i].continuous_variables();
      SurrogateDataResp& lf_hat_sdr  = lf_hat_sdr_array[i];
      short              lf_hat_bits = lf_hat_sdr.active_bits();
      // start from emulation of lowest fidelity QoI (LF-hat)
      fn_val_j = stored_value(c_vars, lf0_key); // coarsest fn
      if (lf_hat_bits & 2)                      // coarsest grad
	fn_grad_j = stored_gradient_nonbasis_variables(c_vars, lf0_key);
      // augment w/ emulation of discrepancies (Delta-hat) preceding active_key
      for (cit = discrep_resp_map.begin(), j=0;
	   cit->first != active_key; ++cit, ++j) {
	stored_val = stored_value(c_vars, cit->first); // Delta-hat
	if (lf_hat_bits & 2) { // recurse using levels j and j-1
	  const RealVector& stored_grad   // discrepancy gradient-hat
	    = stored_gradient_nonbasis_variables(c_vars, cit->first);
	  fn_val_jm1 = fn_val_j;  fn_grad_jm1 = fn_grad_j;
	  for (k=0; k<num_deriv_vars; ++k) // grad corrected to level j
	    fn_grad_j[k] = ( fn_grad_jm1[k] * stored_val +
			     fn_val_jm1 * stored_grad[k] );
	}
	fn_val_j *= stored_val; // fn corrected to level j
      }
      if (lf_hat_bits & 1)
	lf_hat_sdr.response_function(fn_val_j);
      if (lf_hat_bits & 2)
	lf_hat_sdr.response_gradient(fn_grad_j);
    }
    break;
  }
  default: { //case ADD_COMBINE: (correction specification not required)
    Real sum_val;  RealVector sum_grad;
    for (i=0; i<num_pts; ++i) {
      const RealVector&      c_vars  = sdv_array[i].continuous_variables();
      SurrogateDataResp& lf_hat_sdr  = lf_hat_sdr_array[i];
      short              lf_hat_bits = lf_hat_sdr.active_bits();
      if (lf_hat_bits & 1) {
	sum_val = stored_value(c_vars, lf0_key);
	for (cit = discrep_resp_map.begin(); cit->first != active_key; ++cit)
	  sum_val += stored_value(c_vars, cit->first);
	lf_hat_sdr.response_function(sum_val);
      }
      if (lf_hat_bits & 2) {
	sum_grad = stored_gradient_nonbasis_variables(c_vars, lf0_key);
	for (cit = discrep_resp_map.begin(); cit->first != active_key; ++cit)
	  sum_grad += stored_gradient_nonbasis_variables(c_vars, cit->first);
	lf_hat_sdr.response_gradient(sum_grad);
      }
    }
    break;
  }
  }

  surr_data.active_key(active_key); // restore
}


void PolynomialApproximation::combined_to_active(bool clear_combined)
{
  allocate_component_sobol(); // size sobolIndices from shared sobolIndexMap

  // migrate moments
  primaryMeanIter->second = combinedMeanBits;
  primaryVarIter->second  = combinedVarBits;
  std::shared_ptr<SharedPolyApproxData> data_rep =
    std::static_pointer_cast<SharedPolyApproxData>(sharedDataRep);
  if (!data_rep->nonRandomIndices.empty()) {
    const ActiveKey& key = data_rep->activeKey;
    xPrevMean[key] = xPrevCombMean;
    xPrevVar[key]  = xPrevCombVar;
  }
  if (clear_combined) {
    primaryMomIter->second.swap(combinedMoments);
    combinedMoments.resize(0);
    clear_combined_bits();
  }
  else
    primaryMomIter->second = combinedMoments; // deep copy
}


void PolynomialApproximation::
compute_moments(bool full_stats, bool combined_stats)
{
  // default for standard variables mode (2 moments for both full and
  // intermediate stats) is specialized by {Interp,ProjectOrthog}PolyApprox

  if (combined_stats) {
    if (combinedMoments.length() != 2) combinedMoments.resize(2);
    combined_mean(); combined_variance();
  }
  else {
    RealVector& mom1 = primaryMomIter->second;
    if (mom1.length() != 2) mom1.sizeUninitialized(2);
    mean(); variance();
    //standardize_moments(mom1);

    if (!full_stats && !secondaryMoments.empty()) secondaryMoments.resize(0);
  }
}


void PolynomialApproximation::
compute_moments(const RealVector& x, bool full_stats, bool combined_stats)
{
  // default for all variables mode (2 moments) is specialized by ...

  if (combined_stats) {
    if (combinedMoments.length() != 2) combinedMoments.resize(2);
    combined_mean(x); combined_variance(x);
  }
  else {
    RealVector& mom1 = primaryMomIter->second;
    if (mom1.length() != 2)            mom1.sizeUninitialized(2);
    mean(x);          variance(x);
    //standardize_moments(mom1);

    if (!full_stats && !secondaryMoments.empty()) secondaryMoments.resize(0);
  }
}


void PolynomialApproximation::
integrate_moments(const RealVector& coeffs, const RealVector& t1_wts,
		  RealVector& moments)
{
  // computes and stores the following moments:
  // > mean     (1st raw moment)
  // > variance (2nd central moment)
  // > skewness (3rd standardized moment)
  // > kurtosis (4th standardized moment with offset to eliminate "excess")

  // current support for this implementation: can't be open-ended since we
  // employ a specific combination of raw, central, and standardized moments
  size_t num_moments = moments.length();
  if (num_moments < 1 || num_moments > 4) {
    PCerr << "Error: unsupported number of moments requested in Polynomial"
	  << "Approximation::integrate_moments()" << std::endl;
    abort_handler(-1);
  }
  size_t i, j, num_pts = coeffs.length();
  if (t1_wts.length() != num_pts) {
    PCerr << "Error: mismatch in array lengths between integration driver "
	  << "weights (" << t1_wts.length() << ") and coefficients (" << num_pts
	  << ") in PolynomialApproximation::integrate_moments()." << std::endl;
    abort_handler(-1);
  }

#ifdef DEBUG
  PCout <<  "Coeffs in integrate_moments():\n" << coeffs
	<< "Weights in integrate_moments():\n" << t1_wts;
#endif

  // estimate 1st raw moment (mean)
  moments = 0.;
  Real& mean = moments[0];
  for (i=0; i<num_pts; ++i)
    mean += t1_wts[i] * coeffs[i];

  // estimate central moments 2 through num_moments
  if (num_moments > 1) {
    Real centered_fn, pow_fn;
    for (i=0; i<num_pts; ++i) {
      pow_fn = centered_fn = coeffs[i] - mean;
      for (j=1; j<num_moments; ++j) {
	pow_fn     *= centered_fn;
	moments[j] += t1_wts[i] * pow_fn;
      }
    }
  }

  // standardize third and higher central moments, if present
  //standardize_moments(moments);
}


void PolynomialApproximation::
integrate_moments(const RealVector& t1_coeffs, const RealMatrix& t2_coeffs,
		  const RealVector& t1_wts, const RealMatrix& t2_wts,
		  RealVector& moments)
{
  // computes and stores the following moments:
  // > mean     (1st raw moment)
  // > variance (2nd central moment)
  // > skewness (3rd standardized moment)
  // > kurtosis (4th standardized moment with offset to eliminate "excess")

  // current support for this implementation: can't be open-ended since we
  // employ a specific combination of raw, central, and standardized moments
  size_t num_moments = moments.length();
  if (num_moments < 1 || num_moments > 4) {
    PCerr << "Error: unsupported number of moments requested in Polynomial"
	  << "Approximation::integrate_moments()" << std::endl;
    abort_handler(-1);
  }
  size_t i, j, k, num_pts = t1_coeffs.length(), num_v = sharedDataRep->numVars;
  if (t1_wts.length() != num_pts || t2_wts.numCols() != num_pts ||
      t2_coeffs.numCols() != num_pts) {
    PCerr << "Error: mismatch in array lengths among integration driver "
	  << "weights ("  << t1_wts.length() << ", " << t2_wts.numCols()
	  << ") and coefficients (" << num_pts << ", " << t2_coeffs.numCols()
	  << ") in PolynomialApproximation::integrate_moments()." << std::endl;
    abort_handler(-1);
  }

  // estimate 1st raw moment (mean)
  moments = 0.;
  Real& mean = moments[0];
  for (i=0; i<num_pts; ++i) {
    mean += t1_wts[i] * t1_coeffs[i];
    const Real* coeff2_i = t2_coeffs[i];
    const Real*  t2_wt_i = t2_wts[i];
    for (j=0; j<num_v; ++j)
      mean += coeff2_i[j] * t2_wt_i[j];
  }

  // estimate central moments 2 through num_moments
  if (num_moments > 1) {
    Real centered_fn, pow_fn;
    for (i=0; i<num_pts; ++i) {
      pow_fn = centered_fn = t1_coeffs[i] - mean;
      const Real* coeff2_i = t2_coeffs[i];
      const Real*  t2_wt_i = t2_wts[i];
      for (j=1; j<num_moments; ++j) {
	Real& moment_j = moments[j];
	// type2 interpolation of (R - \mu)^n
	// --> interpolated gradients are n(R - \mu)^{n-1} dR/dx
	for (k=0; k<num_v; ++k)
	  moment_j += (j+1) * pow_fn * coeff2_i[k] * t2_wt_i[k];
	// type 1 interpolation of (R - \mu)^n
	pow_fn   *= centered_fn;
	moment_j += t1_wts[i] * pow_fn;
      }
    }
  }

  // convert central moments to std deviation/skewness/kurtosis
  //standardize_moments(moments);
}


void PolynomialApproximation::
standardize_moments(const RealVector& central_moments, RealVector& std_moments)
{
  size_t num_moments = central_moments.length();
  std_moments.sizeUninitialized(num_moments);
  if (num_moments >= 1) std_moments[0] = central_moments[0]; // mean
  if (num_moments <  2) return;

  const Real& var = central_moments[1];
  Real&   std_dev = std_moments[1];
  if (var > 0.) {
    // standardized moment k is E[((X-mu)/sigma)^k] = E[(X-mu)^k]/sigma^k
    std_dev = std::sqrt(var); // not standardized (2nd standardized moment is 1)
    Real pow_fn = var;
    for (size_t i=2; i<num_moments; ++i)
      { pow_fn *= std_dev; std_moments[i] = central_moments[i] / pow_fn; }
    // offset the fourth standardized moment to eliminate excess kurtosis
    if (num_moments > 3)
      std_moments[3] -= 3.;
  }
  else {
    // don't leave uninitialized, even if undefined
    for (size_t i=1; i<num_moments; ++i)
      std_moments[i] = 0.;
    // special case of zero variance is OK for num_moments == 2, but not higher
    if ( !(num_moments == 2 && var == 0.) ) // std_dev OK for var == 0.
      PCerr << "Warning: moments cannot be standardized due to non-positive "
	    << "variance.\n         Skipping standardization." << std::endl;
  }
}


void PolynomialApproximation::allocate_component_sobol()
{
  std::shared_ptr<SharedPolyApproxData> data_rep =
    std::static_pointer_cast<SharedPolyApproxData>(sharedDataRep);
  size_t sobol_len = data_rep->sobolIndexMap.size();
  if (sobolIndices.length() != sobol_len)
    sobolIndices.sizeUninitialized(sobol_len);
}


void PolynomialApproximation::allocate_total_sobol()
{
  // number of total indices independent of number of component indices
  std::shared_ptr<SharedPolyApproxData> data_rep =
    std::static_pointer_cast<SharedPolyApproxData>(sharedDataRep);
  if (totalSobolIndices.empty() && expansionCoeffFlag &&
      data_rep->expConfigOptions.vbdFlag)
    totalSobolIndices.sizeUninitialized(sharedDataRep->numVars);
}


void PolynomialApproximation::
initialize_covariance(PolynomialApproximation* poly_approx_2)
{ } // default is no-op


void PolynomialApproximation::clear_covariance_pointers()
{ } // default is no-op


void PolynomialApproximation::initialize_products()
{ } // default is no-op


bool PolynomialApproximation::product_interpolants()
{ return false; } // default


ULongULongMap PolynomialApproximation::sparse_sobol_index_map() const
{ return ULongULongMap(); } // default is empty map


size_t PolynomialApproximation::expansion_terms() const
{
  PCerr << "Error: expansion_terms() not defined for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return _NPOS;
}


Real PolynomialApproximation::combined_mean()
{
  PCerr << "Error: combined_mean() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::combined_mean(const RealVector& x)
{
  PCerr << "Error: combined_mean() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
combined_covariance(PolynomialApproximation* poly_approx_2)
{
  PCerr << "Error: combined_covariance() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
combined_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  PCerr << "Error: combined_covariance() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::beta(bool cdf_flag, Real z_bar)
{
  PCerr << "Error: beta() not available for this polynomial approximation type."
	<< std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
beta(const RealVector& x, bool cdf_flag, Real z_bar)
{
  PCerr << "Error: beta(x) not available for this polynomial approximation "
	<< "type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::combined_beta(bool cdf_flag, Real z_bar)
{
  PCerr << "Error: combined_beta() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
combined_beta(const RealVector& x, bool cdf_flag, Real z_bar)
{
  PCerr << "Error: combined_beta(x) not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_mean()
{
  PCerr << "Error: delta_mean() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_mean(const RealVector& x)
{
  PCerr << "Error: delta_mean(x) not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_combined_mean()
{
  PCerr << "Error: delta_combined_mean() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_combined_mean(const RealVector& x)
{
  PCerr << "Error: delta_combined_mean(x) not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_std_deviation()
{
  PCerr << "Error: delta_std_deviation() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_std_deviation(const RealVector& x)
{
  PCerr << "Error: delta_std_deviation(x) not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_combined_std_deviation()
{
  PCerr << "Error: delta_combined_std_deviation() not available for this "
	<< "polynomial approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_combined_std_deviation(const RealVector& x)
{
  PCerr << "Error: delta_combined_std_deviation(x) not available for this "
	<< "polynomial approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_covariance(PolynomialApproximation* poly_approx_2)
{
  PCerr << "Error: delta_covariance() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  PCerr << "Error: delta_covariance() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_combined_covariance(PolynomialApproximation* poly_approx_2)
{
  PCerr << "Error: delta_combined_covariance() not available for this "
	<< "polynomial approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_combined_covariance(const RealVector& x,
			  PolynomialApproximation* poly_approx_2)
{
  PCerr << "Error: delta_combined_covariance() not available for this "
	<< "polynomial approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_beta(bool cdf_flag, Real z_bar)
{
  PCerr << "Error: delta_beta() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_beta(const RealVector& x, bool cdf_flag, Real z_bar)
{
  PCerr << "Error: delta_beta(x) not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_z(bool cdf_flag, Real beta_bar)
{
  PCerr << "Error: delta_z() not available for this polynomial approximation "
	<< "type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_z(const RealVector& x, bool cdf_flag, Real beta_bar)
{
  PCerr << "Error: delta_z(x) not available for this polynomial approximation "
	<< "type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_combined_beta(bool cdf_flag, Real z_bar)
{
  PCerr << "Error: delta_combined_beta() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_combined_beta(const RealVector& x, bool cdf_flag, Real z_bar)
{
  PCerr << "Error: delta_combined_beta(x) not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::delta_combined_z(bool cdf_flag, Real beta_bar)
{
  PCerr << "Error: delta_combined_z() not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}


Real PolynomialApproximation::
delta_combined_z(const RealVector& x, bool cdf_flag, Real beta_bar)
{
  PCerr << "Error: delta_combined_z(x) not available for this polynomial "
	<< "approximation type." << std::endl;
  abort_handler(-1);
  return 0.;
}

} // namespace Pecos
