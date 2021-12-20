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

  // update modSurrData from surrData
  synchronize_surrogate_data();

  // For testing of anchor point logic:
  //modSurrData.anchor_index(0); // treat 1st SDV,SDR as anchor point

  // anchor point, if present, is handled differently for different
  // expCoeffsSolnApproach settings:
  //   SAMPLING:   treat it as another data point
  //   QUADRATURE/CUBATURE/*_SPARSE_GRID: error
  //   LEAST_SQ_REGRESSION: use equality-constrained least squares
  if (!modSurrData.points()) {
    PCerr << "Error: nonzero number of sample points required in Polynomial"
	  << "Approximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }
}


void PolynomialApproximation::synchronize_surrogate_data()
{
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.discrepancyType) {
  case RECURSIVE_DISCREP:
    response_data_to_surplus_data();     break;
  case DISTINCT_DISCREP:
    response_data_to_discrepancy_data(); break;
  default: // allow use of approxData or modSurrData, even if no modifications
    // depending on ctor for modSurrData, may be NULL or just empty
    if (modSurrData.is_null()) // shared rep (linking once is sufficient)
      modSurrData = surrData;
    else if (modSurrData.points() == 0) // retain independence: shallow data cp
      modSurrData.copy_active(surrData, SHALLOW_COPY, SHALLOW_COPY);
    break;
  }
}


/** When using a recursive approximation, subtract current polynomial approx
    prediction from the modSurrData so that we form expansion on the surplus */
void PolynomialApproximation::response_data_to_surplus_data()
{
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;
  if (modSurrData.is_null()) modSurrData = SurrogateData(key);
  else                       modSurrData.active_key(key);
  if (key != surrData.active_key()) {
    PCerr << "Error: active key mismatch in PolynomialApproximation::"
	  << "response_data_to_surplus_data()." << std::endl;
    abort_handler(-1);
  }

  // Keep surrData and modSurrData distinct (don't share reps since copy() will
  // not restore independence when it is needed).  Rather use distinct arrays
  // with shallow copies of SDV,SDR instances when they will not be modified.
  // > level 0: both SDVs and SDRs can be shallow copies
  // > shallow copies of popped arrays support replicated pops on modSurrData
  //   (applied to distinct arrays of shared instances)
  const std::map<UShortArray, SDRArray>& resp_data_map
    = surrData.response_data_map();
  std::map<UShortArray, SDRArray>::const_iterator r_cit = resp_data_map.begin();
  if (key == r_cit->first) { // first entry -> no discrepancy/surplus
    modSurrData.copy_active(surrData, SHALLOW_COPY, SHALLOW_COPY);
    return;
  }

  // levels 1 -- L: modSurrData computes hierarchical surplus between level l
  // data (approxData[0] from Dakota::PecosApproximation) and the level l-1
  // approximation.
  // > SDR and popped SDR instances are distinct; only pop counts are copied
  modSurrData.copy_active_sdv(surrData, SHALLOW_COPY);
  modSurrData.size_active_sdr(surrData);
  modSurrData.anchor_index(surrData.anchor_index());
  modSurrData.pop_count_stack(surrData.pop_count_stack());
  // TO DO: do this more incrementally as data sets evolve across levels

  // More efficient to roll up contributions from each level expansion than
  // to combine expansions and then eval once.  Approaches are equivalent for
  // additive roll up.
  size_t i, num_pts = surrData.points();
  const SDVArray&   sdv_array =    surrData.variables_data();
  const SDRArray&   sdr_array =    surrData.response_data();
  SDRArray& surplus_sdr_array = modSurrData.response_data();
  Real delta_val; RealVector delta_grad;
  switch (data_rep->expConfigOptions.combineType) {
  case MULT_COMBINE: {
    Real orig_fn_val, stored_val, fn_val_j, fn_val_jm1;
    RealVector orig_fn_grad, fn_grad_j, fn_grad_jm1;
    size_t j, k, num_deriv_vars = modSurrData.num_derivative_variables();
    std::map<UShortArray, SDRArray>::const_iterator look_ahead_cit;
    for (i=0; i<num_pts; ++i) {
      const RealVector&          c_vars  = sdv_array[i].continuous_variables();
      const SurrogateDataResp& orig_sdr  = sdr_array[i];
      SurrogateDataResp&    surplus_sdr  = surplus_sdr_array[i];
      short                 surplus_bits = surplus_sdr.active_bits();
      delta_val = orig_fn_val = orig_sdr.response_function();
      if (surplus_bits & 2)
	copy_data(orig_sdr.response_gradient(), orig_fn_grad);
      for (r_cit=resp_data_map.begin(), j=0; r_cit->first!=key; ++r_cit, ++j) {
	stored_val = stored_value(c_vars, r_cit->first);
	delta_val /= stored_val;
	if (surplus_bits & 2) { // recurse using levels j and j-1
	  const RealVector& stored_grad
	    = stored_gradient_nonbasis_variables(c_vars, r_cit->first);
	  if (j == 0)
	    { fn_val_j = stored_val; fn_grad_j = stored_grad; }
	  else {
	    fn_val_j = fn_val_jm1 * stored_val;
	    for (k=0; k<num_deriv_vars; ++k)
	      fn_grad_j[k]  = ( fn_grad_jm1[k] * stored_val +
				fn_val_jm1 * stored_grad[j] );
	  }
	  look_ahead_cit = r_cit; ++look_ahead_cit;
	  if (look_ahead_cit->first == key)
	    for (k=0; k<num_deriv_vars; ++k)
	      delta_grad[k] = ( orig_fn_grad[k] - fn_grad_j[k] * delta_val )
	                    / fn_val_j;
	  else
	    { fn_val_jm1 = fn_val_j; fn_grad_jm1 = fn_grad_j; }
	}
      }
      if (surplus_bits & 1)
	surplus_sdr.response_function(delta_val);
      if (surplus_bits & 2)
	surplus_sdr.response_gradient(delta_grad);
    }
    break;
  }
  default: //case ADD_COMBINE: (correction specification not required)
    for (i=0; i<num_pts; ++i) {
      const RealVector&          c_vars  = sdv_array[i].continuous_variables();
      const SurrogateDataResp& orig_sdr  = sdr_array[i];
      SurrogateDataResp&    surplus_sdr  = surplus_sdr_array[i];
      short                 surplus_bits = surplus_sdr.active_bits();
      if (surplus_bits & 1) {
	delta_val = orig_sdr.response_function();
	for (r_cit = resp_data_map.begin(); r_cit->first != key; ++r_cit)
	  delta_val -= stored_value(c_vars, r_cit->first);
	surplus_sdr.response_function(delta_val);
      }
      if (surplus_bits & 2) {
	copy_data(orig_sdr.response_gradient(), delta_grad);
	for (r_cit = resp_data_map.begin(); r_cit->first != key; ++r_cit)
	  delta_grad -=
	    stored_gradient_nonbasis_variables(c_vars, r_cit->first);
	surplus_sdr.response_gradient(delta_grad);
      }
    }
    break;
  }
}


void PolynomialApproximation::response_data_to_discrepancy_data()
{
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;
  if (modSurrData.is_null()) modSurrData = SurrogateData(key);
  else                       modSurrData.active_key(key);
  if (key != surrData.active_key()) {
    PCerr << "Error: active key mismatch in PolynomialApproximation::"
	  << "response_data_to_discrepancy_data()." << std::endl;
    abort_handler(-1);
  }

  const std::map<UShortArray, SDRArray>& resp_data_map
    = surrData.response_data_map();
  std::map<UShortArray, SDRArray>::const_iterator r_cit = resp_data_map.begin();

  // Keep surrData and modSurrData distinct (don't share reps since copy() will
  // not restore independence when it is needed).  Rather use distinct arrays
  // with shallow copies of SDV,SDR instances when they will not be modified.
  // > level 0 mode is BYPASS_SURROGATE: both SDVs & SDRs can be shallow copies
  // > shallow copies of popped arrays support replicated pops on modSurrData
  //   (applied to distinct arrays of shared instances)
  if (key == r_cit->first) { // first entry -> no discrepancy/surplus
    modSurrData.copy_active(surrData, SHALLOW_COPY, SHALLOW_COPY);
    return;
  }

  // levels 1 -- L use AGGREGATED_MODELS mode: modSurrData computes discrepancy
  // between LF data (approxData[0] from Dakota::PecosApproximation; receives
  // level l-1 data) and HF data (approxData[1] from Dakota::PecosApproximation;
  // receives level l data).
  // > SDR and popped SDR instances are distinct; only pop counts are copied
  modSurrData.copy_active_sdv(surrData, SHALLOW_COPY);
  modSurrData.size_active_sdr(surrData);
  modSurrData.anchor_index(surrData.anchor_index());
  modSurrData.pop_count_stack(surrData.pop_count_stack());
  // TO DO: do this more incrementally as data sets evolve across levels

  // More efficient to roll up contributions from each level expansion than
  // to combine expansions and then eval once.  Approaches are equivalent for
  // additive roll up.
  size_t i, num_pts = surrData.points();
  const SDVArray&    sdv_array = surrData.variables_data();
  const SDRArray& hf_sdr_array = surrData.response_data();

  // ***************************************************************************
  // TO DO: improve encapsulation by passing surrogate key
  // (currently replicates logic in NonDExpansion::configure_indices())
  UShortArray lf_key;  paired_lf_key(key, lf_key);
  r_cit = resp_data_map.find(lf_key);
  const SDRArray& lf_sdr_array = r_cit->second;
  // ***************************************************************************

  DiscrepancyCalculator discrepCalc;
  SDRArray& delta_sdr_array = modSurrData.response_data();
  switch (data_rep->expConfigOptions.combineType) {
  case MULT_COMBINE: {
    size_t num_deriv_vars = surrData.num_derivative_variables();
    for (i=0; i<num_pts; ++i) {
      const SurrogateDataResp& lf_sdr  =    lf_sdr_array[i];
      const SurrogateDataResp& hf_sdr  =    hf_sdr_array[i];
      SurrogateDataResp&    delta_sdr  = delta_sdr_array[i];
      short                 delta_bits = delta_sdr.active_bits();
      short                 corr_order = (delta_bits & 2) ? 1 : 0;
      if (discrepCalc.check_multiplicative(hf_sdr.response_function(),
					   lf_sdr.response_function(),
					   corr_order)) {
	PCerr << "Error: numerical FPE in computing multiplicative discrepancy."
	      << "\n       Please change to additive discrepancy." << std::endl;
	abort_handler(-1);
      }
      if (delta_bits & 1)
	discrepCalc.compute_multiplicative(hf_sdr.response_function(),
	  lf_sdr.response_function(), delta_sdr.response_function_view());
      if (delta_bits & 2) {
	RealVector delta_grad(delta_sdr.response_gradient_view());	
	discrepCalc.compute_multiplicative(hf_sdr.response_function(),
	  hf_sdr.response_gradient(), lf_sdr.response_function(),
	  lf_sdr.response_gradient(), delta_grad);
      }
    }
    break;
  }
  default: //case ADD_COMBINE: (correction specification not required)
    for (i=0; i<num_pts; ++i) {
      const SurrogateDataResp& lf_sdr  =    lf_sdr_array[i];
      const SurrogateDataResp& hf_sdr  =    hf_sdr_array[i];
      SurrogateDataResp&    delta_sdr  = delta_sdr_array[i];
      short                 delta_bits = delta_sdr.active_bits();
      if (delta_bits & 1)
	discrepCalc.compute_additive(hf_sdr.response_function(),
	  lf_sdr.response_function(), delta_sdr.response_function_view());
      if (delta_bits & 2) {
	RealVector delta_grad(delta_sdr.response_gradient_view());
	discrepCalc.compute_additive(hf_sdr.response_gradient(),
	  lf_sdr.response_gradient(), delta_grad);
      }
    }
    break;
  }

  modSurrData.data_checks(); // from scratch (aggregate LF,HF failures)
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
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  size_t sobol_len = data_rep->sobolIndexMap.size();
  if (sobolIndices.length() != sobol_len)
    sobolIndices.sizeUninitialized(sobol_len);
}


void PolynomialApproximation::allocate_total_sobol()
{
  // number of total indices independent of number of component indices
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
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
