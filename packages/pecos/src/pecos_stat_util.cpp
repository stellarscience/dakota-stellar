/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "pecos_stat_util.hpp"
#include "RandomVariable.hpp"
#include "SurrogateData.hpp"

static const char rcsId[]="@(#) $Id: pecos_stat_util.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


void int_range_to_xy_pdf(int l_bnd,        int u_bnd,
			 RealArray& x_val, RealArray& y_val)
{
  int i, num_params = u_bnd - l_bnd + 1;
  x_val.resize(num_params);  y_val.assign(num_params, 1.);
  for (i=0; i<num_params; ++i)
    x_val[i] = (Real)(l_bnd + i);
}


void bins_to_xy_cdf(const RealRealMap& h_bin_prs,
		    RealArray& x_val, RealArray& y_val)
{
  size_t i, num_params = h_bin_prs.size(), last_index = num_params - 1;
  RRMCIter cit = h_bin_prs.begin();

  // Note: LHS continuous linear accumulates CDF with first y=0 and last y=1
  x_val.resize(num_params);  y_val.resize(num_params);
  for (i=0; i<num_params; ++i, ++cit)
    x_val[i] = cit->first;
  y_val[0] = 0.;  cit = h_bin_prs.begin();
  for (i=0; i<last_index; ++i, ++cit)
    y_val[i+1] = y_val[i] + cit->second * (x_val[i+1] - x_val[i]);
  // normalize if necessary (h_bin_pairs should have normalized counts)
  Real& cdf_last = y_val[last_index];
  if (cdf_last != 1.) {
    for (i=1; i<last_index; ++i)
      y_val[i] /= cdf_last;
    cdf_last = 1.;
  }
#ifdef DEBUG
  for (i=0; i<num_params; ++i)
    PCout << "hbuv cdf: x_val[" << i << "] is " << x_val[i]
	  << " y_val[" << i << "] is " << y_val[i] << '\n';
#endif // DEBUG
}


void intervals_to_xy_cdf(const RealRealPairRealMap& ci_bpa,
			 RealArray& x_val, RealArray& y_val)
{
  RealArray prob_dens;
  intervals_to_xy_pdf(ci_bpa, x_val, prob_dens);
  size_t i, num_params = x_val.size(), last_index = num_params - 1;

  // LHS continuous linear accumulates CDF with first y=0 and last y=1:
  // put the densities in a cumulative format necessary for LHS histograms.
  y_val.resize(num_params);
  y_val[0] = 0.;
  Real pdf;
  for (i=1; i<num_params; ++i) {
    pdf = prob_dens[i-1];
    y_val[i] = (pdf > 0.) ? y_val[i-1] + pdf * (x_val[i] - x_val[i-1]) :
                            y_val[i-1] + 1.e-4; // handle case of a gap
  }
  // normalize if necessary
  Real& cdf_last = y_val[last_index];
  if (cdf_last != 1.) {
    for (i=1; i<last_index; ++i)
      y_val[i] /= cdf_last;
    cdf_last = 1.;
  }
#ifdef DEBUG
  for (i=0; i<num_params; ++i)
    PCout << "ciuv cdf: x_val[" << i << "] is " << x_val[i]
	  << " y_val[" << i << "] is " << y_val[i] << '\n';
#endif // DEBUG
}

/// SAMPLE MOMENT ESTIMATORS FROM SAMPLES

void accumulate_mean(const RealVectorArray& fn_samples, size_t q,
		     size_t& num_samp, Real& mean)
{
  Real sum = 0., sample;
  size_t s, num_obs = fn_samples.size();
  for (s=0, num_samp=0; s<num_obs; ++s) {
    sample = fn_samples[s][q];
    if (std::isfinite(sample)) { // neither NaN nor +/-Inf
      sum += sample;
      ++num_samp;
    }
  }

  mean = (num_samp) ? sum / (Real)num_samp : 0.;
}


void accumulate_mean(const SDRArray& sdr_samples, size_t& num_samp, Real& mean)
{
  Real sum = 0., sample;
  size_t s, num_obs = sdr_samples.size();
  for (s=0, num_samp=0; s<num_obs; ++s) {
    sample = sdr_samples[s].response_function();
    if (std::isfinite(sample)) { // neither NaN nor +/-Inf
      sum += sample;
      ++num_samp;
    }
  }

  mean = (num_samp) ? sum / (Real)num_samp : 0.;
}


void accumulate_variance(const RealVectorArray& fn_samples, Real mean,
			 size_t q, size_t& num_samp, Real& var)
{
  // accumulate central moments (e.g., variance)
  size_t s, num_obs = fn_samples.size();
  Real sample, centered_fn, sum2 = 0.;
  for (s=0, num_samp=0; s<num_obs; ++s) {
    sample = fn_samples[s][q];
    if (std::isfinite(sample)) { // neither NaN nor +/-Inf
      centered_fn = sample - mean;
      sum2 += centered_fn * centered_fn;
      ++num_samp;
    }
  }

  // unbiased central moment estimator
  var = (num_samp > 1) ? sum2 / ((Real)num_samp - 1.) : 0.;
}


void accumulate_variance(const SDRArray& sdr_samples, Real mean,
			 size_t& num_samp, Real& var)
{
  // accumulate central moments (e.g., variance)
  size_t s, num_obs = sdr_samples.size();
  Real sample, centered_fn, sum2 = 0.;
  for (s=0, num_samp=0; s<num_obs; ++s) {
    sample = sdr_samples[s].response_function();
    if (std::isfinite(sample)) { // neither NaN nor +/-Inf
      centered_fn = sample - mean;
      sum2 += centered_fn * centered_fn;
      ++num_samp;
    }
  }

  // unbiased central moment estimator
  var = (num_samp > 1) ? sum2 / ((Real)num_samp - 1.) : 0.;
}


void accumulate_moments(const RealVectorArray& fn_samples, size_t q,
			short moments_type, Real* moments)
{
  // accumulate central moments (e.g., variance)
  size_t s, num_obs = fn_samples.size(), num_samp = 0;
  Real& mean = moments[0]; // already computed in accumulate_mean()
  Real sample, centered_fn, pow_fn, sum2 = 0., sum3 = 0., sum4 = 0.;
  for (s=0; s<num_obs; ++s) {
    sample = fn_samples[s][q];
    if (std::isfinite(sample)) { // neither NaN nor +/-Inf
      pow_fn  = centered_fn = sample - mean;
      pow_fn *= centered_fn; sum2 += pow_fn; // variance
      pow_fn *= centered_fn; sum3 += pow_fn; // 3rd central moment
      pow_fn *= centered_fn; sum4 += pow_fn; // 4th central moment
      ++num_samp;
    }
  }
  Real ns = (Real)num_samp, nm1 = ns - 1., nm2 = ns - 2., ns_sq = ns * ns;
  // biased central moment estimators (bypass and use raw sums below):
  //biased_cm2 = sum2 / ns; biased_cm3 = sum3 / ns; biased_cm4 = sum4 / ns;

  // unbiased moment estimators (central and standardized):
  bool central = (moments_type == CENTRAL_MOMENTS), pos_var = (sum2 > 0.);
  Real cm2 = sum2 / nm1;
  if (num_samp > 1 && pos_var)
    moments[1] = (central) ? cm2 : std::sqrt(cm2); // unbiased central : std dev
  else moments[1] = 0.;

  if (num_samp > 2 && pos_var)
    moments[2] = (central) ?
      // unbiased central:
      sum3 * ns / (nm1 * nm2) :
      // unbiased standard: cm3 / sigma^3 = N / (nm1 nm2) sum3 / (cm2)^1.5
      sum3 * ns / (nm1 * nm2 * std::pow(cm2, 1.5));
  else moments[2] = 0.;

  if (num_samp > 3 && pos_var)
    moments[3] = (central) ?
      // unbiased central (non-excess) from "Modeling with Data," Klemens 2009
      // (Appendix M):  unbiased_cm4 =
      // ( N^3 biased_cm4 / (N-1) - (6N - 9) unbiased_cm2^2 ) / (N^2 - 3N + 3)
      //(ns * ns * sum4 / nm1 - (6.*ns - 9.) * cm2 * cm2) / (ns*(ns - 3.) + 3.) :
      //[fm] above estimator is not unbiased since cm2 * cm2 is biased, unbiased correction:
      (ns_sq * sum4 / nm1 - (6.*ns - 9.)*(ns_sq - ns)/(ns_sq - 2. * ns + 3) * cm2 * cm2) / ( (ns*(ns - 3.) + 3.) - ( (6.*ns - 9.)*(ns_sq - ns) )/( ns * (ns_sq - 2.*ns + 3) ) ) :
      // unbiased standard (excess kurtosis) from Wikipedia ("Estimators of
      // population kurtosis")
      nm1 * ((ns + 1.) * ns * sum4 / (sum2*sum2) - 3.*nm1) / (nm2*(ns - 3.));
  else moments[3] = (central) ? 0. : -3.;
}


void accumulate_moment_gradients(const RealVectorArray& fn_samples,
				 const RealMatrixArray& grad_samples, size_t q,
				 short moments_type, Real mean, Real mom2,
				 Real* mean_grad, Real* mom2_grad)
{
  size_t s, v, num_obs = std::min(fn_samples.size(), grad_samples.size()),
    num_deriv_vars = (num_obs) ? grad_samples[0].numRows() : 0; // grads = V x Q

  for (v=0; v<num_deriv_vars; ++v)
    mean_grad[v] = mom2_grad[v] = 0.;

  SizetArray num_samp(num_deriv_vars, 0);
  for (s=0; s<num_obs; ++s) {
    // manage faults hierarchically as in Pecos::SurrogateData::response_check()
    Real fn = fn_samples[s][q];
    if (std::isfinite(fn)) {          // neither NaN nor +/-Inf
      const Real* grad = grad_samples[s][q];
      for (v=0; v<num_deriv_vars; ++v)
	if (std::isfinite(grad[v])) { // neither NaN nor +/-Inf
	  mean_grad[v] += grad[v];
	  mom2_grad[v] += fn * grad[v];
	  ++num_samp[v];
	}
    }
  }

  Real ns, nm1; size_t nsv;
  bool central_mom = (moments_type == CENTRAL_MOMENTS);
  for (v=0; v<num_deriv_vars; ++v) {
    nsv = num_samp[v];
    if (nsv) {
      ns = (Real)nsv;
      // dMean/ds = E[dQ/ds] --> unbiased estimator 1/N Sum(dQ/ds)
      mean_grad[v] /= ns;
    }
    if (nsv > 1) {
      nm1 = ns - 1.;
      // dVar/ds = 2 E[(Q - Mean)(dQ/ds - dMean/ds)] --> unbiased var estimator:
      // = 2/(N-1) [ Sum(Q dQ/ds) - Mean Sum(dQ/ds) -
      //             dMean/ds Sum(Q) + N Mean Mean/ds ]
      // = 2/(N-1) [ Sum(Q dQ/ds) - Mean dMean/ds (N + N - N) ]
      // = 2/(N-1) [ Sum(Q dQ/ds) - N Mean dMean/ds ]
      // dVar/ds = 2 Sigma dSigma/ds -->
      // dSigma/ds = [ Sum(Q dQ/ds) - N Mean dMean/ds ] / (Sigma (N-1))
      mom2_grad[v]  = (central_mom) ?
	2. * ( mom2_grad[v] - ns * mean * mean_grad[v] ) / nm1 :
	     ( mom2_grad[v] - ns * mean * mean_grad[v] ) / (mom2 * nm1);
    }
  }
}


/*
Note: these are not active since RandomVariable type is insufficient to
      define the source variable type in the presence of probability
      transformations (e.g., from design/epistemic/state to std uniform).

void design_state_subset(const std::vector<RandomVariable>& random_vars,
			 BitArray& subset, size_t start_set, size_t num_set)
{
  size_t i, num_rv = random_vars.size(), end_set = start_set + num_set;
  subset.resize(num_rv, false); // init bits to false
  for (i=start_set; i<end_set; ++i)
    // activate design + state vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_RANGE:    case DISCRETE_RANGE: case DISCRETE_SET_INT:
    case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
      subset.set(i); break;
    }
}


void uncertain_subset(const std::vector<RandomVariable>& random_vars,
		      BitArray& subset)
{
  size_t i, num_rv = random_vars.size();
  subset.resize(num_rv, true); // init bits to true
  for (i=0; i<num_rv; ++i)
    // deactivate complement of uncertain vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_RANGE:    case DISCRETE_RANGE: case DISCRETE_SET_INT:
    case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
      subset.reset(i); break;
    }
}


void aleatory_uncertain_subset(const std::vector<RandomVariable>& random_vars,
			       BitArray& subset)
{
  size_t i, num_rv = random_vars.size();
  subset.resize(num_rv, true); // init bits to true
  for (i=0; i<num_rv; ++i)
    // deactivate complement of aleatory uncertain vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_RANGE:    case DISCRETE_RANGE: case DISCRETE_SET_INT:
    case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
    case CONTINUOUS_INTERVAL_UNCERTAIN: case DISCRETE_INTERVAL_UNCERTAIN:
    case DISCRETE_UNCERTAIN_SET_INT:    case DISCRETE_UNCERTAIN_SET_STRING:
    case DISCRETE_UNCERTAIN_SET_REAL:
      subset.reset(i); break;
    }
}


void epistemic_uncertain_subset(const std::vector<RandomVariable>& random_vars,
				BitArray& subset)
{
  size_t i, num_rv = random_vars.size();
  subset.resize(num_rv, false); // init bits to false
  for (i=0; i<num_rv; ++i)
    // activate epistemic uncertain vars
    switch (random_vars[i].type()) {
    case CONTINUOUS_INTERVAL_UNCERTAIN: case DISCRETE_INTERVAL_UNCERTAIN:
    case DISCRETE_UNCERTAIN_SET_INT:    case DISCRETE_UNCERTAIN_SET_STRING:
    case DISCRETE_UNCERTAIN_SET_REAL:
      subset.set(i); break;
    }
}
*/

} // namespace Pecos
