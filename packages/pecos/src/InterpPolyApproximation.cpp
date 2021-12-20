/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        InterpPolyApproximation
//- Description:  Implementation code for InterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "InterpPolyApproximation.hpp"
#include "SharedInterpPolyApproxData.hpp"

//#define DEBUG

namespace Pecos {


int InterpPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface of multiple dimensions
  return (expansionCoeffFlag || expansionCoeffGradFlag) ? 1 : 0;
}


void InterpPolyApproximation::allocate_arrays()
{
  std::shared_ptr<SharedInterpPolyApproxData> data_rep =
    std::static_pointer_cast<SharedInterpPolyApproxData>(sharedDataRep);
  update_active_iterators(data_rep->activeKey);

  allocate_total_sobol();
  allocate_component_sobol();

  size_t num_moments = (data_rep->nonRandomIndices.empty()) ? 4 : 2;
  RealVector& num_int_mom = primaryMomIter->second;
  if (num_int_mom.length() != num_moments)
    num_int_mom.sizeUninitialized(num_moments);
}


void InterpPolyApproximation::combined_to_active(bool clear_combined)
{
  // emulate new surrogate data following promotion of combined expansion
  // coeffs to active: simplifies stats computation for FINAL_RESULTS
  // (integration, VBD, etc.) for the combined-now-active coeffs
  synthetic_surrogate_data(surrData); // overwrite data for activeKey

  // migrate current moments
  PolynomialApproximation::combined_to_active(clear_combined);
}


void InterpPolyApproximation::test_interpolation()
{
  // SC should accurately interpolate the collocation data for TPQ and
  // SSG with fully nested rules, but will exhibit interpolation error
  // for SSG with other rules.
  if (expansionCoeffFlag) {
    std::shared_ptr<SharedPolyApproxData> data_rep =
      std::static_pointer_cast<SharedPolyApproxData>(sharedDataRep);
    bool use_derivs = data_rep->basisConfigOptions.useDerivs;

    const SDVArray& sdv_array = surrData.variables_data();
    const SDRArray& sdr_array = surrData.response_data();

    size_t i, index = 0, num_colloc_pts = surrData.points(),
      num_v = sharedDataRep->numVars, w7 = WRITE_PRECISION+7;
    Real interp_val, err, val_max_err = 0., grad_max_err = 0.,
      val_rmse = 0., grad_rmse = 0.;
    PCout << std::scientific << std::setprecision(WRITE_PRECISION);
    for (i=0; i<num_colloc_pts; ++i, ++index) {
      const RealVector& c_vars = sdv_array[index].continuous_variables();
      const SurrogateDataResp& sdr = sdr_array[index];
      Real resp_fn = sdr.response_function();
      interp_val = value(c_vars);
      err = (std::abs(resp_fn) > DBL_MIN) ? std::abs(1. - interp_val/resp_fn) :
	                                    std::abs(resp_fn - interp_val);
      PCout << "Colloc pt " << std::setw(3) << i+1
	    << ": truth value  = "  << std::setw(w7) << resp_fn
	    << " interpolant = "    << std::setw(w7) << interp_val
	    << " relative error = " << std::setw(w7) << err <<'\n';
      if (err > val_max_err) val_max_err = err;
      val_rmse += err * err;
      if (use_derivs) {
	const RealVector& resp_grad   = sdr.response_gradient();
	const RealVector& interp_grad = gradient_basis_variables(c_vars);
	for (size_t j=0; j<num_v; ++j) {
	  err = (std::abs(resp_grad[j]) > DBL_MIN) ?
	    std::abs(1. - interp_grad[j]/resp_grad[j]) :
	    std::abs(resp_grad[j] - interp_grad[j]);
	  PCout << "               " << "truth grad_" << j+1 << " = "
		<< std::setw(w7) << resp_grad[j]   << " interpolant = "
		<< std::setw(w7) << interp_grad[j] << " relative error = "
		<< std::setw(w7) << err << '\n';
	  if (err > grad_max_err) grad_max_err = err;
	  grad_rmse += err * err;
	}
      }
    }
    val_rmse = std::sqrt(val_rmse/num_colloc_pts);
    PCout << "\nValue interpolation errors:    " << std::setw(w7) << val_max_err
	  << " (max) " << std::setw(w7) << val_rmse << " (RMS)\n";
    if (use_derivs) {
      grad_rmse = std::sqrt(grad_rmse/num_colloc_pts/num_v);
      PCout << "Gradient interpolation errors: " << std::setw(w7)
	    << grad_max_err << " (max) " << std::setw(w7) << grad_rmse
	    << " (RMS)\n";
    }
  }
}


void InterpPolyApproximation::compute_component_sobol()
{
  // initialize partialVariance
  size_t sobol_len = sobolIndices.length();
  if (partialVariance.length() != sobol_len) partialVariance.size(sobol_len);
  else                                       partialVariance = 0.;

  // Compute the total expansion mean & variance.  For standard mode, these are
  // likely already available, as managed by computed{Mean,Variance} data reuse
  // trackers.  For all vars mode, they are computed without partial integration
  // restricted to the random indices.
  Real total_variance = variance();
  if (total_variance > SMALL_NUMBER) { // Solve for partial variances
    Real total_mean = mean();
    // 0th term gets subtracted as child in compute_partial_variance()
    partialVariance[0] = total_mean * total_mean;
    // compute the partial variances corresponding to Sobol' indices
    std::shared_ptr<SharedInterpPolyApproxData> data_rep =
      std::static_pointer_cast<SharedInterpPolyApproxData>(sharedDataRep);
    const BitArrayULongMap& index_map = data_rep->sobolIndexMap;
    for (BAULMCIter cit=index_map.begin(); cit!=index_map.end(); ++cit) {
      unsigned long index = cit->second;
      if (index) { // partialVariance[0] stores mean^2 offset
	compute_partial_variance(cit->first);
	sobolIndices[index] = partialVariance[index] / total_variance;
      }
    }
  }
  else // don't perform variance attribution for zero/negligible variance
    sobolIndices = 0.;
#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_component_sobol(), "
	<< "sobolIndices =\n" << sobolIndices;
#endif // DEBUG
}


void InterpPolyApproximation::compute_total_sobol()
{
  totalSobolIndices = 0.; // init total indices

  std::shared_ptr<SharedInterpPolyApproxData> data_rep =
    std::static_pointer_cast<SharedInterpPolyApproxData>(sharedDataRep);
  if (data_rep->expConfigOptions.vbdOrderLimit)
    // all component indices may not be available, so compute total indices
    // independently.  This approach parallels partial_variance_integral()
    // where the algorithm is separated by integration approach.
    compute_total_sobol_indices(); // virtual
  else {
    // all component indices are available, so add them up: totalSobolIndices
    // simply parses the bit sets of each of the sobolIndices and adds them to
    // each matching variable bin
    size_t k, num_v = sharedDataRep->numVars;
    const BitArrayULongMap& index_map = data_rep->sobolIndexMap;
    for (BAULMCIter cit=index_map.begin(); cit!=index_map.end(); ++cit)
      for (k=0; k<num_v; ++k)
        if (cit->first[k]) // var k is present in this Sobol' index
          totalSobolIndices[k] += sobolIndices[cit->second];
    // ensure non-negativity of indices
    //for (k=0; k<num_v; ++k)
    //  totalSobolIndices[k] = std::abs(totalSobolIndices[k]);
  }

#ifdef DEBUG
  PCout << "In InterpPolyApproximation::compute_total_sobol(), "
	<< "totalSobolIndices =\n" << totalSobolIndices;
#endif // DEBUG
}


/** Computes the variance of component functions.  Assumes that partial
    variances of all subsets of set_value have been computed in advance:
    compute_component_sobol() calls compute_partial_variance() using
    the ordered set_value's in sobolIndexMap. */
void InterpPolyApproximation::
compute_partial_variance(const BitArray& set_value)
{
  // derived classes override to compute partialVariance and then invoke
  // base version for post-processing of proper subsets

  // compute child subsets.  An alternate approach would be to iterate
  // over sobolIndexMap using it->is_proper_subset_of(set_value).
  BitArraySet children;
  proper_subsets(set_value, children);

  std::shared_ptr<SharedInterpPolyApproxData> data_rep =
    std::static_pointer_cast<SharedInterpPolyApproxData>(sharedDataRep);
  BitArrayULongMap& index_map = data_rep->sobolIndexMap;

  // index of parent set within sobolIndices and partialVariance
  unsigned long set_index;
  if (!children.empty())
    set_index = index_map[set_value];

  // subtract the contributions from child subsets.  partialVariance
  // calculations are computed by ordered traversal of sobolIndexMap, which
  // uses BitArray keys that sort lexicographically (in ascending order of
  // the corresponding unsigned integer) --> all proper subsets are available.
  for (BASIter it=children.begin(); it!=children.end(); ++it) {
    unsigned long subset_index  = index_map[*it];
    partialVariance[set_index] -= partialVariance[subset_index];
  }
}


/** For input parent set, recursively finds constituent child subsets
    with one fewer element */
void InterpPolyApproximation::
proper_subsets(const BitArray& parent_set, BitArraySet& children)
{
  size_t k, num_v = sharedDataRep->numVars;
  for (k=0; k<num_v; ++k)
    if (parent_set[k]) { // check for membership of variable k in parent set
      // remove var k from parent set to create child set
      BitArray child_set = parent_set; child_set.flip(k);
      // if child set has not been stored previously, insert it and recurse
      if (children.find(child_set) == children.end()) {
	children.insert(child_set);
	proper_subsets(child_set, children); // recurse until {0} is reached
      }
    }
}

} // namespace Pecos
