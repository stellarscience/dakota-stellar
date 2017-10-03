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
//#define INTERPOLATION_TEST

namespace Pecos {


int InterpPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface of multiple dimensions
  return (expansionCoeffFlag || expansionCoeffGradFlag) ? 1 : 0;
}


void InterpPolyApproximation::allocate_arrays()
{
  allocate_total_sobol();
  allocate_component_sobol();

  if (numericalMoments.empty()) {
    SharedInterpPolyApproxData* data_rep
      = (SharedInterpPolyApproxData*)sharedDataRep;
    size_t num_moments = (data_rep->nonRandomIndices.empty()) ? 4 : 2;
    numericalMoments.sizeUninitialized(num_moments);
  }
}


void InterpPolyApproximation::compute_coefficients()
{
  if (!expansionCoeffFlag && !expansionCoeffGradFlag) {
    PCerr << "Warning: neither expansion coefficients nor expansion "
	  << "coefficient gradients\n         are active in "
	  << "InterpPolyApproximation::compute_coefficients().\n         "
	  << "Bypassing approximation construction." << std::endl;
    return;
  }

  // For testing of anchor point logic:
  //size_t index = surrData.points() - 1;
  //surrData.anchor_point(surrData.variables_data()[index],
  //                      surrData.response_data()[index]);
  //surrData.pop(1);

  size_t num_colloc_pts = surrData.points();
  if (surrData.anchor()) // anchor point, if present, is first expansionSample
    ++num_colloc_pts;
  if (!num_colloc_pts) {
    PCerr << "Error: nonzero number of sample points required in "
	  << "InterpPolyApproximation::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }

  allocate_arrays();
  compute_expansion_coefficients();

#ifdef INTERPOLATION_TEST
  // SC should accurately interpolate the collocation data for TPQ and
  // SSG with fully nested rules, but will exhibit interpolation error
  // for SSG with other rules.
  if (expansionCoeffFlag) {
    size_t i, index = 0, offset = (surrData.anchor()) ? 1 : 0,
      w7 = WRITE_PRECISION+7, num_v = sharedDataRep->numVars;
    Real interp_val, err, val_max_err = 0., grad_max_err = 0.,
      val_rmse = 0., grad_rmse = 0.;
    PCout << std::scientific << std::setprecision(WRITE_PRECISION);
    for (i=offset; i<num_colloc_pts; ++i, ++index) {
      const RealVector& c_vars = surrData.continuous_variables(index);
      Real      resp_fn = surrData.response_function(index);
      interp_val = value(c_vars);
      err = (std::abs(resp_fn) > DBL_MIN) ? std::abs(1. - interp_val/resp_fn) :
	                                    std::abs(resp_fn - interp_val);
      PCout << "Colloc pt " << std::setw(3) << i+1
	    << ": truth value  = "  << std::setw(w7) << resp_fn
	    << " interpolant = "    << std::setw(w7) << interp_val
	    << " relative error = " << std::setw(w7) << err <<'\n';
      if (err > val_max_err) val_max_err = err;
      val_rmse += err * err;
      if (basisConfigOptions.useDerivs) {
	const RealVector& resp_grad   = surrData.response_gradient(index);
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
    val_rmse = std::sqrt(val_rmse/(num_colloc_pts-offset));
    PCout << "\nValue interpolation errors:    " << std::setw(w7) << val_max_err
	  << " (max) " << std::setw(w7) << val_rmse << " (RMS)\n";
    if (basisConfigOptions.useDerivs) {
      grad_rmse = std::sqrt(grad_rmse/(num_colloc_pts-offset)/num_v);
      PCout << "Gradient interpolation errors: " << std::setw(w7)
	    << grad_max_err << " (max) " << std::setw(w7) << grad_rmse
	    << " (RMS)\n";
    }
  }
#endif // INTERPOLATION_TEST
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
    SharedInterpPolyApproxData* data_rep
      = (SharedInterpPolyApproxData*)sharedDataRep;
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
	<< "sobolIndices =\n"; write_data(PCout, sobolIndices);
#endif // DEBUG
}


void InterpPolyApproximation::compute_total_sobol()
{
  totalSobolIndices = 0.; // init total indices

  SharedInterpPolyApproxData* data_rep
    = (SharedInterpPolyApproxData*)sharedDataRep;
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
	<< "totalSobolIndices =\n"; write_data(PCout, totalSobolIndices);
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

  SharedInterpPolyApproxData* data_rep
    = (SharedInterpPolyApproxData*)sharedDataRep;
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
