/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "MarginalsCorrDistribution.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: MarginalsCorrDistribution.cpp 4768 2007-12-17 17:49:32Z mseldre $";

//#define DEBUG


namespace Pecos {

void MarginalsCorrDistribution::
initialize_types(const ShortArray& rv_types, const BitArray& active_vars)
{
  ranVarTypes = rv_types;
  activeVars  = active_vars;

  // construction of x-space random variables occurs once (updates to
  // distribution parameters can occur repeatedly)
  size_t i, num_v = rv_types.size();
  randomVars.resize(num_v);
  for (i=0; i<num_v; ++i)
    randomVars[i] = RandomVariable(rv_types[i]);

  // assign globalBndsFlag
  check_global_bounds();
}


void MarginalsCorrDistribution::
initialize_correlations(const RealSymMatrix& corr, const BitArray& active_corr)
{
  corrMatrix = corr;
  activeCorr = active_corr; // active RV subset for correlation matrix

  initialize_correlations();
}


void MarginalsCorrDistribution::initialize_correlations()
{
  correlationFlag = false;
  size_t num_corr = corrMatrix.numRows();
  if (num_corr == 0) return;

  size_t num_rv = randomVars.size();
  bool no_mask = activeCorr.empty();
  if (no_mask) {
    if (num_corr != num_rv) {
      PCerr << "Error: correlation matrix size (" << num_corr
	    << ") inconsistent with number of random variables (" << num_rv
	    << ")." << std::endl;
      abort_handler(-1);
    }
  }
  else {
    if (num_corr != activeCorr.count()) {
      PCerr << "Error: correlation matrix size (" << num_corr
	    << ") inconsistent with active correlation subset ("
	    << activeCorr.count() << ")." << std::endl;
      abort_handler(-1);
    }
  }

  size_t i, j, cntr_i, cntr_j;
  for (i=0, cntr_i=0; i<num_rv; ++i) {
    if (no_mask || activeCorr[i]) {
      for (j=0, cntr_j=0; j<i; ++j) {
	if (no_mask || activeCorr[j]) {
	  if (std::abs(corrMatrix(cntr_i, cntr_j)) > SMALL_NUMBER)
	    correlationFlag = true;
	  ++cntr_j;
	}
      }
      ++cntr_i;
    }
    if (correlationFlag) break;
  }
}


/** For general random variable ordering (e.g., NestedModel mappings). */
void MarginalsCorrDistribution::
pull_distribution_parameters
(const std::shared_ptr<MultivariateDistribution> pull_mvd_rep,
 size_t pull_index, size_t push_index)
{
  RandomVariable&       push_rv = randomVars[push_index];
  const RandomVariable& pull_rv = pull_mvd_rep->random_variable(pull_index);
  short push_type = ranVarTypes[push_index],
        pull_type = pull_mvd_rep->random_variable_type(pull_index);
  switch (push_type) {

  // push RV are fully standardized: no updates to perform
  case STD_NORMAL:  case STD_UNIFORM:  case STD_EXPONENTIAL:  break;

  // push RV have standardized scale params; copy shape params
  case STD_BETA: {
    Real alpha;  pull_rv.pull_parameter(BE_ALPHA, alpha);
    Real beta;   pull_rv.pull_parameter(BE_BETA,  beta);
    push_rv.push_parameter(BE_ALPHA, alpha);
    push_rv.push_parameter(BE_BETA,  beta);
    break;
  }
  case STD_GAMMA: {
    Real alpha;  pull_rv.pull_parameter(GA_ALPHA, alpha);
    push_rv.push_parameter(GA_ALPHA, alpha);
    break;
  }

  // push RV are non-standard, pull all non-standard data from rv_in
  default:
    switch (pull_type) {
    // pull RV are fully standardized: no data to pull
    case STD_NORMAL:  case STD_UNIFORM:  case STD_EXPONENTIAL:  break;

    // pull RV have standardized scale params; pull shape params
    case STD_BETA: {
      Real alpha;  pull_rv.pull_parameter(BE_ALPHA, alpha);
      Real beta;   pull_rv.pull_parameter(BE_BETA,  beta);
      push_rv.push_parameter(BE_ALPHA, alpha);
      push_rv.push_parameter(BE_BETA,  beta);
      break;
    }
    case STD_GAMMA: {
      Real alpha;  pull_rv.pull_parameter(GA_ALPHA, alpha);
      push_rv.push_parameter(GA_ALPHA, alpha);
      break;
    }

    // pull and push RV are non-standardized; copy all params
    default:
      push_rv.copy_parameters(pull_rv);                         break;
    }
    break;
  }
}


void MarginalsCorrDistribution::
copy_rep(std::shared_ptr<MultivariateDistribution> source_rep)
{
  // copy base class data
  MultivariateDistribution::copy_rep(source_rep);
  // specialization for marginals + corr
  std::shared_ptr<MarginalsCorrDistribution> mcd_rep =
    std::static_pointer_cast<MarginalsCorrDistribution>(source_rep);
  initialize_types(mcd_rep->ranVarTypes, mcd_rep->activeVars);
  initialize_correlations(mcd_rep->corrMatrix, mcd_rep->activeCorr);
  pull_distribution_parameters(source_rep);
}


/* Reshaping the correlation matrix should no longer be required
   (subsets now supported with activeCorr BitArray)
void MarginalsCorrDistribution::
expand_correlation_matrix(size_t num_lead_v, size_t num_prob_v,
			  size_t num_trail_v)
{
  if (!correlationFlag)
    return;
  size_t num_prev_v  = corrMatrix.numRows(),
         num_total_v = num_lead_v + num_prob_v + num_trail_v;
  if (num_prev_v == num_total_v)
    return;
  else if (num_prev_v != num_prob_v) { // old->new subset assumed below
    PCerr << "\nError: unexpected matrix size (" << num_prev_v
	  << ") in MarginalsCorrDistribution::expand_correlation_matrix()."
	  << std::endl;
    abort_handler(-1);
  }

  // expand from num_prob_v to num_total_v
  // Note: a reshape is not helpful due to num_lead_v

  size_t i, j, offset;
  RealSymMatrix old_corr_matrix(corrMatrix); // copy
  corrMatrix.shape(num_total_v); // init to zero
  for (i=0; i<num_lead_v; i++)
    corrMatrix(i,i) = 1.;
  offset = num_lead_v;
  for (i=0; i<num_prob_v; i++)
    for (j=0; j<=i; j++)
      corrMatrix(i+offset, j+offset) = old_corr_matrix(i,j);
  offset += num_prob_v;
  for (i=0; i<num_trail_v; i++)
    corrMatrix(i+offset, i+offset) = 1.;
}


void MarginalsCorrDistribution::
contract_correlation_matrix(size_t num_lead_v, size_t num_prob_v,
			    size_t num_trail_v)
{
  if (!correlationFlag)
    return;
  size_t num_prev_v  = corrMatrix.numRows(),
         num_total_v = num_lead_v + num_prob_v + num_trail_v;
  if (num_prev_v == num_prob_v)
    return;
  else if (num_prev_v != num_total_v) {
    PCerr << "\nError: unexpected matrix size (" << num_prev_v
	  << ") in MarginalsCorrDistribution::contract_correlation_matrix()."
	  << std::endl;
    abort_handler(-1);
  }

  // contract from num_total_v to num_prob_v
  // Note: a reshape is not helpful due to num_lead_v

  size_t i, j;
  RealSymMatrix old_corr_matrix(corrMatrix); // copy
  corrMatrix.shape(num_prob_v); // init to zero
  for (i=0; i<num_prob_v; i++)
    for (j=0; j<=i; j++)
      corrMatrix(i, j) = old_corr_matrix(i+num_lead_v,j+num_lead_v);
}
*/

} // namespace Pecos
