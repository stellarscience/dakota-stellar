/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "pecos_stat_util.hpp"
#include "RandomVariable.hpp"

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
