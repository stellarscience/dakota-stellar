/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 SparseGridDriver
//- Description: Implementation code for SparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "SparseGridDriver.hpp"
#include "PolynomialApproximation.hpp"
#include "sandia_sgmga.hpp"
#include "sandia_sgmgg.hpp"
#include "pecos_stat_util.hpp"
#include "MultivariateDistribution.hpp"

static const char rcsId[]="@(#) $Id: SparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void SparseGridDriver::assign_1d_collocation_points_weights()
{
  // resize arrays
  unsigned short ssg_lev = ssgLevIter->second;
  size_t i, num_levels = ssg_lev + 1, curr_lev;
  curr_lev = collocPts1D.size();
  if (num_levels > curr_lev) {
    collocPts1D.resize(num_levels);
    for (i=curr_lev; i<num_levels; ++i)
      collocPts1D[i].resize(numVars);
  }
  curr_lev = type1CollocWts1D.size();
  if (num_levels > curr_lev) {
    type1CollocWts1D.resize(num_levels);
    for (i=curr_lev; i<num_levels; ++i)
      type1CollocWts1D[i].resize(numVars);
  }
  curr_lev = type2CollocWts1D.size();
  if (computeType2Weights && num_levels > curr_lev) {
    type2CollocWts1D.resize(num_levels);
    for (i=curr_lev; i<num_levels; ++i)
      type2CollocWts1D[i].resize(numVars);
  }
  // assign values
  // level_index (j indexing) range is 0:w, level (i indexing) range is 1:w+1
  unsigned short l_index, q_order;
  for (i=0; i<numVars; i++)
    for (l_index=0; l_index<num_levels; ++l_index) {
      level_to_order(i, l_index, q_order);
      IntegrationDriver::
	assign_1d_collocation_points_weights(i, q_order, l_index);
    }
}


void SparseGridDriver::dimension_preference(const RealVector& dim_pref)
{
  RealVector aniso_wts;
  if (!dim_pref.empty()) {
    size_t num_pref = dim_pref.length();
    aniso_wts.sizeUninitialized(num_pref);
    webbur::sandia_sgmga_importance_to_aniso(num_pref, dim_pref.values(),
					     aniso_wts.values());
#ifdef DEBUG
    PCout << "dimension preference:\n" << dim_pref << "anisotropic weights "
	  << "after sandia_sgmga_importance_to_aniso():\n" << aniso_wts;
#endif
  }
  anisotropic_weights(aniso_wts);
}


void SparseGridDriver::anisotropic_weights(const RealVector& aniso_wts)
{
  RealVector& curr_aniso_wts = anisoWtsIter->second;
  if (aniso_wts.empty()) {
    if (!curr_aniso_wts.empty()) { // change from current
      curr_aniso_wts.sizeUninitialized(0);
      clear_grid(); // clear state to mandate a grid / grid size update
    }
  }
  else {
    if (aniso_wts.length() != numVars) {
      PCerr << "Error: length of sparse grid anisotropic weights specification "
	    << "is inconsistent with\n       number of variables in SparseGrid"
	    << "Driver::anisotropic_weights()." << std::endl;
      abort_handler(-1);
    }

    // detect anisotropy
    bool dim_iso = true;  size_t i;
    const Real& wt0 = aniso_wts[0];
    for (i=1; i<numVars; ++i)
      if (std::abs(aniso_wts[i] - wt0) > DBL_EPSILON)
	{ dim_iso = false; break; }
    // update active anisoLevelWts and grid update indicator
    if (dim_iso) {
      if (!curr_aniso_wts.empty()) {
	curr_aniso_wts.sizeUninitialized(0);
	clear_grid(); // clear state to mandate a grid / grid size update
      }
    }
    else {
      RealVector prev_aniso_wts = curr_aniso_wts; // for update indicator
      // truncate any negative values
      curr_aniso_wts.resize(numVars);
      for (i=0; i<numVars; ++i)
	curr_aniso_wts[i] = std::max(aniso_wts[i], 0.);
      // normalize and enforce axis lower bounds/weight upper bounds
      int option = 1; // weights scaled so that minimum nonzero entry is 1
      webbur::sandia_sgmga_aniso_normalize(option, numVars,
					   curr_aniso_wts.values());
#ifdef DEBUG
      PCout << "anisoLevelWts after sandia_sgmga_aniso_normalize():\n"
	    << curr_aniso_wts;
#endif
      // enforce axis lower bounds, if present, for current ssgLevel.  An axis
      // lower bound defines a weight upper bound based on the current ssgLevel:
      // LB_i = level*wt_min/wt_i --> wt_i = level*wt_min/LB_i and wt_min=1.
      // Catch special case of dim_pref_i = 0 --> wt_i = LB_i = 0.
      const RealVector& axis_l_bnds
	= axisLowerBounds[activeKey];//axisLBndsIter->second;
      if (!axis_l_bnds.empty()) {
	Real ssg_lev = (Real)(ssgLevIter->second);
	for (i=0; i<numVars; ++i)
	  if (axis_l_bnds[i] > SMALL_NUMBER) {                     // nonzero LB
	    Real wt_u_bnd = ssg_lev / axis_l_bnds[i];
	    curr_aniso_wts[i] = (curr_aniso_wts[i] > SMALL_NUMBER) // nonzero wt
	      ? std::min(wt_u_bnd, curr_aniso_wts[i]) : wt_u_bnd;
	  }
#ifdef DEBUG
	PCout << "anisoLevelWts after axisLowerBounds enforcement:\n"
	      << curr_aniso_wts;
#endif
      }
      // indicate need for numCollocPts update
      if (curr_aniso_wts != prev_aniso_wts)
	clear_grid(); // clear state to mandate a grid / grid size update
    }
  }
}


void SparseGridDriver::update_axis_lower_bounds()
{
  const RealVector& aniso_wts = anisoWtsIter->second;
  RealVector& axis_l_bnds = axisLowerBounds[activeKey];//axisLBndsIter->second;
  if (axis_l_bnds.empty())
    axis_l_bnds.sizeUninitialized(numVars);
  // An axis lower bound is the maximum index coverage achieved on a coordinate
  // axis (when all other indices are zero); it defines a constraint for
  // minimum coordinate coverage in future refinements.  The linear index set
  // constraint is level*wt_min-|wt| < j.wt <= level*wt_min, which becomes
  // level-|wt| < j_i w_i <= level for wt_min=1 and all other indices=0.
  // The max feasible j_i is then level/w_i (except for special case w_i=0).
  Real ssg_lev = (Real)(ssgLevIter->second);
  if (aniso_wts.empty())
    axis_l_bnds = ssg_lev; // all weights = 1
  else // min nonzero weight scaled to 1 --> just catch special case w_i=0
    for (size_t i=0; i<numVars; ++i)
      axis_l_bnds[i] = (aniso_wts[i] > SMALL_NUMBER) ? // nonzero wt
	ssg_lev / aniso_wts[i] : 0.;
}


void SparseGridDriver::
initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		const MultivariateDistribution& u_dist,
		const ExpansionConfigOptions& ec_options,
		BasisConfigOptions& bc_options, short growth_rate)
{
  growthRate             = growth_rate;
  //refineType           = ec_options.refineType;
  refineControl          = ec_options.refineControl;

  // For unrestricted exponential growth, use of nested rules is restricted
  // to uniform/normal in order to enforce similar growth rates:
  if (bc_options.nestedRules && growthRate == UNRESTRICTED_GROWTH) {
    const ShortArray&   u_types = u_dist.random_variable_types();
    const BitArray& active_vars = u_dist.active_variables();
    size_t i, num_u_types = u_types.size(); // numVars not yet defined
    bool no_mask = active_vars.empty();
    for (i=0; i<num_u_types; ++i)
      if ( ( no_mask || active_vars[i] ) &&
	   u_types[i] != STD_UNIFORM && u_types[i] != STD_NORMAL )
	{ bc_options.nestedRules = false; break; }
  }
  // For MODERATE and SLOW restricted exponential growth, nested rules
  // can be used heterogeneously and synchronized with STANDARD and SLOW
  // linear growth, respectively.

  IntegrationDriver::initialize_grid(u_dist, ec_options, bc_options);

  level(ssg_level);
  dimension_preference(dim_pref);
}


void SparseGridDriver::precompute_rules()
{
  unsigned short l, m, ssg_lev = ssgLevIter->second;
  if (isotropic())
    for (size_t i=0; i<numVars; ++i) {
      level_to_order(i, ssg_lev, m); // max order is full level in this dim
      polynomialBasis[i].precompute_rules(m);
    }
  else {
    const RealVector& aniso_wts = anisoWtsIter->second;
    Real wt_i;
    for (size_t i=0; i<numVars; ++i) {
      wt_i = aniso_wts[i];
      l = (wt_i > 0.) ? (unsigned short)((Real)ssg_lev / wt_i) : 0;
      level_to_order(i, l, m); // max order is full aniso level[dim]
      polynomialBasis[i].precompute_rules(m);
    }
  }
}


void SparseGridDriver::initialize_sets()
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "initialize_sets()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::increment_smolyak_multi_index(const UShortArray& set)
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "increment_smolyak_multi_index()." << std::endl;
  abort_handler(-1);
}


bool SparseGridDriver::
push_trial_available(const UShortArray& key, const UShortArray& tr_set)
{ return false; }


bool SparseGridDriver::push_trial_available(const UShortArray& key)
{ return false; }


bool SparseGridDriver::push_trial_available()
{ return false; }


size_t SparseGridDriver::
push_trial_index(const UShortArray& key, const UShortArray& tr_set)
{ return _NPOS; }


size_t SparseGridDriver::push_trial_index(const UShortArray& key)
{ return _NPOS; }


size_t SparseGridDriver::push_trial_index()
{ return _NPOS; }


size_t SparseGridDriver::push_index(const UShortArray& key) const
{ return _NPOS; }


size_t SparseGridDriver::restore_index(const UShortArray& key) const
{ return push_index(key); } // default for identity mapping (flat to flat)


size_t SparseGridDriver::finalize_index(size_t i, const UShortArray& key) const
{ return i; } // default is an identity mapping


void SparseGridDriver::push_set()
{
  PCerr << "Error: no default implementation for SparseGridDriver::push_set()."
	<< std::endl;
  abort_handler(-1);
}


void SparseGridDriver::pop_set()
{
  PCerr << "Error: no default implementation for SparseGridDriver::pop_set()."
	<< std::endl;
  abort_handler(-1);
}


void SparseGridDriver::
finalize_sets(bool output_sets, bool converged_within_tol, bool reverted)
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "finalize_sets()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::update_reference()
{
  // Not needed for HierarchSparseGridDriver, so use no-op as default

  /*
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "update_reference()." << std::endl;
  abort_handler(-1);
  */
}


void SparseGridDriver::compute_trial_grid(RealMatrix& var_sets)
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "compute_trial_grid()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::compute_increment(RealMatrix& var_sets)
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "compute_increment()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::push_increment()
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "push_increment()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::pop_increment()
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "pop_increment()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::merge_unique()
{ } // not needed for HierarchSparseGridDriver, so use no-op as default


const UShortArray& SparseGridDriver::trial_set(const UShortArray& key) const
{
  PCerr << "Error: no default implementation for SparseGridDriver::trial_set()."
	<< std::endl;
  abort_handler(-1);
  return key; // dummy UShortArray
}


const UShortArray& SparseGridDriver::trial_set() const
{ return trial_set(activeKey); } // default implementation


int SparseGridDriver::unique_trial_points() const
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "unique_trial_points()." << std::endl;
  abort_handler(-1);
  return 0;
}


void SparseGridDriver::update_smolyak_arrays()
{
  PCerr << "Error: no default implementation for SparseGridDriver::"
	<< "update_smolyak_arrays()." << std::endl;
  abort_handler(-1);
}


void SparseGridDriver::update_sets(const UShortArray& set_star)
{
  // set_star is passed as *cit_star from the best entry in activeMultiIndex.
  // Therefore, we must use caution in updates to activeMultiIndex that can
  // invalidate cit_star.

  // update smolyakMultiIndex (permanently, will not be popped)
  increment_smolyak_multi_index(set_star);
  push_set();     // calls increment_unique()      --> INC2
  merge_unique(); // promotes increment to new ref --> INC3

  // use trial set rather than incoming set_star due to iterator invalidation
  const UShortArray&    tr_set = trial_set();
  UShortArrayDeque& pop_trials =  poppedTrialSets[activeKey];
  UShortArraySet&    active_mi = activeMultiIndex[activeKey];
  UShortArraySet&       old_mi =    oldMultiIndex[activeKey];

  // update set O by adding the trial set to oldMultiIndex:
  old_mi.insert(tr_set);
  // remove the trial set from set A by erasing from activeMultiIndex:
  active_mi.erase(tr_set); // invalidates cit_star -> set_star
  // update subset of A that have been evaluated as trial sets but not selected
  UShortArrayDeque::iterator tr_it
    = std::find(pop_trials.begin(), pop_trials.end(), tr_set);
  if (tr_it != pop_trials.end()) pop_trials.erase(tr_it);

  // update set A (activeMultiIndex) based on neighbors of trial set
  add_active_neighbors(tr_set, false);//, isotropic());

  // Note: pruning irrelevant sets that have Coeff = 0 would be tricky,
  //       since a 0 close to the frontier can become nonzero

#ifdef DEBUG
  PCout << "Sets updated: (Smolyak,Old,Active,Trial) = (" << smolyak_size()
	<< ',' << old_mi.size() << ',' << active_mi.size() << ','
	<< pop_trials.size() << ')' << std::endl;
#endif // DEBUG
}


void SparseGridDriver::
add_active_neighbors(const UShortArray& set, bool frontier)
{
  UShortArray     trial_set =    set;
  UShortArraySet&    old_mi =    oldMultiIndex[activeKey];
  UShortArraySet& active_mi = activeMultiIndex[activeKey];
  UShortArraySet::const_iterator cit;
  size_t i, j, num_v = set.size();
  for (i=0; i<num_v; ++i) {
    // i^{th} candidate for set A (active) computed from forward neighbor:
    // increment by 1 in dimension i
    unsigned short& trial_set_i = trial_set[i];
    ++trial_set_i;
    // if !frontier, then candidates could exist in oldMultiIndex
    if (frontier || old_mi.find(trial_set) == old_mi.end()) {
      // test all backwards neighbors for membership in set O (old)
      bool backward_old = true;
      for (j=0; j<num_v; ++j) {
	unsigned short& trial_set_j = trial_set[j];
	if (trial_set_j) { // if 0, then admissible by default
	  --trial_set_j;
	  cit = old_mi.find(trial_set);
	  ++trial_set_j; // restore
	  if (cit == old_mi.end())
	    { backward_old = false; break; }
	}
      }
      if (backward_old) // std::set<> will discard any active duplicates
	active_mi.insert(trial_set);
    }
    --trial_set_i; // restore
  }
}


void SparseGridDriver::clear_inactive()
{
  // These are always defined in update_active_iterators()
  std::map<UShortArray, unsigned short>::iterator sg_it = ssgLevel.begin();
  std::map<UShortArray, RealVector>::iterator     aw_it = anisoLevelWts.begin();
  std::map<UShortArray, int>::iterator            cp_it = numCollocPts.begin();
  while (sg_it != ssgLevel.end())
    if (sg_it == ssgLevIter) // preserve active
      { ++sg_it; ++aw_it; ++cp_it; }
    else { // clear inactive: postfix increments manage iterator invalidations
      ssgLevel.erase(sg_it++);
      anisoLevelWts.erase(aw_it++);
      numCollocPts.erase(cp_it++);
    }

  // Generalized sparse grid sets may be active
  if (!oldMultiIndex.empty()) {
    std::map<UShortArray, UShortArraySet>::iterator
      om_it = oldMultiIndex.begin(), om_act_it = oldMultiIndex.find(activeKey);
    std::map<UShortArray, UShortArraySet>::iterator am_it
      = activeMultiIndex.begin();
    std::map<UShortArray, UShortArrayDeque>::iterator pt_it
      = poppedTrialSets.begin();
    while (om_it != oldMultiIndex.end())
      if (om_it == om_act_it) // preserve active
	{ ++om_it; ++am_it; ++pt_it; }
      else { // clear inactive: postfix increments manage iterator invalidations
	oldMultiIndex.erase(om_it++);
	activeMultiIndex.erase(am_it++);
	poppedTrialSets.erase(pt_it++);
      }
  }

  // Anisotropic refinement bounds may be active
  if (!axisLowerBounds.empty()) {
    std::map<UShortArray, RealVector>::iterator ab_it = axisLowerBounds.begin(),
      ab_act_it = axisLowerBounds.find(activeKey);
    while (ab_it != axisLowerBounds.end())
      if (ab_it == ab_act_it) // preserve active
	++ab_it;
      else // clear inactive: postfix increments manage iterator invalidations
	axisLowerBounds.erase(ab_it++);
  }
}


int SparseGridDriver::level_to_order_exp_hgk_interp(int level, int growth)
{
  if (level == 0) return 1;

  switch (growth) {
  case SLOW_RESTRICTED_GROWTH: {
    unsigned short level = 0, max_level = 5;
    int m = 1, m_goal = level + 1;
    while (level <= max_level && m < m_goal)
      { ++level; m = orderGenzKeister[level]; }
    return m; break;
  }
  case MODERATE_RESTRICTED_GROWTH: {
    unsigned short level = 0, max_level = 5;
    int m = 1, m_goal = 2 * level + 1;
    while (level <= max_level && m < m_goal)
      { ++level; m = orderGenzKeister[level]; }
    return m; break;
  }
  case UNRESTRICTED_GROWTH:
    return orderGenzKeister[std::min(level, 5)]; break;
  }
}


int SparseGridDriver::level_to_order_exp_closed_interp(int level, int growth)
{
  if (level == 0) return 1;

  switch (growth) {
  case SLOW_RESTRICTED_GROWTH: {
    int m = 1, lev_pow = 1, m_goal = level + 1;
    while (m < m_goal)
      { lev_pow *= 2; m = lev_pow + 1; } //std::pow(2.,level) + 1;
    return m; break;
  }
  case MODERATE_RESTRICTED_GROWTH: {
    int m = 1, lev_pow = 1, m_goal = 2 * level + 1;
    while (m < m_goal)
      { lev_pow *= 2; m = lev_pow + 1; } //std::pow(2.,level) + 1;
    return m; break;
  }
  case UNRESTRICTED_GROWTH:
    return (int)std::pow(2., level) + 1; break;
  }
}


int SparseGridDriver::level_to_order_exp_open_interp(int level, int growth)
{
  if (level == 0) return 1;

  switch (growth) {
  case SLOW_RESTRICTED_GROWTH: {
    int m = 1, lev_pow = 2, m_goal = level + 1;
    while (m < m_goal)
      { lev_pow *= 2; m = lev_pow - 1; } //std::pow(2.,level+1) - 1;
    return m; break;
  }
  case MODERATE_RESTRICTED_GROWTH: {
    int m = 1, lev_pow = 2, m_goal = 2 * level + 1;
    while (m < m_goal)
      { lev_pow *= 2; m = lev_pow - 1; } //std::pow(2.,level+1) - 1;
    return m; break;
  }
  case UNRESTRICTED_GROWTH:
    return (int)std::pow(2., level+1) - 1; break;
  }
}

} // namespace Pecos
