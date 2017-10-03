/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 HierarchSparseGridDriver
//- Description: Implementation code for HierarchSparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "HierarchSparseGridDriver.hpp"
#include "SharedPolyApproxData.hpp"
#include "sandia_sgmga.hpp"
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: HierarchSparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void HierarchSparseGridDriver::
initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		const ShortArray& u_types,
		const ExpansionConfigOptions& ec_options,
		BasisConfigOptions& bc_options, short growth_rate,
		bool track_colloc_indices)
{
  SparseGridDriver::initialize_grid(ssg_level, dim_pref, u_types, ec_options,
				    bc_options, growth_rate);
  trackCollocIndices = track_colloc_indices;
}


void HierarchSparseGridDriver::store_grid(size_t index)
{
  size_t stored_len = storedType1WeightSets.size();
  if (index == _NPOS || index == stored_len) { // append
    storedLevMultiIndex.push_back(smolyakMultiIndex);
    storedCollocKey.push_back(collocKey);
    //storedCollocIndices.push_back(collocIndices);
    storedType1WeightSets.push_back(type1WeightSets);
    storedType2WeightSets.push_back(type2WeightSets);
  }
  else if (index < stored_len) { // replace
    storedLevMultiIndex[index]   = smolyakMultiIndex;
    storedCollocKey[index]       = collocKey;
    //storedCollocIndices[index] = collocIndices;
    storedType1WeightSets[index] = type1WeightSets;
    storedType2WeightSets[index] = type2WeightSets;
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in "
	  << "HierarchSparseGridDriver::store_grid()" << std::endl;
    abort_handler(-1);
  }
}


void HierarchSparseGridDriver::restore_grid(size_t index)
{
  size_t stored_len = storedType1WeightSets.size();
  if (index == _NPOS) {
    smolyakMultiIndex = storedLevMultiIndex.back();
    collocKey       = storedCollocKey.back();
    //collocIndices = storedCollocIndices.back();
    type1WeightSets = storedType1WeightSets.back();
    type2WeightSets = storedType2WeightSets.back();
  }
  else if (index < stored_len) {
    smolyakMultiIndex = storedLevMultiIndex[index];
    collocKey       = storedCollocKey[index];
    //collocIndices = storedCollocIndices[index];
    type1WeightSets = storedType1WeightSets[index];
    type2WeightSets = storedType2WeightSets[index];
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in "
	  << "HierarchSparseGridDriver::restore_grid()" << std::endl;
    abort_handler(-1);
  }
}


void HierarchSparseGridDriver::remove_stored_grid(size_t index)
{
  size_t stored_len = storedType1WeightSets.size();
  if (index == _NPOS || index == stored_len) {
    storedLevMultiIndex.pop_back();
    storedCollocKey.pop_back();
    //storedCollocIndices.pop_back();
    storedType1WeightSets.pop_back();
    storedType2WeightSets.pop_back();
  }
  else if (index < stored_len) {
    UShort4DArray::iterator u4it = storedLevMultiIndex.begin();
    std::advance(u4it, index); storedLevMultiIndex.erase(u4it);
    UShort5DArray::iterator u5it = storedCollocKey.begin();
    std::advance(u5it, index); storedCollocKey.erase(u5it);
    //Sizet4DArray::iterator s4it = storedCollocIndices.begin();
    //std::advance(s4it, index); storedCollocIndices.erase(s4it);
    RealVector3DArray::iterator vit = storedType1WeightSets.begin();
    std::advance(vit, index); storedType1WeightSets.erase(vit);
    RealMatrix3DArray::iterator mit = storedType2WeightSets.begin();
    std::advance(mit, index); storedType2WeightSets.erase(mit);
  }
}


void HierarchSparseGridDriver::clear_stored()
{
  storedLevMultiIndex.clear(); storedCollocKey.clear();
  //storedCollocIndices.clear();

  storedType1WeightSets.clear();
  storedType2WeightSets.clear();
}


size_t HierarchSparseGridDriver::maximal_grid() const
{
  size_t i, j, k, num_stored = storedType1WeightSets.size(), max_index = _NPOS,
    num_max_wts = 0, num_stored_wts, num_lev = type1WeightSets.size(), num_sets;
  for (j=0; j<num_lev; ++j) {
    const RealVectorArray& curr_wts_j = type1WeightSets[j];
    num_sets = curr_wts_j.size();
    for (k=0; k<num_sets; ++k)
      num_max_wts += curr_wts_j[k].length();
  }
  for (i=0; i<num_stored; ++i) {
    const RealVector2DArray& stored_wts = storedType1WeightSets[i];
    num_lev = stored_wts.size(); num_stored_wts = 0;
    for (j=0; j<num_lev; ++j) {
      const RealVectorArray& stored_wts_j = stored_wts[j];
      num_sets = stored_wts_j.size();
      for (k=0; k<num_sets; ++k)
	num_stored_wts += stored_wts_j[k].length();
    }
    if (num_stored_wts > num_max_wts)
      { max_index = i; num_max_wts = num_stored_wts; }
  }
  return max_index;
}


void HierarchSparseGridDriver::swap_grid(size_t index)
{
  std::swap(storedLevMultiIndex[index], smolyakMultiIndex);
  std::swap(storedCollocKey[index],     collocKey);
  //std::swap(storedCollocIndices[index], collocIndices);

  std::swap(storedType1WeightSets[index], type1WeightSets);
  std::swap(storedType2WeightSets[index], type2WeightSets);
}


int HierarchSparseGridDriver::grid_size()
{
  if (updateGridSize) {
    numCollocPts = 0;
    if (collocKey.size() == ssgLevel+1) { // collocKey already up to date
      for (unsigned short i=0; i<=ssgLevel; ++i) {
	const UShort3DArray& key_i = collocKey[i];
	size_t j, num_sets = key_i.size();
	for (j=0; j<num_sets; ++j)
	  numCollocPts += key_i[j].size(); // hierarchical point increments
      }
    }
    else {
      update_smolyak_multi_index();
      // rather than full collocKey update, just sum grid sizes:
      UShortArray delta_sizes(numVars);
      unsigned short lev, set, num_sets;
      for (lev=0; lev<=ssgLevel; ++lev) {
	UShort2DArray& sm_mi_l = smolyakMultiIndex[lev];
	num_sets = sm_mi_l.size();
	for (set=0; set<num_sets; ++set) {
	  levels_to_delta_sizes(sm_mi_l[set], delta_sizes);
	  numCollocPts +=
	    SharedPolyApproxData::tensor_product_terms(delta_sizes, false);
	}
      }
    }
    updateGridSize = false;
  }
  return numCollocPts;
}


/*
size_t HierarchSparseGridDriver::smolyak_size()
{
  size_t num_sm_mi = 0;
  unsigned short i, num_lev = smolyakMultiIndex.size();
  for (i=0; i<num_lev; ++i)
    num_sm_mi += smolyakMultiIndex[i].size();
  return num_sm_mi;
}
*/


void HierarchSparseGridDriver::update_smolyak_multi_index(bool clear_sm_mi)
{
  if (clear_sm_mi)
    smolyakMultiIndex.clear();

  size_t prev_sm_len = smolyakMultiIndex.size();
  // anisotropic weight updates always accompanied by a level increment, so
  // we are already up to date for both iso and aniso cases if equal:
  if (prev_sm_len == ssgLevel+1)
    return;

  // this function is for use with isotropic/anisotropic grids, including
  // the initial starting point for a generalized sparse grid adaptation
  bool from_scratch = (prev_sm_len == 0);
  if (!from_scratch && refineControl == DIMENSION_ADAPTIVE_CONTROL_GENERALIZED){
    PCerr << "Error: HierarchSparseGridDriver::update_smolyak_multi_index() "
	  << "intended for use with isotropic and anisotropic grid refinements."
	  << std::endl;
    abort_handler(-1);
  }

  // Populate smolyakMultiIndex based on {iso,aniso}tropic index set constraint:
  // w = q - N = dimension-independent level.  For isotropic,
  //   w + 1 <= |i| <= w + N for i starts at 1 (used for index set defn.)
  //   w - N + 1 <= |j| <= w for j = i - 1 starts at 0 (used for generation)
  // For anisotropic, a weighted linear index set constraint is used.

  size_t lev;
  smolyakMultiIndex.resize(ssgLevel+1);
  if (dimIsotropic)
    for (lev=prev_sm_len; lev<=ssgLevel; ++lev)
      SharedPolyApproxData::total_order_multi_index(lev, numVars,
						    smolyakMultiIndex[lev]);
  else { // utilize webbur::sandia_sgmga_vcn_ordered

    // With scaling alpha_min = 1: q_min < |alpha . j| <= q_max.
    IntArray x(numVars), x_max(numVars); //x_max = ssgLevel;
    Real q_min = -1., q_max = ssgLevel; // no lower bound for hierarchical
    for (size_t i=0; i<numVars; ++i) {
      const Real& wt_i = anisoLevelWts[i];
      // minimum nonzero weight is scaled to 1, so just catch special case of 0
      x_max[i] = (wt_i > 1.e-10) ? (int)std::ceil(q_max/wt_i) : 0;
    }

    // We take the simple approach and restart vcn iteration from scratch,
    // checking each index_set for current inclusion in smolyakMultiIndex.
    // Warm starting the vcn iteration from the initial smolyakMultiIndex
    // state would be a more efficient option, although also more complex.
    bool more = false;
    UShortArray index_set(numVars);
    Real *aniso_wts = anisoLevelWts.values();
    int  *x0 = &x[0], *xm0 = &x_max[0];
    webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0, q_min,
				     q_max, &more);
    while (more) {
      for (size_t i=0; i<numVars; ++i)
	index_set[i] = (unsigned short)x[i];
      lev = l1_norm(index_set);
      UShort2DArray& sm_mi_l = smolyakMultiIndex[lev];
      if (from_scratch ||
	  std::find(sm_mi_l.begin(), sm_mi_l.end(), index_set) == sm_mi_l.end())
	sm_mi_l.push_back(index_set);
      webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0, q_min,
				       q_max, &more);
    }
  }

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::update_smolyak_multi_index():\n";
  size_t set, num_sets;
  for (lev=0; lev<=ssgLevel; ++lev) {
    num_sets = smolyakMultiIndex[lev].size();
    for (set=0; set<num_sets; ++set)
      PCout << "Smolyak multi_index[" << lev << "][" << set << "]:\n"
	    << smolyakMultiIndex[lev][set];
  }
#endif // DEBUG
}


void HierarchSparseGridDriver::assign_collocation_key()
{
  if (collocKey.size() == ssgLevel+1)
    return;

  collocKey.resize(ssgLevel+1);
  if (nestedGrid) {
    size_t lev, set, num_sets;
    UShort2DArray delta_keys(numVars);
    for (lev=0; lev<=ssgLevel; ++lev) {
      const UShort2DArray& sm_mi_l = smolyakMultiIndex[lev];
      UShort3DArray&         key_l = collocKey[lev];
      num_sets = sm_mi_l.size();
      key_l.resize(num_sets);
      for (set=0; set<num_sets; ++set) {
	levels_to_delta_keys(sm_mi_l[set], delta_keys);
	SharedPolyApproxData::hierarchical_tensor_product_multi_index(
	  delta_keys, key_l[set]);
      }
    }
  }
  //else
  //  SparseGridDriver::assign_collocation_key();

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::assign_collocation_key():\n";
  size_t lev, set, pt, num_sets, num_tp_pts;
  for (lev=0; lev<=ssgLevel; ++lev) {
    num_sets = collocKey[lev].size();
    for (set=0; set<num_sets; ++set) {
      num_tp_pts = collocKey[lev][set].size();
      for (pt=0; pt<num_tp_pts; ++pt)
	PCout << "Collocation key[" << lev << "][" << set << "][" << pt
	      << "]:\n" << collocKey[lev][set][pt];
    }
  }
#endif // DEBUG
}


void HierarchSparseGridDriver::update_collocation_key()
{
  size_t sm_mi_len = smolyakMultiIndex.size();
  if (collocKey.size() < sm_mi_len)
    collocKey.resize(sm_mi_len);
  UShort2DArray delta_keys(numVars), key_ls;
  if (refineControl == DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
    levels_to_delta_keys(trial_set(), delta_keys);

    UShort3DArray& key_l = collocKey[trialLevel];
    size_t set = key_l.size();
    key_l.push_back(key_ls); // update in place
    SharedPolyApproxData::hierarchical_tensor_product_multi_index(delta_keys,
								  key_l[set]);

#ifdef DEBUG
    PCout << "HierarchSparseGridDriver::update_collocation_key():\n";
    size_t pt, num_tp_pts = key_l[set].size();
    for (pt=0; pt<num_tp_pts; ++pt)
      PCout << "Collocation key[" << trialLevel << "][" << set << "][" << pt
	    << "]:\n" << key_l[set][pt];
#endif // DEBUG
  }
  else { // isotropic and anisotropic grid refinements
    // define incrementSets to track iso/aniso grid refinement increment
    size_t lev, set, start_set, num_sets;
    incrementSets.resize(ssgLevel+1);
    for (lev=0; lev<=ssgLevel; ++lev)
      incrementSets[lev] = collocKey[lev].size();
    // update collocKey to correspond to smolyakMultiIndex
    for (lev=0; lev<=ssgLevel; ++lev) {
      UShort2DArray& sm_mi_l = smolyakMultiIndex[lev];
      UShort3DArray&   key_l = collocKey[lev];
      start_set = incrementSets[lev]; num_sets = sm_mi_l.size();
      for (set=start_set; set<num_sets; ++set) {
	levels_to_delta_keys(sm_mi_l[set], delta_keys);
	key_l.push_back(key_ls); // update in place
	SharedPolyApproxData::hierarchical_tensor_product_multi_index(
	  delta_keys, key_l[set]);
      }
    }

#ifdef DEBUG
    PCout << "HierarchSparseGridDriver::update_collocation_key():\n";
    size_t pt, num_tp_pts;
    for (lev=0; lev<=ssgLevel; ++lev) {
      start_set = incrementSets[lev]; num_sets = collocKey[lev].size();
      for (set=start_set; set<num_sets; ++set) {
	num_tp_pts = collocKey[lev][set].size();
	for (pt=0; pt<num_tp_pts; ++pt)
	  PCout << "Collocation key[" << lev << "][" << set << "][" << pt
		<< "]:\n" << collocKey[lev][set][pt];
      }
    }
#endif // DEBUG
  }
}


void HierarchSparseGridDriver::assign_collocation_indices()
{
  if (collocIndices.size() == ssgLevel+1)
    return;

  collocIndices.resize(ssgLevel+1);
  size_t lev, set, pt, cntr = 0, num_sets, num_tp_pts;
  for (lev=0; lev<=ssgLevel; ++lev) {
    const UShort3DArray& key_l = collocKey[lev];
    num_sets = key_l.size();
    Sizet2DArray& indices_l = collocIndices[lev];
    indices_l.resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      const UShort2DArray& key_ls = key_l[set];
      num_tp_pts = key_ls.size();
      SizetArray& indices_ls = indices_l[set];
      indices_ls.resize(num_tp_pts);
      for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
	indices_ls[pt] = cntr; // simple sequential ordering for unique pts
    }
  }
  numCollocPts = cntr;

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::assign_collocation_indices():\n"
	<< "numCollocPts = " << numCollocPts << '\n';
  for (lev=0; lev<=ssgLevel; ++lev) {
    num_sets = collocIndices[lev].size();
    for (set=0; set<num_sets; ++set)
      PCout << "Collocation indices[" << lev << "][" << set << "]:\n"
	    << collocIndices[lev][set];
  }
#endif // DEBUG
}


void HierarchSparseGridDriver::update_collocation_indices()
{
  size_t cntr = numCollocPts, num_lev = collocKey.size();
  if (collocIndices.size() < num_lev)
    collocIndices.resize(num_lev);

  if (refineControl == DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
    size_t pt, num_tp_pts = collocKey[trialLevel].back().size();
    Sizet2DArray& indices_l = collocIndices[trialLevel];
    SizetArray indices; indices_l.push_back(indices); // update in place
    SizetArray& trial_indices = indices_l.back();
    trial_indices.resize(num_tp_pts);
    for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
      trial_indices[pt] = cntr;
    numCollocPts += num_tp_pts;

#ifdef DEBUG
    PCout << "HierarchSparseGridDriver::update_collocation_indices():\n"
	  << "numCollocPts = " << numCollocPts << '\n';
    size_t set = collocIndices[trialLevel].size() - 1;
    PCout << "Collocation indices[" << trialLevel << "][" << set << "]:\n"
	  << trial_indices;
#endif // DEBUG
  }
  else { // isotropic and anisotropic grid refinements
    size_t lev, set, start_set, num_sets, pt, num_tp_pts; SizetArray indices;
    for (lev=0; lev<num_lev; ++lev) {
      UShort2DArray&  sm_mi_l = smolyakMultiIndex[lev];
      UShort3DArray&    key_l = collocKey[lev];
      Sizet2DArray& indices_l = collocIndices[lev];
      start_set = incrementSets[lev]; num_sets = sm_mi_l.size();
      for (set=start_set; set<num_sets; ++set) {
	indices_l.push_back(indices); // update in place
	SizetArray& trial_indices = indices_l.back();
	num_tp_pts = collocKey[lev][set].size();
	trial_indices.resize(num_tp_pts);
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
	  trial_indices[pt] = cntr;
	numCollocPts += num_tp_pts;
      }
    }

#ifdef DEBUG
    PCout << "HierarchSparseGridDriver::update_collocation_indices():\n"
	  << "numCollocPts = " << numCollocPts << '\n';
    for (lev=0; lev<=ssgLevel; ++lev) {
      start_set = incrementSets[lev]; num_sets = collocIndices[lev].size();
      for (set=start_set; set<num_sets; ++set)
	PCout << "Collocation indices[" << lev << "][" << set << "]:\n"
	      << collocIndices[lev][set];
    }
#endif // DEBUG
  }
}


unsigned short HierarchSparseGridDriver::
level_to_delta_size(size_t i, unsigned short level)
{
  switch (level) { // growth restriction should not occur for lev_i = 0 or 1
  case 0: return 1;   break; // 1 new pt
  case 1: return 2;   break; // new ends of 3 pt (open or closed)
  default: { // difference point counts for level & level-1
    unsigned short num_delta, ord_lm1;
    level_to_order(i, level, num_delta); level_to_order(i, level-1, ord_lm1);
    num_delta -= ord_lm1; // Note: num_delta = 0 in case of growth restriction
    return num_delta; break;
  }
  }
}


void HierarchSparseGridDriver::
level_to_delta_key(size_t i, unsigned short lev_i, UShortArray& delta_key_i)
{
  unsigned short num_delta = level_to_delta_size(i, lev_i);
  delta_key_i.resize(num_delta);
  switch(collocRules[i]) {
  case GAUSS_PATTERSON: // open nested
    for (size_t j=0; j<num_delta; ++j)
      delta_key_i[j] = 2*j; // new ends + new interior: 0,2,4,6,8,...
    break;
  case NEWTON_COTES: case CLENSHAW_CURTIS: // closed nested
    switch (lev_i) { // growth restriction should not occur for lev_i = 0 or 1
    case 0: delta_key_i[0] = 0;                      break; // center of 1 pt
    case 1: delta_key_i[0] = 0; delta_key_i[1] = 2; break; // ends of 3 pt
    default:
      for (size_t j=0; j<num_delta; ++j)
	delta_key_i[j] = 2*j+1; // new interior: 1,3,5,7,9,...
      break;
    }
    break;
  case GENZ_KEISTER: // open nested table lookup
    // switch on num_delta i/o lev_i due to possibility of growth restriction
    switch (num_delta) {
    case 1: delta_key_i[0] = 0;                     break; // center of 1 pt
    case 2: delta_key_i[0] = 0; delta_key_i[1] = 2; break; // ends of 3 pt
    case 6:
      delta_key_i[0] = 0; delta_key_i[1] = 1; delta_key_i[2] = 3;
      delta_key_i[3] = 5; delta_key_i[4] = 7; delta_key_i[5] = 8;
      break; // 9 pt rule reusing 2, 4 (center), 6
    case 10:
      delta_key_i[0] =  0; delta_key_i[1] =  1; delta_key_i[2] =  3;
      delta_key_i[3] =  5; delta_key_i[4] =  7; delta_key_i[5] = 11;
      delta_key_i[6] = 13; delta_key_i[7] = 15; delta_key_i[8] = 17;
      delta_key_i[9] = 18;
      break; // 19 pt rule reusing 
    case 16:
      delta_key_i[0]  =  0; delta_key_i[1]  =  1; delta_key_i[2]  =  2;
      delta_key_i[3]  =  4; delta_key_i[4]  =  6; delta_key_i[5]  =  8;
      delta_key_i[6]  = 12; delta_key_i[7]  = 16; delta_key_i[8]  = 18;
      delta_key_i[9]  = 22; delta_key_i[10] = 26; delta_key_i[11] = 28;
      delta_key_i[12] = 30; delta_key_i[13] = 32; delta_key_i[14] = 33;
      delta_key_i[15] = 34;
      break; // 35 pt rule reusing
      // 3,5,7,9,10,11,13,14,15,17,19,20,21,23,24,25,27,29,31
    //case 5:  // 43 pt rule augments 19 pt rule, not 35 pt rule
    //  break; // disallow for hierarchical interpolation
    default:
      PCerr << "Error: out of range for hierarchical Genz-Keister rules in "
	    << "HierarchSparseGridDriver::level_to_delta_key()"
	    << std::endl;
      abort_handler(-1);
      break;
    }
    break;
  default:
    PCerr << "Error: bad rule type in level_to_delta_key()" << std::endl;
    abort_handler(-1);
    break;
  }
}


UShortUShortPair HierarchSparseGridDriver::
level_to_delta_pair(size_t i, unsigned short level)
{
  switch (level) { // growth restriction should not occur for level = 0 or 1
  case 0: // +1 pt,  max index = 0 for 1 pt rule
    return UShortUShortPair(1,0); break;
  case 1: // +2 pts, max index = 2 for right end of 3 pt rule (open or closed)
    return UShortUShortPair(2,2); break;
  default: {
    unsigned short max_key, num_delta = level_to_delta_size(i, level);
    switch(collocRules[i]) {
    case GAUSS_PATTERSON:                    // open nested
      max_key = 2 * num_delta - 2; break; // new exterior
    case NEWTON_COTES: case CLENSHAW_CURTIS: // closed nested
      max_key = 2 * num_delta - 1; break; // new interior
    case GENZ_KEISTER: // open nested table lookup
      // switch on num_delta i/o level due to possibility of growth restriction
      switch (num_delta) {
      case  6: max_key =  8; break; // 9 pt rule
      case 10: max_key = 18; break; // 19 pt rule
      case 16: max_key = 34; break; // 35 pt rule
      //case 5:  // 43 pt rule augments 19 pt rule, not 35 pt rule
      //  break; // disallow for hierarchical interpolation
      default:
	PCerr << "Error: out of range for hierarchical Genz-Keister rules in "
	      << "HierarchSparseGridDriver::level_to_delta_pair()" << std::endl;
	abort_handler(-1);
	break;
      }
      break;
    default:
      PCerr << "Error: bad rule type in level_to_delta_pair()" << std::endl;
      abort_handler(-1);
      break;
    }
    return UShortUShortPair(num_delta, max_key);
    break;
  }
  }
}


void HierarchSparseGridDriver::compute_grid(RealMatrix& var_sets)
{
  bool clear = (refineControl != NO_CONTROL); // restore prev state if refined
  update_smolyak_multi_index(clear);          // compute smolyakMultiIndex
  assign_collocation_key();                   // compute collocKey
  assign_1d_collocation_points_weights();     // define 1-D point/weight sets

  // For efficiency reasons, incremental sparse grid definition uses different
  // point orderings than sgmg/sgmga.  Therefore, the reference grid
  // computations are kept completely separate.

  if (nestedGrid) {
    compute_points_weights(var_sets, type1WeightSets, type2WeightSets);
    if (trackCollocIndices)
      assign_collocation_indices();
  }
  /*
  else {
    // TO DO: hierarchical interpolation must difference interpolants among
    // full point sets rather than evaluating surpluses at point increments
    reference_unique(var_sets); // define reference grid

#ifdef DEBUG
    PCout << "HierarchSparseGridDriver::compute_grid() results:\n"
	  << "uniqueIndexMapping:\n" << uniqueIndexMapping << "\nvar_sets:\n";
    write_data(PCout, var_sets, false, true, true);
    if (trackUniqueProdWeights) {
      PCout << "\ntype1WeightSets:\n"; write_data(PCout, type1WeightSets);
      if (computeType2Weights) {
	PCout << "\ntype2WeightSets:\n";
	write_data(PCout, type2WeightSets, false, true, true);
      }
    }
#endif
  }
  */
}


void HierarchSparseGridDriver::compute_trial_grid(RealMatrix& var_sets)
{
  // track trial sets that have been evaluated (do here since
  // push_trial_set() used for both new trials and restorations)
  computedTrialSets.insert(trial_set());

  // update collocKey and compute trial variable/weight sets
  update_collocation_key();
  if (nestedGrid) {
    if (type1WeightSets.size() <= trialLevel ||
	type2WeightSets.size() <= trialLevel) {
      type1WeightSets.resize(trialLevel+1);
      type2WeightSets.resize(trialLevel+1);
    }
    RealVectorArray& t1_wts_l = type1WeightSets[trialLevel];
    RealMatrixArray& t2_wts_l = type2WeightSets[trialLevel];
    size_t set = t1_wts_l.size();
    RealVector t1_wts_ls; t1_wts_l.push_back(t1_wts_ls); // update in place
    RealMatrix t2_wts_ls; t2_wts_l.push_back(t2_wts_ls); // update in place
    compute_points_weights(var_sets, t1_wts_l[set], t2_wts_l[set]);
    if (trackCollocIndices)
      update_collocation_indices();
  }
  /*
  else {
    compute_points_weights(a2Points, a2Type1Weights, a2Type2Weights);
    // update collocIndices, uniqueIndexMapping, and var_sets,
    // but don't recompute a2 data
    increment_unique(false, true, var_sets);
  }
  */

#ifdef DEBUG
  PCout << "compute_trial_grid() increment:\nunique variable sets:\n"
	<< var_sets;
#endif // DEBUG
}


void HierarchSparseGridDriver::compute_grid_increment(RealMatrix& var_sets)
{
  // update collocKey and compute trial variable/weight sets
  update_smolyak_multi_index();
  update_collocation_key();
  size_t lev, set, num_lev = incrementSets.size(), num_sets;
  if (nestedGrid) {
    if (type1WeightSets.size() < num_lev || type2WeightSets.size() < num_lev)
      { type1WeightSets.resize(num_lev); type2WeightSets.resize(num_lev); }
    // compute total increment evaluations and size var_sets
    size_t num_colloc_pts = 0, start_set, cntr = 0, num_tp_pts;
    for (lev=0; lev<num_lev; ++lev) {
      const UShort3DArray& key_l = collocKey[lev];
      start_set = incrementSets[lev]; num_sets = key_l.size();
      for (set=start_set; set<num_sets; ++set)
	num_colloc_pts += key_l[set].size();
    }
    if (var_sets.numCols() != num_colloc_pts)
      var_sets.shapeUninitialized(numVars, num_colloc_pts);
    // update type1/2 weights and subset view of points
    for (lev=0; lev<num_lev; ++lev) {
      const UShort2DArray& sm_mi_l = smolyakMultiIndex[lev];
      const UShort3DArray&   key_l = collocKey[lev];
      RealVectorArray&    t1_wts_l = type1WeightSets[lev];
      RealMatrixArray&    t2_wts_l = type2WeightSets[lev];
      start_set = incrementSets[lev]; num_sets = sm_mi_l.size();
      for (set=start_set; set<num_sets; ++set) {
	RealVector t1_wts_ls; t1_wts_l.push_back(t1_wts_ls); // update in place
	RealMatrix t2_wts_ls; t2_wts_l.push_back(t2_wts_ls); // update in place
	const UShort2DArray& key_ls = key_l[set];
	num_tp_pts = key_ls.size();
	RealMatrix pts_ls(Teuchos::View, var_sets, numVars, num_tp_pts, 0,cntr);
	compute_points_weights(pts_ls, t1_wts_l[set], t2_wts_l[set],
			       sm_mi_l[set], key_ls);
	cntr += num_tp_pts;
      }
    }
    if (trackCollocIndices)
      update_collocation_indices();
  }
  /*
  else {
    compute_points_weights(a2Points, a2Type1Weights, a2Type2Weights);
    // update collocIndices, uniqueIndexMapping, and var_sets,
    // but don't recompute a2 data
    increment_unique(false, true, var_sets);
  }
  */

#ifdef DEBUG
  PCout << "compute_trial_grid() increment:\nunique variable sets:\n"
	<< var_sets;
#endif // DEBUG
}


void HierarchSparseGridDriver::
compute_points_weights(RealMatrix& pts, RealVector& t1_wts, RealMatrix& t2_wts,
		       const UShortArray& sm_index,
		       const UShort2DArray& colloc_key)
{
  size_t k, l, m, num_tp_pts = colloc_key.size();
  if (pts.numCols() != num_tp_pts)
    pts.shapeUninitialized(numVars, num_tp_pts);
  if (t1_wts.length() != num_tp_pts)
    t1_wts.sizeUninitialized(num_tp_pts);
  if (computeType2Weights && t2_wts.numCols() != num_tp_pts)
    t2_wts.shapeUninitialized(numVars, num_tp_pts);

  // update collocPts1D, type1CollocWts1D, and type2CollocWts1D
  UShortArray total_order;
  level_to_order(sm_index, total_order);
  update_1d_collocation_points_weights(total_order, sm_index);

  // define points and type 1/2 weights; weights are products of 1D weights
  for (k=0; k<num_tp_pts; ++k) {
    const UShortArray& key_k = colloc_key[k];
    Real* pt    =    pts[k]; // column vector
    Real& t1_wt = t1_wts[k]; t1_wt = 1.;
    for (l=0; l<numVars; ++l) {
      pt[l]  =      collocPts1D[sm_index[l]][l][key_k[l]];
      t1_wt *= type1CollocWts1D[sm_index[l]][l][key_k[l]];
    }
    if (computeType2Weights) {
      Real* t2_wt = t2_wts[k]; // column vector
      for (l=0; l<numVars; ++l) {
	Real& t2_wt_l = t2_wt[l]; t2_wt_l = 1.;
	for (m=0; m<numVars; ++m)
	  t2_wt_l *= (m==l) ? type2CollocWts1D[sm_index[m]][m][key_k[m]] :
	                      type1CollocWts1D[sm_index[m]][m][key_k[m]];
      }
    }
  }

#ifdef DEBUG
  PCout << "Tensor product points =\n"; write_data(PCout,pts,false,true,true);
  PCout << "Tensor product weights =\ntype1:\n"; write_data(PCout, t1_wts);
  PCout << "type2:\n"; write_data(PCout, t2_wts, false, true, true);
#endif // DEBUG
}


void HierarchSparseGridDriver::
compute_points_weights(RealMatrix& pts, RealVector& t1_wts, RealMatrix& t2_wts)
{
  compute_points_weights(pts, t1_wts, t2_wts,
			 smolyakMultiIndex[trialLevel].back(),
			 collocKey[trialLevel].back());
}


/** Points are collapsed as required for compute_grid(var_sets), but t1/t2
    weights are hierarchical 2D arrays. */
void HierarchSparseGridDriver::
compute_points_weights(RealMatrix& pts, RealVector2DArray& t1_wts,
		       RealMatrix2DArray& t2_wts)
{
  size_t i, j, cntr = 0, num_colloc_pts = 0, num_tp_pts,
    num_lev = collocKey.size(), num_sets;
  if (t1_wts.size() != num_lev) t1_wts.resize(num_lev);
  if (t2_wts.size() != num_lev) t2_wts.resize(num_lev);
  // define num_colloc_pts
  for (i=0; i<num_lev; ++i) {
    const UShort3DArray& key_i = collocKey[i];
    num_sets = key_i.size();
    if (t1_wts[i].size() != num_sets) t1_wts[i].resize(num_sets);
    if (t2_wts[i].size() != num_sets) t2_wts[i].resize(num_sets);
    for (j=0; j<num_sets; ++j)
      num_colloc_pts += key_i[j].size();
  }
  if (pts.numCols() != num_colloc_pts)
    pts.shapeUninitialized(numVars, num_colloc_pts);

  // define points and type 1/2 weights; weights are products of 1D weights
  for (i=0; i<num_lev; ++i) {
    const UShort3DArray& key_i = collocKey[i];
    num_sets = key_i.size();
    for (j=0; j<num_sets; ++j) {
      const UShort2DArray& key_ij = key_i[j];
      num_tp_pts = key_ij.size();
      // take pts_ij sub-matrix view of full sample matrix pts
      RealMatrix pts_ij(Teuchos::View, pts, numVars, num_tp_pts, 0, cntr);
      compute_points_weights(pts_ij, t1_wts[i][j], t2_wts[i][j],
			     smolyakMultiIndex[i][j], key_ij);
      cntr += num_tp_pts;
    }
  }
}


void HierarchSparseGridDriver::initialize_sets()
{
  // define set O (old) from smolyakMultiIndex and smolyakCoeffs:
  //oldMultiIndex = smolyakMultiIndex;
  oldMultiIndex.clear();
  for (unsigned short lev=0; lev<=ssgLevel; ++lev)
    oldMultiIndex.insert(smolyakMultiIndex[lev].begin(),
			 smolyakMultiIndex[lev].end());

  // computedTrialSets no longer cleared in finalize_sets(), so do on init
  computedTrialSets.clear();

  // compute initial set A (active) by applying add_active_neighbors()
  // to the frontier of smolyakMultiIndex:
  if (dimIsotropic) {
    const UShort2DArray& sm_mi_l = smolyakMultiIndex[ssgLevel];
    size_t i, num_old_sets = sm_mi_l.size();
    for (i=0; i<num_old_sets; ++i)
      add_active_neighbors(sm_mi_l[i], true); // on frontier
  }
  else { // TO DO
    // For anisotropic, need to compute Pareto set.
  }

#ifdef DEBUG
  PCout << "HierarchSparseGridDriver::initialize_sets():\nold sets:\n"
	<< oldMultiIndex << "active sets:\n" << activeMultiIndex << std::endl;
#endif // DEBUG
}


void HierarchSparseGridDriver::push_trial_set(const UShortArray& set)
{
  trialLevel = l1_norm(set);
  if (smolyakMultiIndex.size() <= trialLevel)
    smolyakMultiIndex.resize(trialLevel+1);
  smolyakMultiIndex[trialLevel].push_back(set);

  // collocKey, collocIndices, and uniqueIndexMapping updated within
  // either restore_set() or compute_trial_grid()
}


void HierarchSparseGridDriver::restore_set()
{
  // recompute collocKey from trial set
  update_collocation_key();
  if (nestedGrid) {
    if (trackCollocIndices)
      update_collocation_indices();
    // This approach stores less history than WeightSetsRef approach
    const UShortArray& tr_set = trial_set();
    type1WeightSets[trialLevel].push_back(poppedT1WtSets[tr_set]);
    poppedT1WtSets.erase(tr_set);
    if (computeType2Weights) {
      type2WeightSets[trialLevel].push_back(poppedT2WtSets[tr_set]);
      poppedT2WtSets.erase(tr_set);
    }
  }
  /*
  else { // compute a2
    // update collocIndices and uniqueIndexMapping, but don't update pt/wt sets
    RealMatrix dummy_set;
    increment_unique(true, false, dummy_set);
    merge_unique(); // reset a1 --> INC3
  }
  */
}


void HierarchSparseGridDriver::pop_trial_set()
{
  if (nestedGrid)
    numCollocPts -= collocKey[trialLevel].back().size(); // subtract # trial pts
  /*
  else {
    numCollocPts -= numUnique2; // subtract number of trial points
    uniqueIndexMapping.resize(numCollocPts); // prune trial set from end
  }
  */

  // migrate weights from popped to active status
  const UShortArray& tr_set = trial_set(); // valid prior to smolyakMI pop
  poppedT1WtSets[tr_set] = type1WeightSets[trialLevel].back();
  type1WeightSets[trialLevel].pop_back();
  if (computeType2Weights) {
    poppedT2WtSets[tr_set] = type2WeightSets[trialLevel].back();
    type2WeightSets[trialLevel].pop_back();
  }
  // pop trailing set from smolyakMultiIndex, collocKey, collocIndices
  smolyakMultiIndex[trialLevel].pop_back(); // tr_set no longer valid
  collocKey[trialLevel].pop_back();
  if (trackCollocIndices)
    collocIndices[trialLevel].pop_back();
}


/*
void HierarchSparseGridDriver::merge_set()
{
  if (nestedGrid) {
    // no-op
  }
  //else
  //  merge_unique();
}
*/


void HierarchSparseGridDriver::
finalize_sets(bool output_sets, bool converged_within_tol)
{
  if (output_sets && converged_within_tol) {
    size_t i, j, num_lev = smolyakMultiIndex.size();
    PCout << "Above tolerance index sets:\n";
    for (i=0; i<num_lev; ++i) {
      const UShort2DArray& sm_mi_i = smolyakMultiIndex[i];
      size_t num_sets = sm_mi_i.size();
      if (i==trialLevel) --num_sets; // omit trial set
      for (j=0; j<num_sets; ++j)
	print_index_set(PCout, sm_mi_i[j]);
    }
    PCout << "Below tolerance index sets:\n";
    print_index_set(PCout, smolyakMultiIndex[trialLevel].back());
  }

  // For final answer, push all evaluated sets into old and clear active.
  // Multiple trial insertion approach must be compatible with bookkeeping
  // elsewhere (e.g., Dakota::Approximation), i.e., inc2/inc3 set insertions
  // occur one at a time without mixing.

  // don't insert activeMultiIndex, as this may include sets which have not
  // been evaluated (due to final update_sets() call); use computedTrialSets
  if (nestedGrid) {
    UShortArraySet::iterator it;
    for (it=computedTrialSets.begin(); it!=computedTrialSets.end(); ++it) {
      const UShortArray& tr_set = *it;
      trialLevel = l1_norm(tr_set);
      smolyakMultiIndex[trialLevel].push_back(tr_set);
      update_collocation_key();       // update collocKey
      if (trackCollocIndices)
	update_collocation_indices(); // update collocIndices and numCollocPts
      type1WeightSets[trialLevel].push_back(poppedT1WtSets[tr_set]);
      if (computeType2Weights)
	type2WeightSets[trialLevel].push_back(poppedT2WtSets[tr_set]);
      if (output_sets && converged_within_tol) // print trials below tol
	print_index_set(PCout, tr_set);
    }
  }
  /*
  else {
    // ...
    // update a2 data, uniqueIndexMapping, collocIndices, numCollocPts
    finalize_unique(start_index);// assure no mixing of discrete a2's
    //merge_unique(); // a1 reference update not needed, no addtnl increments
    //update_reference();
  }
  */

  if (output_sets && !converged_within_tol) { // print all together in order
    PCout << "Final index sets:\n";
    size_t i, j, num_lev = smolyakMultiIndex.size();
    for (i=0; i<num_lev; ++i) {
      const UShort2DArray& sm_mi_i = smolyakMultiIndex[i];
      size_t num_sets = sm_mi_i.size();
      for (j=0; j<num_sets; ++j)
	print_index_set(PCout, sm_mi_i[j]);
    }
  }

  activeMultiIndex.clear(); poppedT1WtSets.clear(); poppedT2WtSets.clear();
  // defer since needed for SharedPolyApproxData::finalization_index()
  //computedTrialSets.clear();
}


void HierarchSparseGridDriver::
partition_keys(UShort2DArray& reference_set_range,
	       UShort2DArray& increment_set_range) const
{
  size_t lev, num_lev = smolyakMultiIndex.size(), num_sets;
  reference_set_range.resize(num_lev); increment_set_range.resize(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    UShortArray&  ref_l = reference_set_range[lev];  ref_l.resize(2);
    UShortArray& incr_l = increment_set_range[lev]; incr_l.resize(2);
    const UShort2DArray& sm_mi_l = smolyakMultiIndex[lev];
    num_sets = sm_mi_l.size(); ref_l[0] = 0; incr_l[1] = num_sets;
    if (refineControl == DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
      if (lev == trialLevel) ref_l[1] = incr_l[0] = num_sets - 1;
      else                   ref_l[1] = incr_l[0] = num_sets;
    }
    else
      ref_l[1] = incr_l[0] = incrementSets[lev];
  }
}


void HierarchSparseGridDriver::
partition_keys(UShort3DArray& reference_pt_range,
	       UShort3DArray& increment_pt_range) const
{
  if (refineControl != DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
    PCerr << "Error: point set partitioning only supported in HierarchSparse"
	  << "GridDriver::partition_keys() for generalized sparse grids."
	  << std::endl;
    abort_handler(-1);
  }

  size_t lev, num_lev = collocKey.size(), set, num_sets, num_tp_pts;
  reference_pt_range.resize(num_lev); increment_pt_range.resize(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    num_sets = collocKey[lev].size();
    reference_pt_range[lev].resize(num_sets);
    increment_pt_range[lev].resize(num_sets);
    for (set=0; set<num_sets; ++set) {
      const UShortArray& sm_mi_ls = smolyakMultiIndex[lev][set];
      UShortArray&  ref_ls = reference_pt_range[lev][set];
      UShortArray& incr_ls = increment_pt_range[lev][set];
      ref_ls.resize(2); incr_ls.resize(2);
      num_tp_pts = collocKey[lev][set].size();
      ref_ls[0] = 0; incr_ls[1] = num_tp_pts;
      /*
      if (set == trial_set())
	ref_ls[1] = incr_ls[0] = num_tp_pts-1;
      else
      */
	ref_ls[1] = incr_ls[0] = num_tp_pts;
    }
  }
}

} // namespace Pecos
