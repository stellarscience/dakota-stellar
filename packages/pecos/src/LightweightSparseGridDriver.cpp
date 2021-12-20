/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 LightweightSparseGridDriver
//- Description: Implementation code for LightweightSparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "LightweightSparseGridDriver.hpp"
#include "SharedPolyApproxData.hpp"

static const char rcsId[]="@(#) $Id: LightweightSparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void LightweightSparseGridDriver::
initialize_grid(size_t num_v, unsigned short ssg_level)
{
  numVars             = num_v;
  ssgLevel[activeKey] = ssg_level;
  // leave trackUniqueProdWeights as false
  // leave active anisoLevelWts as empty (isotropic)

  UShortArray levels(numVars, ssg_level);
  SharedPolyApproxData::total_order_multi_index(levels, smolyakMultiIndex);
  // TO DO: verify that a lower_bound_offset should NOT be provided
  // (we want ALL index sets to maximize granularity in restriction)
}


void LightweightSparseGridDriver::initialize_sets()
{
  // define set O (old) from smolyakMultiIndex
  UShortArraySet& old_mi = oldMultiIndex[activeKey];
  old_mi.clear();
  old_mi.insert(smolyakMultiIndex.begin(), smolyakMultiIndex.end());

  // compute initial set A (active) by applying add_active_neighbors()
  // to the frontier of smolyakMultiIndex:
  activeMultiIndex[activeKey].clear();
  UShortArraySet::const_iterator cit;
  unsigned short ssg_lev = ssgLevIter->second;
  //bool dim_iso = isotropic();
  for (cit=old_mi.begin(); cit!=old_mi.end(); ++cit)
    if ( /*!dim_iso ||*/ l1_norm(*cit) == ssg_lev )
      add_active_neighbors(*cit, true);//, dim_iso);
}


void LightweightSparseGridDriver::prune_sets(const SizetSet& save_tp)
{
  // prune smolyakMultiIndex (as part of a sparse restriction operation)
  size_t save_index, new_index = 0; StSCIter sit;
  for (sit=save_tp.begin(); sit!=save_tp.end(); ++sit, ++new_index) {
    save_index = *sit;
    if (save_index != new_index)
      smolyakMultiIndex[new_index] = smolyakMultiIndex[save_index];
  }
  smolyakMultiIndex.resize(new_index); // prune trailing

  // define oldMultiIndex from smolyakMultiIndex
  UShortArraySet& old_mi = oldMultiIndex[activeKey];
  old_mi.clear();
  old_mi.insert(smolyakMultiIndex.begin(), smolyakMultiIndex.end());

  // redefine set A (activeMultiIndex) based on admissible forward
  // neighbors for all oldMultiIndex terms (allow for gaps)
  activeMultiIndex[activeKey].clear();
  UShortArraySet::const_iterator oit;
  for (oit=old_mi.begin(); oit!=old_mi.end(); ++oit)
    add_active_neighbors(*oit, false); // not exclusively frontier
}

} // namespace Pecos
