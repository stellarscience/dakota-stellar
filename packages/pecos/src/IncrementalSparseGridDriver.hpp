/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 SparseGridDriver
//- Description: Wrapper class for C++ code from packages/quadrature/sparse_grid
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef INCREMENTAL_SPARSE_GRID_DRIVER_HPP
#define INCREMENTAL_SPARSE_GRID_DRIVER_HPP

#include "CombinedSparseGridDriver.hpp"


namespace Pecos {

/// Derived integration driver class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis, Newton-Cotes, and Gaussian quadrature rules
    within Smolyak sparse grids. */

class IncrementalSparseGridDriver: public CombinedSparseGridDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  IncrementalSparseGridDriver();
  /// constructor
  IncrementalSparseGridDriver(unsigned short ssg_level,
			      const RealVector& dim_pref = RealVector(),
			      short growth_rate = MODERATE_RESTRICTED_GROWTH,
			      short refine_control = NO_CONTROL);
  /// destructor
  ~IncrementalSparseGridDriver();

  //
  //- Heading: Virtual function redefinitions
  //

  /// update {smolMI,smolCoeffs,collocKey,collocInd}Iter from activeKey
  void update_active_iterators();

  void compute_grid();
  void compute_grid(RealMatrix& var_sets);
  int grid_size();

  /// initialize all sparse grid settings (distribution params already
  /// set within poly_basis)
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis);

  void clear_inactive();
  void clear_keys();

  void initialize_sets();
  void increment_smolyak_multi_index(const UShortArray& set);
  bool push_trial_available(const ActiveKey& key, const UShortArray& tr_set);
  bool push_trial_available(const ActiveKey& key);
  bool push_trial_available();
  size_t push_trial_index(const ActiveKey& key, const UShortArray& tr_set);
  size_t push_trial_index(const ActiveKey& key);
  size_t push_trial_index();
  size_t push_index(const ActiveKey& key) const;
  //size_t finalize_index(size_t i, const ActiveKey& key) const;

  void compute_trial_grid(RealMatrix& var_sets);
  void push_set();
  void pop_set();
  void finalize_sets(bool output_sets, bool converged_within_tol,
		     bool reverted);

  void compute_increment(RealMatrix& var_sets);
  void push_increment();
  void pop_increment();
  void merge_unique();
  
  /// return smolyakCoeffsRef
  const IntArray& smolyak_coefficients_reference() const;
  /// update smolyakCoeffsRef and type{1,2}WeightSetsRef for use within the
  /// adaptive grid refinement procedures
  void update_reference();

  /// return active trial set under evaluation as a refinement candidate
  const UShortArray& trial_set() const;
  /// return trial set corresponding to key
  const UShortArray& trial_set(const ActiveKey& key) const;

  /// return num_unique2
  int unique_trial_points() const;

  //
  //- Heading: Member functions
  //

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
    const MultivariateDistribution& u_dist,
    const ExpansionConfigOptions& ec_options, BasisConfigOptions& bc_options,
    short growth_rate = MODERATE_RESTRICTED_GROWTH,
    bool track_uniq_prod_wts = true);

  /// update smolyakMultiIndex and smolyakCoeffs
  void update_smolyak_arrays();
  /// update collocKey for the trailing index sets within smolyakMultiIndex
  void update_collocation_key();

  /// define a1{Points,Type1Weights,Type2Weights} based on the reference grid
  void reference_unique(bool update_1d_pts_wts = true);
  /// define a2Points and update collocIndices and uniqueIndexMapping
  /// for trailing index sets within smolyakMultiIndex
  void increment_unique(size_t start_index, bool update_1d_pts_wts = true);
  /// apply all remaining trial sets
  void finalize_unique(size_t start_index);

private:

  //
  //- Heading: Convenience functions
  //

  /// modular helper for public increment_unique(size_t, bool)
  void increment_unique_points_weights(size_t start_index,
    const UShort2DArray& sm_mi, const IntArray& sm_coeffs,
    const IntArray& sm_coeffs_ref, const UShort3DArray& colloc_key,
    Sizet2DArray& colloc_ind, int& num_colloc_pts, RealMatrix& a1_pts,
    RealVector& a1_t1w, RealMatrix& a1_t2w,  RealMatrix& a2_pts,
    RealVector& a2_t1w, RealMatrix& a2_t2w, RealVector& zv, RealVector& r1v,
    RealVector& r2v, IntArray& sind1, BitArray& isu1, IntArray& uind1,
    IntArray& uset1, int& num_u1, IntArray& sind2, BitArray& isu2,
    IntArray& uind2, IntArray& uset2, int& num_u2, IntArray& unique_index_map,
    bool update_1d_pts_wts, RealMatrix& pts, RealVector& t1_wts,
    RealMatrix& t2_wts);
  /// modular helper for public merge_unique()
  void merge_unique_points_weights(const UShort2DArray& sm_mi,
    const IntArray& sm_coeffs, const IntArray& sm_coeffs_ref,
    const UShort3DArray& colloc_key, Sizet2DArray& colloc_ind,
    int& num_colloc_pts, RealMatrix& a1_pts, RealVector& a1_t1w,
    RealMatrix& a1_t2w, RealMatrix& a2_pts, RealVector& a2_t1w,
    RealMatrix& a2_t2w, RealVector& r1v, RealVector& r2v,
    IntArray& sind1, BitArray& isu1, IntArray& uind1, IntArray& uset1,
    int& num_u1, IntArray& sind2, BitArray& isu2, IntArray& uind2,
    IntArray& uset2, int& num_u2, IntArray& unique_index_map,
    RealMatrix& pts, RealVector& t1_wts, RealMatrix& t2_wts);

  /// updates sm_mi from sm_coeffs after uniform/isotropic refinement
  void update_smolyak_arrays(UShort2DArray& sm_mi, IntArray& sm_coeffs);
  /// updates sm_mi from sm_coeffs after anisotropic refinement
  void update_smolyak_arrays_aniso(UShort2DArray& sm_mi, IntArray& sm_coeffs);
  /// increment sm_{mi,coeffs} to sync with new_sm_{mi,coeffs}
  void increment_smolyak_arrays(const UShort2DArray& new_sm_mi,
				const IntArray& new_sm_coeffs,
				UShort2DArray& sm_mi, IntArray& sm_coeffs);
  /// decrement sm_{mi,coeffs} to sync with new_sm_{mi,coeffs}
  void decrement_smolyak_arrays(const UShort2DArray& new_sm_mi,
				const IntArray& new_sm_coeffs,
				UShort2DArray& sm_mi, IntArray& sm_coeffs);

  /// define an increment to the collocation indices
  void update_unique_indices(size_t start_index, int num_uniq1,
			     const IntArray& xdnu1, const IntArray& undx1,
			     const BitArray& isu2,  const IntArray& xdnu2,
			     const IntArray& undx2, IntArray& unique_index_map);

  /// process raw a2 points to create a unique point set increment
  void increment_sparse_points(const Sizet2DArray& colloc_ind,
			       size_t start_index,
			       const BitArray& raw_is_unique,
			       size_t colloc_index_offset,
			       const RealMatrix& raw_incr_pts,
			       RealMatrix& unique_incr_pts);

  /// define the reference type{1,2}WeightSets
  void assign_sparse_weights();
  /// update type{1,2}WeightSets based on a grid increment
  void update_sparse_weights(size_t start_index);
  /// convenience function for updating sparse weights from two sets of
  /// tensor weights and updated coefficients
  void update_sparse_weights(size_t start_index,
			     const UShort3DArray& colloc_key,
			     const Sizet2DArray& colloc_ind, int num_colloc_pts,
			     const IntArray& sm_coeffs,
			     const IntArray& sm_coeffs_ref,
			     const RealVector& a1_t1_wts,
			     const RealMatrix& a1_t2_wts,
			     const RealVector& a2_t1_wts,
			     const RealMatrix& a2_t2_wts,
			     RealVector& unique_t1_wts,
			     RealMatrix& unique_t2_wts);
  /// restore type{1,2}WeightSets to reference values
  void pop_weights();

  //
  //- Heading: Data
  //

  /// reference values for the Smolyak combinatorial coefficients;
  /// used in incremental approaches that update smolyakCoeffs
  std::map<ActiveKey, IntArray> smolyakCoeffsRef;
  /// reference values for the type1 weights corresponding to the current
  /// reference grid; used in incremental approaches that update type1WeightSets
  std::map<ActiveKey, RealVector> type1WeightSetsRef;
  /// reference values for the type2 weights corresponding to the current
  /// reference grid; used in incremental approaches that update type2WeightSets
  std::map<ActiveKey, RealMatrix> type2WeightSetsRef;

  /// index into poppedTrialSets for data to be restored
  std::map<ActiveKey, size_t> pushIndex;
  // indices into poppedTrialSets indicating finalization order
  //std::map<ActiveKey, SizetArray> finalizeIndex;

  /// number of unique points in set 1 (reference)
  std::map<ActiveKey, int> numUnique1;
  /// active entry within numUnique1
  std::map<ActiveKey, int>::iterator numUniq1Iter;
  /// number of unique points in set 2 (increment)
  std::map<ActiveKey, int> numUnique2;
  /// active entry within numUnique2
  std::map<ActiveKey, int>::iterator numUniq2Iter;

  /// random vector used within sgmgg for sorting
  std::map<ActiveKey, RealVector> zVec;
  /// distance values for sorting in set 1 (reference)
  std::map<ActiveKey, RealVector> r1Vec;
  /// distance values for sorting in set 2 (increment)
  std::map<ActiveKey, RealVector> r2Vec;

  /// array of collocation points in set 1 (reference)
  std::map<ActiveKey, RealMatrix> a1Points;
  /// active entry within a1Points
  std::map<ActiveKey, RealMatrix>::iterator a1PIter;
  /// vector of type1 weights in set 1 (reference)
  std::map<ActiveKey, RealVector> a1Type1Weights;
  /// active entry within a1Type1Weights
  std::map<ActiveKey, RealVector>::iterator a1T1WIter;
  /// matrix of type2 weights in set 1 (reference)
  std::map<ActiveKey, RealMatrix> a1Type2Weights;
  /// active entry within a1Type2Weights
  std::map<ActiveKey, RealMatrix>::iterator a1T2WIter;

  /// array of collocation points in set 2 (increment)
  std::map<ActiveKey, RealMatrix> a2Points;
  /// active entry within a2Points
  std::map<ActiveKey, RealMatrix>::iterator a2PIter;
  /// vector of type1 weights in set 2 (increment)
  std::map<ActiveKey, RealVector> a2Type1Weights;
  /// active entry within a2Type1Weights
  std::map<ActiveKey, RealVector>::iterator a2T1WIter;
  /// matrix of type2 weights in set 2 (increment)
  std::map<ActiveKey, RealMatrix> a2Type2Weights;
  /// active entry within a2Type2Weights
  std::map<ActiveKey, RealMatrix>::iterator a2T2WIter;

  /// ascending sort index for set 1 (reference)
  std::map<ActiveKey, IntArray> sortIndex1;
  /// ascending sort index for set 2 (increment)
  std::map<ActiveKey, IntArray> sortIndex2;

  /// index within a1 (reference set) of unique points
  std::map<ActiveKey, IntArray> uniqueSet1;
  /// active entry within uniqueSet1
  std::map<ActiveKey, IntArray>::iterator uniqSet1Iter;
  /// index within a2 (increment set) of unique points
  std::map<ActiveKey, IntArray> uniqueSet2;
  /// active entry within uniqueSet2
  std::map<ActiveKey, IntArray>::iterator uniqSet2Iter;

  /// index within uniqueSet1 corresponding to all of a1
  std::map<ActiveKey, IntArray> uniqueIndex1;
  /// active entry within uniqueIndex1
  std::map<ActiveKey, IntArray>::iterator uniqInd1Iter;
  /// index within uniqueSet2 corresponding to all of a2
  std::map<ActiveKey, IntArray> uniqueIndex2;
  /// active entry within uniqueIndex2
  std::map<ActiveKey, IntArray>::iterator uniqInd2Iter;

  /// key to unique points in set 1 (reference)
  std::map<ActiveKey, BitArray> isUnique1;
  /// active entry within isUnique1
  std::map<ActiveKey, BitArray>::iterator isUniq1Iter;
  /// key to unique points in set 2 (increment)
  std::map<ActiveKey, BitArray> isUnique2;
  /// active entry within isUnique2
  std::map<ActiveKey, BitArray>::iterator isUniq2Iter;
};


inline IncrementalSparseGridDriver::IncrementalSparseGridDriver():
  CombinedSparseGridDriver(), a1PIter(a1Points.end())
{ IncrementalSparseGridDriver::update_active_iterators(); }


inline IncrementalSparseGridDriver::
IncrementalSparseGridDriver(unsigned short ssg_level,
			    const RealVector& dim_pref, short growth_rate,
			    short refine_control):
  CombinedSparseGridDriver(ssg_level, dim_pref, growth_rate, refine_control),
  a1PIter(a1Points.end())
{ IncrementalSparseGridDriver::update_active_iterators(); }


inline IncrementalSparseGridDriver::~IncrementalSparseGridDriver()
{ }


inline void IncrementalSparseGridDriver::update_active_iterators()
{
  // Test for change
  if (a1PIter != a1Points.end() && a1PIter->first == activeKey)
    return;

  a1PIter      = a1Points.find(activeKey);
  a1T1WIter    = a1Type1Weights.find(activeKey);
  a1T2WIter    = a1Type2Weights.find(activeKey);
  a2PIter      = a2Points.find(activeKey);
  a2T1WIter    = a2Type1Weights.find(activeKey);
  a2T2WIter    = a2Type2Weights.find(activeKey);
  numUniq1Iter = numUnique1.find(activeKey);
  numUniq2Iter = numUnique2.find(activeKey);
  uniqSet1Iter = uniqueSet1.find(activeKey);
  uniqSet2Iter = uniqueSet2.find(activeKey);
  uniqInd1Iter = uniqueIndex1.find(activeKey);
  uniqInd2Iter = uniqueIndex2.find(activeKey);
  isUniq1Iter  = isUnique1.find(activeKey);
  isUniq2Iter  = isUnique2.find(activeKey);

  /* So long as we only create new keys and avoid modifying existing ones,
     this deep copy is not needed.
  ActiveKey active_copy; // share 1 deep copy of current active key
  if (a1PIter    == a1Points.end()       || a1T1WIter == a1Type1Weights.end() ||
      a1T2WIter  == a1Type2Weights.end() || a2PIter   == a2Points.end()       ||
      a2T1WIter  == a2Type1Weights.end() || a2T2WIter == a2Type2Weights.end() ||
      numUniq1Iter == numUnique1.end()   || numUniq2Iter == numUnique2.end()  ||
      uniqSet1Iter == uniqueSet1.end()   || uniqSet2Iter == uniqueSet2.end()  ||
      uniqInd1Iter == uniqueIndex1.end() ||
      uniqInd2Iter == uniqueIndex2.end() ||
      isUniq1Iter == isUnique1.end()     || isUniq2Iter  == isUnique2.end())
    active_copy = activeKey.copy();
  */

  if (a1PIter == a1Points.end()) {
    std::pair<ActiveKey, RealMatrix>
      rm_pair(activeKey/*active_copy*/, RealMatrix());
    a1PIter = a1Points.insert(rm_pair).first;
  }
  if (a1T1WIter == a1Type1Weights.end()) {
    std::pair<ActiveKey, RealVector>
      rv_pair(activeKey/*active_copy*/, RealVector());
    a1T1WIter = a1Type1Weights.insert(rv_pair).first;
  }
  if (a1T2WIter == a1Type2Weights.end()) {
    std::pair<ActiveKey, RealMatrix>
      rm_pair(activeKey/*active_copy*/, RealMatrix());
    a1T2WIter = a1Type2Weights.insert(rm_pair).first;
  }
  if (a2PIter == a2Points.end()) {
    std::pair<ActiveKey, RealMatrix>
      rm_pair(activeKey/*active_copy*/, RealMatrix());
    a2PIter = a2Points.insert(rm_pair).first;
  }
  if (a2T1WIter == a2Type1Weights.end()) {
    std::pair<ActiveKey, RealVector>
      rv_pair(activeKey/*active_copy*/, RealVector());
    a2T1WIter = a2Type1Weights.insert(rv_pair).first;
  }
  if (a2T2WIter == a2Type2Weights.end()) {
    std::pair<ActiveKey, RealMatrix>
      rm_pair(activeKey/*active_copy*/, RealMatrix());
    a2T2WIter = a2Type2Weights.insert(rm_pair).first;
  }
  if (numUniq1Iter == numUnique1.end()) {
    std::pair<ActiveKey, int> i_pair(activeKey/*active_copy*/, 0);
    numUniq1Iter = numUnique1.insert(i_pair).first;
  }
  if (numUniq2Iter == numUnique2.end()) {
    std::pair<ActiveKey, int> i_pair(activeKey/*active_copy*/, 0);
    numUniq2Iter = numUnique2.insert(i_pair).first;
  }
  if (uniqSet1Iter == uniqueSet1.end()) {
    std::pair<ActiveKey, IntArray> ia_pair(activeKey/*active_copy*/,IntArray());
    uniqSet1Iter = uniqueSet1.insert(ia_pair).first;
  }
  if (uniqSet2Iter == uniqueSet2.end()) {
    std::pair<ActiveKey, IntArray> ia_pair(activeKey/*active_copy*/,IntArray());
    uniqSet2Iter = uniqueSet2.insert(ia_pair).first;
  }
  if (uniqInd1Iter == uniqueIndex1.end()) {
    std::pair<ActiveKey, IntArray> ia_pair(activeKey/*active_copy*/,IntArray());
    uniqInd1Iter = uniqueIndex1.insert(ia_pair).first;
  }
  if (uniqInd2Iter == uniqueIndex2.end()) {
    std::pair<ActiveKey, IntArray> ia_pair(activeKey/*active_copy*/,IntArray());
    uniqInd2Iter = uniqueIndex2.insert(ia_pair).first;
  }
  if (isUniq1Iter == isUnique1.end()) {
    std::pair<ActiveKey, BitArray> ba_pair(activeKey/*active_copy*/,BitArray());
    isUniq1Iter = isUnique1.insert(ba_pair).first;
  }
  if (isUniq2Iter == isUnique2.end()) {
    std::pair<ActiveKey, BitArray> ba_pair(activeKey/*active_copy*/,BitArray());
    isUniq2Iter = isUnique2.insert(ba_pair).first;
  }

  CombinedSparseGridDriver::update_active_iterators();
}


inline void IncrementalSparseGridDriver::clear_keys()
{
  CombinedSparseGridDriver::clear_keys();

  smolyakCoeffsRef.clear();
  type1WeightSetsRef.clear(); type2WeightSetsRef.clear();

  zVec.clear();            r1Vec.clear();          r2Vec.clear();
  sortIndex1.clear();      sortIndex2.clear();

  numUnique1.clear();      numUniq1Iter = numUnique1.end();
  numUnique2.clear();      numUniq2Iter = numUnique2.end();
  a1Points.clear();        a1PIter      = a1Points.end();
  a1Type1Weights.clear();  a1T1WIter    = a1Type1Weights.end();
  a1Type2Weights.clear();  a1T2WIter    = a1Type2Weights.end();
  a2Points.clear();        a2PIter      = a2Points.end();
  a2Type1Weights.clear();  a2T1WIter    = a2Type1Weights.end();
  a2Type2Weights.clear();  a2T2WIter    = a2Type2Weights.end();
  uniqueSet1.clear();      uniqSet1Iter = uniqueSet1.end();
  uniqueSet2.clear();      uniqSet2Iter = uniqueSet2.end();
  uniqueIndex1.clear();    uniqInd1Iter = uniqueIndex1.end();
  uniqueIndex2.clear();    uniqInd2Iter = uniqueIndex2.end();
  isUnique1.clear();       isUniq1Iter  = isUnique1.end();
  isUnique2.clear();       isUniq2Iter  = isUnique2.end();
}


inline void IncrementalSparseGridDriver::compute_grid(RealMatrix& var_sets)
{
  compute_grid();
  var_sets = varSetsIter->second; // copy
}


inline void IncrementalSparseGridDriver::
reference_unique(bool update_1d_pts_wts)
{
  compute_unique_points_weights(smolMIIter->second, smolCoeffsIter->second,
    collocKeyIter->second, collocIndIter->second, numPtsIter->second,
    a1PIter->second, a1T1WIter->second, a1T2WIter->second, zVec[activeKey],
    r1Vec[activeKey], sortIndex1[activeKey], isUniq1Iter->second,
    uniqInd1Iter->second, uniqSet1Iter->second, numUniq1Iter->second,
    uniqIndMapIter->second, update_1d_pts_wts, varSetsIter->second,
    t1WtIter->second, t2WtIter->second);
}


inline void IncrementalSparseGridDriver::
increment_unique(size_t start_index, bool update_1d_pts_wts)
{
  increment_unique_points_weights(start_index, smolMIIter->second,
    smolCoeffsIter->second, smolyakCoeffsRef[activeKey], collocKeyIter->second,
    collocIndIter->second, numPtsIter->second, a1PIter->second,
    a1T1WIter->second, a1T2WIter->second, a2PIter->second, a2T1WIter->second,
    a2T2WIter->second, zVec[activeKey], r1Vec[activeKey], r2Vec[activeKey],
    sortIndex1[activeKey], isUniq1Iter->second, uniqInd1Iter->second,
    uniqSet1Iter->second, numUniq1Iter->second, sortIndex2[activeKey],
    isUniq2Iter->second,  uniqInd2Iter->second, uniqSet2Iter->second,
    numUniq2Iter->second, uniqIndMapIter->second, update_1d_pts_wts,
    varSetsIter->second,  t1WtIter->second, t2WtIter->second);
}


inline void IncrementalSparseGridDriver::merge_unique()
{
  merge_unique_points_weights(smolMIIter->second, smolCoeffsIter->second,
    smolyakCoeffsRef[activeKey], collocKeyIter->second, collocIndIter->second,
    numPtsIter->second, a1PIter->second, a1T1WIter->second, a1T2WIter->second,
    a2PIter->second, a2T1WIter->second, a2T2WIter->second, r1Vec[activeKey],
    r2Vec[activeKey], sortIndex1[activeKey], isUniq1Iter->second,
    uniqInd1Iter->second, uniqSet1Iter->second, numUniq1Iter->second,
    sortIndex2[activeKey], isUniq2Iter->second, uniqInd2Iter->second,
    uniqSet2Iter->second, numUniq2Iter->second, uniqIndMapIter->second,
    varSetsIter->second, t1WtIter->second, t2WtIter->second);
}


inline void IncrementalSparseGridDriver::finalize_unique(size_t start_index)
{
  increment_unique(start_index, false);
  merge_unique();

  // Note: This doesn't address issue of potential point replication changes
  // between initial trial set status and finalization.  Need an improved
  // mechanism for point restore/finalize in Dakota::Approximation.  Could add
  // a virtual fn to interrogate collocation_indices() from the Approximation
  // level.  Perhaps run some performance tests first to verify that this
  // condition is possible (or does structure of admissible indices prevent
  // replication in trial sets that is not first detected in old sets?).
}


inline const UShortArray& IncrementalSparseGridDriver::
trial_set(const ActiveKey& key) const
{
  std::map<ActiveKey, UShort2DArray>::const_iterator cit
    = smolyakMultiIndex.find(key);
  if (cit == smolyakMultiIndex.end()) {
    PCerr << "Error: key not found in IncrementalSparseGridDriver::trial_set()"
	  << std::endl;
    abort_handler(-1);
  }
  return cit->second.back();
}


inline const UShortArray& IncrementalSparseGridDriver::trial_set() const
{ return smolMIIter->second.back(); } // last set appended to active smolyak MI


/** identify if newly-pushed trial set exists within stored data sets */
inline bool IncrementalSparseGridDriver::
push_trial_available(const ActiveKey& key, const UShortArray& tr_set)
{
  const UShortArrayDeque& pop_tr = poppedTrialSets[key];
  return (std::find(pop_tr.begin(), pop_tr.end(), tr_set) != pop_tr.end());
}


/** identify if newly-pushed trial set exists within stored data sets */
inline bool IncrementalSparseGridDriver::
push_trial_available(const ActiveKey& key)
{
  const UShortArrayDeque& pop_tr = poppedTrialSets[key];
  return
    (std::find(pop_tr.begin(), pop_tr.end(), trial_set(key)) != pop_tr.end());
}


/** identify if newly-pushed trial set exists within stored data sets */
inline bool IncrementalSparseGridDriver::push_trial_available()
{
  const UShortArrayDeque& pop_tr = poppedTrialSets[activeKey];
  return (std::find(pop_tr.begin(), pop_tr.end(), trial_set()) != pop_tr.end());
}


/** identify where newly-pushed trial set exists within stored data sets */
inline size_t IncrementalSparseGridDriver::
push_trial_index(const ActiveKey& key, const UShortArray& tr_set)
{ return find_index(poppedTrialSets[key], tr_set); }


/** identify where newly-pushed trial set exists within stored data sets */
inline size_t IncrementalSparseGridDriver::
push_trial_index(const ActiveKey& key)
{ return find_index(poppedTrialSets[key], trial_set(key)); }


/** identify where newly-pushed trial set exists within stored data sets */
inline size_t IncrementalSparseGridDriver::push_trial_index()
{ return find_index(poppedTrialSets[activeKey], trial_set()); }


inline size_t IncrementalSparseGridDriver::
push_index(const ActiveKey& key) const
{
  std::map<ActiveKey, size_t>::const_iterator cit = pushIndex.find(key);
  return (cit == pushIndex.end()) ? _NPOS : cit->second;
}


/*
inline size_t IncrementalSparseGridDriver::
finalize_index(size_t i, const ActiveKey& key) const
{
  std::map<ActiveKey, SizetArray>::const_iterator cit
    = finalizeIndex.find(key);
  return (cit == finalizeIndex.end()) ? _NPOS : cit->second[i];
}
*/


inline void IncrementalSparseGridDriver::update_reference()
{
  smolyakCoeffsRef[activeKey] = smolCoeffsIter->second;
  if (trackUniqueProdWeights) {
    type1WeightSetsRef[activeKey] = t1WtIter->second;
    if (computeType2Weights)
      type2WeightSetsRef[activeKey] = t2WtIter->second;
  }
}


inline void IncrementalSparseGridDriver::pop_weights()
{
  // restore type{1,2}WeightSets to reference values
  // (update_sparse_weights() involves overlays so nontrivial to back out)
  if (trackUniqueProdWeights) {
    t1WtIter->second = type1WeightSetsRef[activeKey];
    if (computeType2Weights)
      t2WtIter->second = type2WeightSetsRef[activeKey];
  }
}


inline const IntArray& IncrementalSparseGridDriver::
smolyak_coefficients_reference() const
{
  std::map<ActiveKey, IntArray>::const_iterator cit
    = smolyakCoeffsRef.find(activeKey);
  if (cit == smolyakCoeffsRef.end()) {
    PCerr << "Error: active key not found in CombinedSparseGridDriver::"
	  << "smolyak_coefficients_reference()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline int IncrementalSparseGridDriver::unique_trial_points() const
{ return numUniq2Iter->second; }


/** Start from scratch rather than incur incremental coefficient update. */
inline void IncrementalSparseGridDriver::update_smolyak_arrays()
{
  if (isotropic())
    update_smolyak_arrays(smolMIIter->second, smolCoeffsIter->second);
  else
    update_smolyak_arrays_aniso(smolMIIter->second, smolCoeffsIter->second);
}

} // namespace Pecos

#endif
