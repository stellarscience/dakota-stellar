/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 HierarchSparseGridDriver
//- Description: Wrapper class for C++ code from packages/quadrature/sparse_grid
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef HIERARCH_SPARSE_GRID_DRIVER_HPP
#define HIERARCH_SPARSE_GRID_DRIVER_HPP

#include "SparseGridDriver.hpp"

namespace Pecos {


/// Derived integration driver class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis, Newton-Cotes, and Gaussian quadrature rules
    within Smolyak sparse grids. */

class HierarchSparseGridDriver: public SparseGridDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  HierarchSparseGridDriver();
  /// constructor
  HierarchSparseGridDriver(unsigned short ssg_level,
			   const RealVector& dim_pref = RealVector(),
			   short growth_rate = MODERATE_RESTRICTED_GROWTH,
			   short refine_control = NO_CONTROL);
  /// destructor
  ~HierarchSparseGridDriver();

  //
  //- Heading: Virtual function redefinitions
  //

  void compute_grid();
  void compute_grid(RealMatrix& var_sets);
  int  grid_size();
  void combine_grid();
  void combined_to_active(bool clear_combined);

  void clear_inactive();
  void clear_keys();

  //const ActiveKey& maximal_grid() const;

  void initialize_sets();
  void increment_smolyak_multi_index(const UShortArray& set);
  bool push_trial_available(const ActiveKey& key, const UShortArray& tr_set);
  bool push_trial_available(const ActiveKey& key);
  bool push_trial_available();
  size_t push_trial_index(const ActiveKey& key, const UShortArray& tr_set);
  size_t push_trial_index(const ActiveKey& key);
  size_t push_trial_index();
  size_t push_index(const ActiveKey& key) const;
  size_t restore_index(const ActiveKey& key) const;
  size_t finalize_index(size_t i, const ActiveKey& key) const;
  void push_set();
  void compute_trial_grid(RealMatrix& var_sets);
  void pop_set();
  void finalize_sets(bool output_sets, bool converged_within_tol,
		     bool reverted);

  const UShortArray& trial_set(const ActiveKey& key) const;
  const UShortArray& trial_set() const;
  unsigned short trial_level() const;
  int unique_trial_points() const;

  void compute_increment(RealMatrix& var_sets);
  void push_increment();
  void pop_increment();
  //void merge_unique();

  void print_smolyak_multi_index() const;

  // concatenate type1WeightSets for use in abstract integration functions
  //const RealVector& type1_weight_sets(); // const;
  // concatenate type2WeightSets for use in abstract integration functions
  //const RealMatrix& type2_weight_sets(); // const;

  UShortUShortPair level_to_delta_pair(size_t i, unsigned short lev_i);
  unsigned short level_to_delta_size(size_t i, unsigned short lev_i);
  void level_to_delta_key(size_t i, unsigned short lev_i,
			  UShortArray& delta_key_i);

  /// convert a Smolyak index set into hierarchical quadrature keys
  void levels_to_delta_keys(const UShortArray& levels,
			    UShort2DArray& delta_keys);
  /// convert a Smolyak index set into the sizes of hierarchical
  /// quadrature increments
  void levels_to_delta_sizes(const UShortArray& levels,
			     UShortArray& delta_sizes);

  //
  //- Heading: Member functions
  //

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		       const MultivariateDistribution& u_dist,
		       const ExpansionConfigOptions& ec_options,
		       BasisConfigOptions& bc_options,
		       short growth_rate = MODERATE_RESTRICTED_GROWTH,
		       bool track_colloc_indices = true);

  void assign_collocation_key();
  void assign_collocation_key(const UShort3DArray& sm_mi,
			      UShort4DArray& colloc_key, bool ordered = true);
  void update_collocation_key_from_trial(const UShortArray& trial_set);
  void update_collocation_key_from_trial(const UShortArray& trial_set,
					 const UShort3DArray& sm_mi,
					 UShort4DArray& colloc_key);
  void update_collocation_key_from_increment(UShortArray& incr_sets);
  void update_collocation_key_from_increment(UShortArray& incr_sets,
					     const UShort3DArray& sm_mi,
					     UShort4DArray& colloc_key);

  void assign_collocation_indices();
  void assign_collocation_indices(const UShort4DArray& colloc_key,
				  Sizet3DArray& colloc_indices,
				  int& num_colloc_pts, bool ordered = true);
  void update_collocation_indices_from_trial(const UShortArray& trial_set);
  void update_collocation_indices_from_trial(const UShortArray& trial_set,
					     const UShort4DArray& colloc_key,
					     Sizet3DArray& colloc_indices,
					     int& num_colloc_pts);
  void update_collocation_indices_from_increment(const UShortArray& incr_sets);
  void update_collocation_indices_from_increment(const UShortArray& incr_sets,
					     const UShort4DArray& colloc_key,
					     Sizet3DArray& colloc_indices,
					     int& num_colloc_pts);

  /// update active numCollocPts from active collocKey (contains unique points)
  void update_collocation_points();
  /// update num_colloc_pts from colloc_key (contains unique points)
  void update_collocation_points(const UShort4DArray& colloc_key,
				 int& num_colloc_pts);

  /// return active entry in incrementSets
  const UShortArray& increment_sets() const;

  // return number of sets within popped{T1,T2}WtSets identified by key
  //size_t popped_sets(const ActiveKey& key) const;
  // return number of sets within popped{T1,T2}WtSets identified by activeKey
  //size_t popped_sets() const;

  /// return active entry in smolyakMultiIndex
  const UShort3DArray& smolyak_multi_index() const;
  /// set active entry in smolyakMultiIndex
  void smolyak_multi_index(const UShort3DArray& sm_mi);
  /// return smolyakMultiIndex[key]
  const UShort3DArray& smolyak_multi_index(const ActiveKey& key) const;
  /// return smolyakMultiIndex
  const std::map<ActiveKey, UShort3DArray>& smolyak_multi_index_map() const;

  /// set trackCollocIndices
  void track_collocation_indices(bool track_colloc_indices);
  /// get trackCollocIndices
  bool track_collocation_indices() const;

  /// return active entry in collocKey
  const UShort4DArray& collocation_key() const;
  /// set active entry in collocKey
  void collocation_key(const UShort4DArray& key);
  /// return collocKey[key]
  const UShort4DArray& collocation_key(const ActiveKey& key) const;
  /// return collocKey
  const std::map<ActiveKey, UShort4DArray>& collocation_key_map() const;

  /// return active entry in collocIndices
  const Sizet3DArray& collocation_indices() const;
  /// set active entry in collocIndices
  void collocation_indices(const Sizet3DArray& indices);
  /// return collocIndices[key]
  const Sizet3DArray& collocation_indices(const ActiveKey& key) const;
  /// return collocIndices
  const std::map<ActiveKey, Sizet3DArray>& collocation_indices_map() const;

  /// convert a one-sided increment as in incrementSets (assumes end is defined
  /// by the total number of sets) to a two-sided key (with beginning and end)
  void increment_sets_to_increment_key(const UShortArray& incr_sets,
				       UShort2DArray& incr_key) const;
  /// define portions of the active level-set hierarchy that are in the
  /// current increment
  void partition_increment_key(UShort2DArray& incr_key) const;
  /// define portions of the active level-set hierarchy that are reference sets
  void partition_reference_key(UShort2DArray& ref_key) const;
  /// discriminate portions of the active level-set hierarchy that are
  /// reference sets from those in the current increment
  void partition_keys(UShort2DArray& ref_key, UShort2DArray& incr_key) const;

  /// discriminate portions of the level-set hierarchy that are in the
  /// current increment for each model key
  void partition_increment_key(
    std::map<ActiveKey, UShort2DArray>& incr_key_map) const;
  /// discriminate portions of the level-set hierarchy that are reference sets
  /// for each model key
  void partition_reference_key(
    std::map<ActiveKey, UShort2DArray>& ref_key_map) const;
  /// discriminate portions of the level-set hierarchy that are reference sets
  /// from those in the current increment for each model key
  void partition_keys(std::map<ActiveKey, UShort2DArray>& ref_key_map,
    std::map<ActiveKey, UShort2DArray>& incr_key_map) const;

  /// discriminate portions of the active level-set-point hierarchy that are
  /// in the reference grid from those in the current increment
  void partition_keys(UShort3DArray&  ref_pt_range,
		      UShort3DArray& incr_pt_range) const;

  /// compute points and weights for all levels and sets of the hierarchical
  /// sparse grid indicated by Smolyak multi-index and collocation key
  void compute_points_weights(const UShort3DArray& sm_mi,
			      const UShort4DArray& colloc_key,
			      RealMatrix2DArray& pts, RealVector2DArray& t1_wts,
			      RealMatrix2DArray& t2_wts);
  /// compute points and weights for a trial set
  void compute_points_weights(RealMatrix& pts, RealVector& t1_wts,
			      RealMatrix& t2_wts);

  // overlay all type{1,2}WeightSets and store in active key
  //void combine_weight_sets(const Sizet3DArray& combined_sm_mi_map,
  //			     RealVector2DArray& comb_t1_wts,
  //			     RealMatrix2DArray& comb_t2_wts);

  /// return active variableSets
  const RealMatrix2DArray& hierarchical_variable_sets() const;
  /// set active variableSets
  void hierarchical_variable_sets(const RealMatrix2DArray& rm2);
  /// return active type1WeightSets
  const RealVector2DArray& type1_hierarchical_weight_sets() const;
  /// set active type1WeightSets
  void type1_hierarchical_weight_sets(const RealVector2DArray& t1_wts);
  /// return active type2WeightSets
  const RealMatrix2DArray& type2_hierarchical_weight_sets() const;
  /// set active type2WeightSets
  void type2_hierarchical_weight_sets(const RealMatrix2DArray& t2_wts);

  /// return complete map for variableSets
  const std::map<ActiveKey, RealMatrix2DArray>& variable_sets_map() const;
  /// return complete map for type1WeightSets
  const std::map<ActiveKey, RealVector2DArray>& type1_weight_sets_map() const;
  /// return complete map for type2WeightSets
  const std::map<ActiveKey, RealMatrix2DArray>& type2_weight_sets_map() const;

  /// return combinedSmolyakMultiIndex
  const UShort3DArray& combined_smolyak_multi_index() const;
  /// return combinedSmolyakMultiIndexMap
  const Sizet3DArray& combined_smolyak_multi_index_map() const;
  /// return combinedCollocKey
  const UShort4DArray& combined_collocation_key() const;
  /// return combinedVarSets
  const RealMatrix2DArray& combined_hierarchical_variable_sets() const;
  /// return combinedT1WeightSets
  const RealVector2DArray& combined_type1_hierarchical_weight_sets() const;
  /// return combinedT2WeightSets
  const RealMatrix2DArray& combined_type2_hierarchical_weight_sets() const;

private:

  //
  //- Heading: Convenience functions
  //

  /// update {smolMI,collocKey,collocInd}Iter from activeKey
  void update_active_iterators();

  /// update active smolyakMultiIndex for change in level and/or aniso weights
  void update_smolyak_multi_index(bool clear_sm_mi = false);

  /// moves all data from popped points/weights to active arrays
  void push_popped_points_weights();

  /// kernel routine used for computing points and weights for a tensor grid
  /// corresponding to a single index set
  void compute_points_weights(const UShortArray& sm_index,
			      const UShort2DArray& colloc_key, RealMatrix& pts,
			      RealVector& t1_wts, RealMatrix& t2_wts);

  //
  //- Heading: Data
  //

  /// flag for use of fully nested 1D rules, allowing formulation using
  /// collocation point increments
  bool nestedGrid;

  /// due to the hierarchical structure, collocation indices only need
  /// to be defined in special cases (e.g., generalized sparse grids
  /// for which index sets can appear in different orders).
  bool trackCollocIndices;

  /// interpolation depth by index set by numVars array for identifying
  /// the index to use within the polynomialBasis for a particular variable
  /** The index sets correspond to j (0-based) for use as indices, which
      are offset from the i indices (1-based) normally used in the Smolyak
      expressions.  The indices correspond to levels, one within each
      anisotropic tensor-product integration of a Smolyak recursion. */
  std::map<ActiveKey, UShort3DArray> smolyakMultiIndex;
  /// iterator for active entry within smolyakMultiIndex
  std::map<ActiveKey, UShort3DArray>::iterator smolMIIter;

  /// level of trial evaluation set passed to increment_smolyak_multi_index()
  /// during an index set-based refinement
  std::map<ActiveKey, unsigned short> trialLevel;
  /// iterator for active entry within trialLevel
  std::map<ActiveKey, unsigned short>::iterator trialLevIter;
  /// identifies the trailing index set increments within smolyakMultiIndex
  /// due to an isotropic/anistropic grid refinement
  std::map<ActiveKey, UShortArray> incrementSets;
  /// iterator for active entry within incrementSets
  std::map<ActiveKey, UShortArray>::iterator incrSetsIter;

  /// levels-by-index sets-by-numDeltaPts-by-numVars array for identifying
  /// the 1-D point indices for sets of tensor-product collocation points
  std::map<ActiveKey, UShort4DArray> collocKey;
  /// iterator for active entry within collocKey
  std::map<ActiveKey, UShort4DArray>::iterator collocKeyIter;

  /// levels-by-index sets-by-numTensorProductPts array for linking the
  /// set of tensor products to the unique collocation points evaluated
  std::map<ActiveKey, Sizet3DArray> collocIndices;
  /// iterator for active entry within collocIndices
  std::map<ActiveKey, Sizet3DArray>::iterator collocIndIter;

  /// the set of points in the sparse grid
  std::map<ActiveKey, RealMatrix2DArray> variableSets;
  /// iterator for active entry within variableSets
  std::map<ActiveKey, RealMatrix2DArray>::iterator varSetsIter;
  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the sparse grid
  std::map<ActiveKey, RealVector2DArray> type1WeightSets;
  /// iterator for active entry within type1WeightSets
  std::map<ActiveKey, RealVector2DArray>::iterator t1WtIter;
  /// the set of type2 weights (for integration of gradient interpolants)
  /// for each derivative component and for each point in the sparse grid
  std::map<ActiveKey, RealMatrix2DArray> type2WeightSets;
  /// iterator for active entry within type2WeightSets
  std::map<ActiveKey, RealMatrix2DArray>::iterator t2WtIter;

  /// multi-index that is the result of combining a set of level expansions
  /// (defined from HierarchSparseGridDriver::smolyakMultiIndex)
  UShort3DArray combinedSmolyakMultiIndex;
  /// mapping of terms when aggregating HierarchSparseGridDriver::
  /// smolyakMultiIndex into combinedSmolyakMultiIndex in pre_combine_data()
  Sizet3DArray combinedSmolyakMultiIndexMap;
  /// collocation key that is the result of combining a set of level expansions
  /// (defined from combinedSmolyakMultiIndex)
  UShort4DArray combinedCollocKey;

  /// combination of grid point sets generated by HierarchSparseGridDriver
  /** Could also be managed within SurrogateData, but would require data
      sharing per HierarchInterpPolyApproximation instance. */
  RealMatrix2DArray combinedVarSets;
  /// combination of HierarchSparseGridDriver::type1WeightSets, consistent
  /// with combination of level expansions
  RealVector2DArray combinedT1WeightSets;
  /// combination of HierarchSparseGridDriver::type2WeightSets, consistent
  /// with combination of level expansions
  RealMatrix2DArray combinedT2WeightSets;

  // concatenation of type1WeightSets RealVector2DArray into a RealVector
  //RealVector concatT1WeightSets;
  // concatenation of type2WeightSets RealMatrix2DArray into a RealMatrix
  //RealMatrix concatT2WeightSets;

  /// popped trial sets that were computed but not selected
  std::map<ActiveKey, UShortArrayDequeArray> poppedLevMultiIndex;
  /// hierarchical index into poppedLevMultiIndex[lev] for data to be pushed
  std::map<ActiveKey, size_t> pushIndex;
  /// flattened index for data to be restored
  std::map<ActiveKey, size_t> restoreIndex;
  /// flattened indices for data to be finalized
  std::map<ActiveKey, SizetArray> finalizeIndex;

  /// point sets popped during decrement for later restoration to variableSets
  std::map<ActiveKey, RealMatrixDequeArray> poppedVarSets;
  /// type 1 weight sets popped during decrement for later restoration to
  /// type1WeightSets
  std::map<ActiveKey, RealVectorDequeArray> poppedT1WtSets;
  /// type 2 weight sets popped during decrement for later restoration to
  /// type2WeightSets
  std::map<ActiveKey, RealMatrixDequeArray> poppedT2WtSets;
};


inline HierarchSparseGridDriver::HierarchSparseGridDriver():
  SparseGridDriver(), nestedGrid(true), trackCollocIndices(true),
  smolMIIter(smolyakMultiIndex.end())
{ HierarchSparseGridDriver::update_active_iterators(); }


inline HierarchSparseGridDriver::
HierarchSparseGridDriver(unsigned short ssg_level, const RealVector& dim_pref,
			 short growth_rate, short refine_control):
  SparseGridDriver(ssg_level, dim_pref, growth_rate, refine_control),
  nestedGrid(true), trackCollocIndices(true),
  smolMIIter(smolyakMultiIndex.end())
{ HierarchSparseGridDriver::update_active_iterators(); }


inline HierarchSparseGridDriver::~HierarchSparseGridDriver()
{ }


inline void HierarchSparseGridDriver::update_active_iterators()
{
  // Test for change
  if (smolMIIter != smolyakMultiIndex.end() && smolMIIter->first == activeKey)
    return;

  smolMIIter    = smolyakMultiIndex.find(activeKey);
  trialLevIter  = trialLevel.find(activeKey);
  incrSetsIter  = incrementSets.find(activeKey);
  collocKeyIter = collocKey.find(activeKey);
  collocIndIter = collocIndices.find(activeKey);
  varSetsIter   = variableSets.find(activeKey);
  t1WtIter      = type1WeightSets.find(activeKey);
  t2WtIter      = type2WeightSets.find(activeKey);

  /* So long as we only create new keys and avoid modifying existing ones,
     this deep copy is not needed.
  ActiveKey active_copy; // share 1 deep copy of current active key
  if (smolMIIter    == smolyakMultiIndex.end() ||
      trialLevIter  == trialLevel.end()        ||
      incrSetsIter  == incrementSets.end()     ||
      collocKeyIter == collocKey.end()         ||
      collocIndIter == collocIndices.end()     ||
      varSetsIter   == variableSets.end()      ||
      t1WtIter      == type1WeightSets.end()   ||
      t2WtIter      == type2WeightSets.end())
    active_copy = activeKey.copy();
  */

  if (smolMIIter == smolyakMultiIndex.end()) {
    std::pair<ActiveKey, UShort3DArray>
      u3a_pair(activeKey/*active_copy*/, UShort3DArray());
    smolMIIter = smolyakMultiIndex.insert(u3a_pair).first;
  }
  if (trialLevIter == trialLevel.end()) {
    std::pair<ActiveKey, unsigned short> us_pair(activeKey/*active_copy*/, 0);
    trialLevIter = trialLevel.insert(us_pair).first;
  }
  if (incrSetsIter == incrementSets.end()) {
    std::pair<ActiveKey, UShortArray>
      ua_pair(activeKey/*active_copy*/, UShortArray());
    incrSetsIter = incrementSets.insert(ua_pair).first;
  }
  if (collocKeyIter == collocKey.end()) {
    std::pair<ActiveKey, UShort4DArray>
      u4a_pair(activeKey/*active_copy*/, UShort4DArray());
    collocKeyIter = collocKey.insert(u4a_pair).first;
  }
  if (collocIndIter == collocIndices.end()) {
    std::pair<ActiveKey, Sizet3DArray>
      s3a_pair(activeKey/*active_copy*/, Sizet3DArray());
    collocIndIter = collocIndices.insert(s3a_pair).first;
  }
  if (varSetsIter == variableSets.end()) {
    std::pair<ActiveKey, RealMatrix2DArray>
      rm2_pair(activeKey/*active_copy*/, RealMatrix2DArray());
    varSetsIter = variableSets.insert(rm2_pair).first;
  }
  if (t1WtIter == type1WeightSets.end()) {
    std::pair<ActiveKey, RealVector2DArray>
      rv2_pair(activeKey/*active_copy*/, RealVector2DArray());
    t1WtIter = type1WeightSets.insert(rv2_pair).first;
  }
  if (t2WtIter == type2WeightSets.end()) {
    std::pair<ActiveKey, RealMatrix2DArray>
      rm2_pair(activeKey/*active_copy*/, RealMatrix2DArray());
    t2WtIter = type2WeightSets.insert(rm2_pair).first;
  }

  SparseGridDriver::update_active_iterators();
}


inline void HierarchSparseGridDriver::assign_collocation_key()
{ assign_collocation_key(smolMIIter->second, collocKeyIter->second); }


inline void HierarchSparseGridDriver::
update_collocation_key_from_trial(const UShortArray& trial_set)
{
  update_collocation_key_from_trial(trial_set, smolMIIter->second,
				    collocKeyIter->second);
}


inline void HierarchSparseGridDriver::
update_collocation_key_from_increment(UShortArray& incr_sets)
{
  update_collocation_key_from_increment(incr_sets, smolMIIter->second,
					collocKeyIter->second);
}


inline void HierarchSparseGridDriver::assign_collocation_indices()
{
  assign_collocation_indices(collocKeyIter->second, collocIndIter->second,
			     numPtsIter->second);
}


inline void HierarchSparseGridDriver::
update_collocation_indices_from_trial(const UShortArray& trial_set)
{
  update_collocation_indices_from_trial(trial_set, collocKeyIter->second,
					collocIndIter->second,
					numPtsIter->second);
}


inline void HierarchSparseGridDriver::
update_collocation_indices_from_increment(const UShortArray& incr_sets)
{
  update_collocation_indices_from_increment(incr_sets, collocKeyIter->second,
					    collocIndIter->second,
					    numPtsIter->second);
}


inline void HierarchSparseGridDriver::update_collocation_points()
{ update_collocation_points(collocKeyIter->second, numPtsIter->second); }


inline void HierarchSparseGridDriver::clear_keys()
{
  SparseGridDriver::clear_keys();

  smolyakMultiIndex.clear();  smolMIIter    = smolyakMultiIndex.end();
  trialLevel.clear();         trialLevIter  = trialLevel.end();
  incrementSets.clear();      incrSetsIter  = incrementSets.end();
  collocKey.clear();          collocKeyIter = collocKey.end();
  collocIndices.clear();      collocIndIter = collocIndices.end();
  variableSets.clear();       varSetsIter   = variableSets.end();
  type1WeightSets.clear();    t1WtIter      = type1WeightSets.end();
  type2WeightSets.clear();    t2WtIter      = type2WeightSets.end();

  poppedLevMultiIndex.clear(); poppedT1WtSets.clear(); poppedT2WtSets.clear();
}


inline const UShortArray& HierarchSparseGridDriver::
trial_set(const ActiveKey& key) const
{
  std::map<ActiveKey, UShort3DArray>::const_iterator sm_cit
    = smolyakMultiIndex.find(key);
  std::map<ActiveKey, unsigned short>::const_iterator tl_cit
    = trialLevel.find(key);
  if (sm_cit == smolyakMultiIndex.end() || tl_cit == trialLevel.end()) {
    PCerr << "Error: key not found in HierarchSparseGridDriver::trial_set()"
	  << std::endl;
    abort_handler(-1);
  }
  return sm_cit->second[tl_cit->second].back();
}


inline const UShortArray& HierarchSparseGridDriver::trial_set() const
{ return smolMIIter->second[trialLevIter->second].back(); }


inline unsigned short HierarchSparseGridDriver::trial_level() const
{ return trialLevIter->second; }


inline int HierarchSparseGridDriver::unique_trial_points() const
{ return collocKeyIter->second[trialLevIter->second].back().size(); }


/** identify if newly-pushed trial set exists within stored data sets */
inline bool HierarchSparseGridDriver::
push_trial_available(const ActiveKey& key, const UShortArray& tr_set)
{
  size_t tr_lev = l1_norm(tr_set);
  const UShortArrayDequeArray& pop_mi = poppedLevMultiIndex[key];
  if (pop_mi.size() <= tr_lev)
    return false;
  else {
    const UShortArrayDeque& pop_mi_l = pop_mi[tr_lev];
    return
      (std::find(pop_mi_l.begin(), pop_mi_l.end(), tr_set) != pop_mi_l.end());
  }
}


/** identify if newly-pushed trial set exists within stored data sets */
inline bool HierarchSparseGridDriver::
push_trial_available(const ActiveKey& key)
{
  const UShortArray& tr_set = trial_set(key);
  size_t tr_lev = l1_norm(tr_set);
  const UShortArrayDequeArray& pop_mi = poppedLevMultiIndex[key];
  if (pop_mi.size() <= tr_lev)
    return false;
  else {
    const UShortArrayDeque& pop_mi_l = pop_mi[tr_lev];
    return
      (std::find(pop_mi_l.begin(), pop_mi_l.end(), tr_set) != pop_mi_l.end());
  }
}


/** identify if newly-pushed trial set exists within stored data sets */
inline bool HierarchSparseGridDriver::push_trial_available()
{ return push_trial_available(activeKey, trial_set()); }


/** identify where newly-pushed trial set exists within stored data sets */
inline size_t HierarchSparseGridDriver::
push_trial_index(const ActiveKey& key, const UShortArray& tr_set)
{
  size_t tr_lev = l1_norm(tr_set);
  const UShortArrayDequeArray& pop_mi = poppedLevMultiIndex[key];
  return (pop_mi.size() <= tr_lev) ? _NPOS : find_index(pop_mi[tr_lev], tr_set);
}


/** identify where newly-pushed trial set exists within stored data sets */
inline size_t HierarchSparseGridDriver::push_trial_index(const ActiveKey& key)
{
  const UShortArray& tr_set = trial_set(key);
  size_t tr_lev = l1_norm(tr_set);
  const UShortArrayDequeArray& pop_mi = poppedLevMultiIndex[key];
  return (pop_mi.size() <= tr_lev) ? _NPOS : find_index(pop_mi[tr_lev], tr_set);
}


/** identify where newly-pushed trial set exists within stored data sets */
inline size_t HierarchSparseGridDriver::push_trial_index()
{ return push_trial_index(activeKey, trial_set()); }


inline size_t HierarchSparseGridDriver::push_index(const ActiveKey& key) const
{
  std::map<ActiveKey, size_t>::const_iterator cit = pushIndex.find(key);
  return (cit == pushIndex.end()) ? _NPOS : cit->second;
}


inline size_t HierarchSparseGridDriver::
restore_index(const ActiveKey& key) const
{
  std::map<ActiveKey, size_t>::const_iterator cit = restoreIndex.find(key);
  return (cit == restoreIndex.end()) ? _NPOS : cit->second;
}


inline size_t HierarchSparseGridDriver::
finalize_index(size_t i, const ActiveKey& key) const
{
  std::map<ActiveKey, SizetArray>::const_iterator cit
    = finalizeIndex.find(key);
  return (cit == finalizeIndex.end()) ? _NPOS : cit->second[i];
}


inline const UShortArray& HierarchSparseGridDriver::increment_sets() const
{ return incrSetsIter->second; }


/*
inline size_t HierarchSparseGridDriver::
popped_sets(const ActiveKey& key) const
{
  // Avoid double lookup since T2 cannot currently exist w/o T1
  //return std::max(poppedT1WtSets[key].size(), poppedT2WtSets[key].size());

  std::map<ActiveKey, RealVectorDequeArray>::const_iterator cit
    = poppedT1WtSets.find(key);
  if (cit == poppedT1WtSets.end())
    return 0;
  else {
    const RealVectorDequeArray& pop_t1w = cit->second;
    size_t lev, num_lev = pop_t1w.size(), num_sets = 0;
    for (lev=0; lev<num_lev; ++lev)
      num_sets += pop_t1w[lev].size();
    return num_sets;
  }
}


inline size_t HierarchSparseGridDriver::popped_sets() const
{ return popped_sets(activeKey); }
*/


inline void HierarchSparseGridDriver::print_smolyak_multi_index() const
{
  const UShort3DArray& sm_mi = smolMIIter->second;
  size_t i, j, k, num_lev = sm_mi.size(), cntr = 1;
  for (i=0; i<num_lev; ++i) {
    const UShort2DArray& sm_mi_i = sm_mi[i];
    size_t num_sets = sm_mi_i.size();
    for (j=0; j<num_sets; ++j, ++cntr) {
      PCout << "Smolyak index set " << cntr << ':';
      print_index_set(PCout, sm_mi_i[j]);
    }
  }
}


inline const UShort3DArray& HierarchSparseGridDriver::
smolyak_multi_index() const
{ return smolMIIter->second; }


inline void HierarchSparseGridDriver::
smolyak_multi_index(const UShort3DArray& sm_mi)
{ smolMIIter->second = sm_mi; }


inline const UShort3DArray& HierarchSparseGridDriver::
smolyak_multi_index(const ActiveKey& key) const
{
  std::map<ActiveKey, UShort3DArray>::const_iterator cit
    = smolyakMultiIndex.find(key);
  if (cit == smolyakMultiIndex.end()) {
    PCerr << "Error: key not found in HierarchSparseGridDriver::"
	  << "smolyak_multi_index()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const std::map<ActiveKey, UShort3DArray>& HierarchSparseGridDriver::
smolyak_multi_index_map() const
{ return smolyakMultiIndex; }


inline void HierarchSparseGridDriver::
track_collocation_indices(bool track_colloc_indices)
{ trackCollocIndices = track_colloc_indices; }


inline bool HierarchSparseGridDriver::track_collocation_indices() const
{ return trackCollocIndices; }


inline const UShort4DArray& HierarchSparseGridDriver::collocation_key() const
{ return collocKeyIter->second; }


inline void HierarchSparseGridDriver::collocation_key(const UShort4DArray& key)
{ collocKeyIter->second = key; }


inline const UShort4DArray& HierarchSparseGridDriver::
collocation_key(const ActiveKey& key) const
{
  std::map<ActiveKey, UShort4DArray>::const_iterator cit
    = collocKey.find(key);
  if (cit == collocKey.end()) {
    PCerr << "Error: key not found in HierarchSparseGridDriver::"
	  << "collocation_key()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const std::map<ActiveKey, UShort4DArray>& HierarchSparseGridDriver::
collocation_key_map() const
{ return collocKey; }


inline const Sizet3DArray& HierarchSparseGridDriver::collocation_indices() const
{ return collocIndIter->second; }


inline void HierarchSparseGridDriver::
collocation_indices(const Sizet3DArray& indices)
{ collocIndIter->second = indices; }


inline const Sizet3DArray& HierarchSparseGridDriver::
collocation_indices(const ActiveKey& key) const
{
  std::map<ActiveKey, Sizet3DArray>::const_iterator cit
    = collocIndices.find(key);
  if (cit == collocIndices.end()) {
    PCerr << "Error: key not found in HierarchSparseGridDriver::"
	  << "collocation_indices()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const std::map<ActiveKey, Sizet3DArray>& HierarchSparseGridDriver::
collocation_indices_map() const
{ return collocIndices; }


inline void HierarchSparseGridDriver::
compute_points_weights(RealMatrix& pts, RealVector& t1_wts, RealMatrix& t2_wts)
{
  unsigned short trial_lev = trialLevIter->second;
  compute_points_weights(smolMIIter->second[trial_lev].back(),
			 collocKeyIter->second[trial_lev].back(),
			 pts, t1_wts, t2_wts);
}


inline void HierarchSparseGridDriver::
compute_points_weights(const UShort3DArray& sm_mi,
		       const UShort4DArray& colloc_key, RealMatrix2DArray& pts,
		       RealVector2DArray& t1_wts, RealMatrix2DArray& t2_wts)
{
  // size consolidated weights according to greatest interpolation depth
  size_t lev, num_lev = sm_mi.size(), set, num_sets;
  pts.resize(num_lev);  t1_wts.resize(num_lev);  t2_wts.resize(num_lev);
  for (lev=0; lev<num_lev; ++lev) {
    const UShort3DArray&   key_l = colloc_key[lev];
    const UShort2DArray& sm_mi_l =  sm_mi[lev];   num_sets = sm_mi_l.size();
    RealMatrixArray&       pts_l =    pts[lev];      pts_l.resize(num_sets);
    RealVectorArray&    t1_wts_l = t1_wts[lev];   t1_wts_l.resize(num_sets);
    RealMatrixArray&    t2_wts_l = t2_wts[lev];   t2_wts_l.resize(num_sets);
    for (set=0; set<num_sets; ++set)
      compute_points_weights(sm_mi_l[set], key_l[set], pts_l[set],
			     t1_wts_l[set], t2_wts_l[set]);
  }
}


inline const RealMatrix2DArray& HierarchSparseGridDriver::
hierarchical_variable_sets() const
{ return varSetsIter->second; }


inline void HierarchSparseGridDriver::
hierarchical_variable_sets(const RealMatrix2DArray& rm2)
{ varSetsIter->second = rm2; }


inline const RealVector2DArray& HierarchSparseGridDriver::
type1_hierarchical_weight_sets() const
{ return t1WtIter->second; }


inline void HierarchSparseGridDriver::
type1_hierarchical_weight_sets(const RealVector2DArray& rv2)
{ t1WtIter->second = rv2; }


inline const RealMatrix2DArray& HierarchSparseGridDriver::
type2_hierarchical_weight_sets() const
{ return t2WtIter->second; }


inline void HierarchSparseGridDriver::
type2_hierarchical_weight_sets(const RealMatrix2DArray& rm2)
{ t2WtIter->second = rm2; }


inline const std::map<ActiveKey, RealMatrix2DArray>&
HierarchSparseGridDriver::variable_sets_map() const
{ return variableSets; }


inline const std::map<ActiveKey, RealVector2DArray>&
HierarchSparseGridDriver::type1_weight_sets_map() const
{ return type1WeightSets; }


inline const std::map<ActiveKey, RealMatrix2DArray>&
HierarchSparseGridDriver::type2_weight_sets_map() const
{ return type2WeightSets; }


inline const UShort3DArray& HierarchSparseGridDriver::
combined_smolyak_multi_index() const
{ return combinedSmolyakMultiIndex; }


inline const Sizet3DArray& HierarchSparseGridDriver::
combined_smolyak_multi_index_map() const
{ return combinedSmolyakMultiIndexMap; }


inline const UShort4DArray& HierarchSparseGridDriver::
combined_collocation_key() const
{ return combinedCollocKey; }


inline const RealMatrix2DArray& HierarchSparseGridDriver::
combined_hierarchical_variable_sets() const
{ return combinedVarSets; }


inline const RealVector2DArray& HierarchSparseGridDriver::
combined_type1_hierarchical_weight_sets() const
{ return combinedT1WeightSets; }


inline const RealMatrix2DArray& HierarchSparseGridDriver::
combined_type2_hierarchical_weight_sets() const
{ return combinedT2WeightSets; }


inline void HierarchSparseGridDriver::
levels_to_delta_sizes(const UShortArray& levels, UShortArray& delta_sizes)
{
  size_t i, num_lev = levels.size();
  if (delta_sizes.size() != num_lev)
    delta_sizes.resize(num_lev);
  for (i=0; i<num_lev; ++i)
    delta_sizes[i] = level_to_delta_size(i, levels[i]);
}


inline void HierarchSparseGridDriver::
levels_to_delta_keys(const UShortArray& levels, UShort2DArray& delta_keys)
{
  size_t i, num_lev = levels.size();
  if (delta_keys.size() != num_lev)
    delta_keys.resize(num_lev);
  for (i=0; i<num_lev; ++i)
    level_to_delta_key(i, levels[i], delta_keys[i]);
}

} // namespace Pecos

#endif
