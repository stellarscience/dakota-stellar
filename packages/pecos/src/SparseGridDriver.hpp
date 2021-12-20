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

#ifndef SPARSE_GRID_DRIVER_HPP
#define SPARSE_GRID_DRIVER_HPP

#include "IntegrationDriver.hpp"
#include "ActiveKey.hpp"
#include "sandia_rules.hpp"

namespace Pecos {


/// Derived integration driver class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis, Newton-Cotes, and Gaussian quadrature rules
    within Smolyak sparse grids. */

class SparseGridDriver: public IntegrationDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  SparseGridDriver();
  /// constructor
  SparseGridDriver(unsigned short ssg_level, const RealVector& dim_pref,
		   short growth_rate, short refine_control);
  /// destructor
  ~SparseGridDriver();

  //
  //- Heading: New virtual functions
  //

  /// initializes old/active/evaluation sets for use within the 
  /// generalized sparse grid procedure
  virtual void initialize_sets();
  /// update smolyakMultiIndex with a new trial set for use within the
  /// generalized sparse grid procedure
  virtual void increment_smolyak_multi_index(const UShortArray& set);

  /// determine whether trial set is available for restoration (by push)
  virtual bool push_trial_available(const ActiveKey& key,
				    const UShortArray& tr_set);
  /// determine whether trial set is available for restoration (by push
  virtual bool push_trial_available(const ActiveKey& key);
  /// determine whether trial set is available for restoration (by push
  virtual bool push_trial_available();

  /// determine index of trial set for restoration (by push)
  virtual size_t push_trial_index(const ActiveKey& key,
				  const UShortArray& tr_set);
  /// determine index of trial set for restoration (by push)
  virtual size_t push_trial_index(const ActiveKey& key);
  /// determine index of trial set for restoration (by push)
  virtual size_t push_trial_index();

  /// return pushIndex (cached lookup result in derived Driver classes),
  /// which may be combined (flattened) or hierarchical (level-specific) index
  virtual size_t push_index(const ActiveKey& key) const;
  /// map pushIndex to consistent (flattened) representation
  virtual size_t restore_index(const ActiveKey& key) const;
  /// return consistent (flattened) index
  virtual size_t finalize_index(size_t i, const ActiveKey& key) const;

  /// update collocKey, collocIndices, and uniqueIndexMapping based on
  /// restoration of previous trial to smolyakMultiIndex
  virtual void push_set();
  /// remove the previously pushed trial set from smolyakMultiIndex
  /// during the course of the generalized sparse grid procedure
  virtual void pop_set();
  /// accept all remaining trial sets within the generalized sparse
  /// grid procedure
  virtual void finalize_sets(bool output_sets, bool converged_within_tol,
			     bool reverted);

  /// computes the tensor grid for the trial index set used in
  /// increment_smolyak_multi_index()
  virtual void compute_trial_grid(RealMatrix& var_sets);

  /// computes a grid increment and evaluates the new parameter sets
  virtual void compute_increment(RealMatrix& var_sets);
  /// restores a previously computed grid increment (no new evaluations)
  virtual void push_increment();
  /// removes a previously computed grid increment
  virtual void pop_increment();
  /// merges a grid increment into the reference grid
  virtual void merge_unique();

  /// update derived reference data, if required
  virtual void update_reference();

  /// return the trial index set used in increment_smolyak_multi_index()
  /// that corresponds to key
  virtual const UShortArray& trial_set(const ActiveKey& key) const;
  /// return the trial index set used in increment_smolyak_multi_index()
  /// that corresponds to activeKey
  virtual const UShortArray& trial_set() const;

  /// return the number of unique collocation points in the trial index set
  virtual int unique_trial_points() const;

  /// update smolyakMultiIndex and smolyakCoeffs while adapting grid
  virtual void update_smolyak_arrays();
  /// print smolyakMultiIndex
  virtual void print_smolyak_multi_index() const = 0;

  /// update active iterators based on activeKey
  virtual void update_active_iterators();

  //
  //- Heading: virtual function redefinitions
  //

  void active_key(const ActiveKey& key);
  void clear_inactive();
  void clear_keys();
  void reset();
  void initialize_grid_parameters(const MultivariateDistribution& mv_dist);

  //
  //- Heading: Member functions
  //

  /// return number of collocation points in active grid
  int collocation_points() const;

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
    const MultivariateDistribution& u_dist,
    const ExpansionConfigOptions& ec_options, BasisConfigOptions& bc_options,
    short growth_rate = MODERATE_RESTRICTED_GROWTH);

  /// set flag indicating that grid details have changed and the size
  /// calculation requires updating
  void clear_size();

  /// return pushIndex (cached lookup result in derived Driver classes)
  size_t push_index() const;

  /// precompute quadrature rules to the maximum current order for each basis
  /// polynomial (efficiency optimization when rules are expensive to compute)
  void precompute_rules();

  // initialize collocPts1D and type{1,2}CollocWts1D
  //void assign_1d_collocation_points_weights();
  /// expand and update collocPts1D and type{1,2}CollocWts1D
  void update_1d_collocation_points_weights();
  /// reset collocPts1D and type{1,2}CollocWts1D for variables with a
  /// distribution parameter change
  void reset_1d_collocation_points_weights();

  /// update axisLowerBounds
  void update_axis_lower_bounds();

  /// accept the best of several trial sets and update old/active
  /// within the generalized sparse grid procedure
  void update_sets(const UShortArray& set_star);
  /// update activeMultiIndex from the passed trial set for use within
  /// the generalized sparse grid procedure
  void add_active_neighbors(const UShortArray& set, bool frontier);

  /// converts an array of sparse grid levels to an array of
  /// quadrature orders based on apiIntegrationRules/apiGrowthRules
  void level_to_order(size_t i, unsigned short level,
		      unsigned short& order);
  /// converts an array of sparse grid levels to an array of
  /// quadrature orders based on apiIntegrationRules/apiGrowthRules
  void level_to_order(const UShortArray& levels, UShortArray& orders);

  /// set active ssgLevel
  void level(unsigned short ssg_level);
  /// return active ssgLevel
  unsigned short level() const;

  /// convert dimension preference and set anisoLevelWts
  void dimension_preference(const RealVector& dim_pref);
  /// set anisoLevelWts
  void anisotropic_weights(const RealVector& aniso_wts);
  /// return anisoLevelWts
  const RealVector& anisotropic_weights() const;
  /// indicates isotropic (empty) weights on sparse grid dimension levels
  bool isotropic() const;

  /// set growthRate
  void growth_rate(short growth_rate);
  /// get growthRate
  short growth_rate() const;

  // get refineType
  //short refinement_type()    const;
  /// set refineControl
  void refinement_control(short cntl);
  /// get refineControl
  short refinement_control() const;

  /// return entry from activeMultiIndex corresponding to key
  const UShortArraySet& active_multi_index(const ActiveKey& key) const;
  /// return entry from activeMultiIndex corresponding to activeKey
  const UShortArraySet& active_multi_index() const;

protected:

  //
  //- Heading: Convenience functions
  //

  /// level to order mapping for interpolation with nested Genz-Keister rules
  static int level_to_order_exp_hgk_interp(int level, int growth);
  /// level to order mapping for interpolation with nested closed rules
  static int level_to_order_exp_closed_interp(int level, int growth);
  /// level to order mapping for interpolation with nested open rules
  static int level_to_order_exp_open_interp(int level, int growth);

  /// print an index set
  void print_index_set(std::ostream& s, const UShortArray& mi) const;

  //
  //- Heading: Data
  //

  /// Smolyak sparse grid levels for each active key
  std::map<ActiveKey, unsigned short> ssgLevel;
  /// iterator to the active Smolyak sparse grid level
  std::map<ActiveKey, unsigned short>::iterator ssgLevIter;

  /// weighting vector for dimension anisotropic grids
  std::map<ActiveKey, RealVector> anisoLevelWts;
  /// weighting vector for dimension anisotropic grids
  std::map<ActiveKey, RealVector>::iterator anisoWtsIter;

  /// enumeration for rate of exponential growth in nested rules
  short growthRate;

  // type of expansion refinement
  //short refineType;
  /// algorithm control governing expansion refinement
  short refineControl;

  /// the number of unique points in each grid
  std::map<ActiveKey, int> numCollocPts;
  /// iterator to the number of unique points in the active grid
  std::map<ActiveKey, int>::iterator numPtsIter;

  /// old reference index sets for generalized sparse grids.  Use std::set
  /// for efficient lookups in add_active_neighbors().
  std::map<ActiveKey, UShortArraySet> oldMultiIndex;
  /// active index sets under current consideration for inclusion in a
  /// generalized sparse grid.  Use std::set for ordering of candidates.
  std::map<ActiveKey, UShortArraySet> activeMultiIndex;
  /// subset of active trial sets that have been evaluated but not selected.
  /// Use std::deque to retain append ordering (mirroring SurrogateData).
  std::map<ActiveKey, UShortArrayDeque> poppedTrialSets;

  /// database key indicating the currently active integration configuration.
  /// the key is a multi-index managing multiple modeling dimensions such as
  /// model form, discretization level, etc.
  ActiveKey activeKey;

private:

  //
  //- Heading: Convenience functions
  //

  /// reset collocPts1D[*][i] and type{1,2}CollocWts1D[*][i] for all levels
  void reset_1d_collocation_points_weights(size_t i);
  /// resize arrays: collocPts1D,type1CollocWts1D,type2CollocWts1D
  void resize_1d_collocation_points_weights();

  //
  //- Heading: Data
  //

  /// refinement constraints that ensure that level/anisotropic weight updates
  /// contain all previous multi-index sets
  std::map<ActiveKey, RealVector> axisLowerBounds;
  // iterator to the active set of axis lower bounds
  //std::map<ActiveKey, RealVector>::iterator axisLBndsIter;
};


inline SparseGridDriver::SparseGridDriver():
  IntegrationDriver(BaseConstructor()), numPtsIter(numCollocPts.end()),
  growthRate(MODERATE_RESTRICTED_GROWTH), //ssgLevel(0), numCollocPts(0),
  refineControl(NO_CONTROL)//, refineType(NO_REFINEMENT)
{
  // initial if-check avoids recursive redundancy
  SparseGridDriver::update_active_iterators();
}


inline SparseGridDriver::
SparseGridDriver(unsigned short ssg_level, const RealVector& dim_pref,
		 short growth_rate, short refine_control):
  IntegrationDriver(BaseConstructor()), numPtsIter(numCollocPts.end()),
  growthRate(growth_rate), //ssgLevel(ssg_level), numCollocPts(0),
  refineControl(refine_control)//, refineType(NO_REFINEMENT)
{
  // ssgLevIter init needs to account for incoming ssg_level
  // So long as we only create new keys and avoid modifying existing ones,
  // this deep copy is not needed.
  std::pair<ActiveKey, unsigned short> us_pair(activeKey/*.copy()*/, ssg_level);
  ssgLevIter = ssgLevel.insert(us_pair).first;
  // now update the rest using default assignments
  SparseGridDriver::update_active_iterators();

  if (!dim_pref.empty()) {
    numVars = dim_pref.length(); // unit length option not supported
    dimension_preference(dim_pref);
  }
}


inline SparseGridDriver::~SparseGridDriver()
{ }


inline void SparseGridDriver::active_key(const ActiveKey& key)
{
  if (activeKey != key) {
    activeKey = key; // shared rep; use deep copy when needed downstream
    update_active_iterators();
  }
}


inline void SparseGridDriver::update_active_iterators()
{
  // Test for change (use numPtsIter due to ctor with ssg_level assignment)
  if (numPtsIter != numCollocPts.end() && numPtsIter->first == activeKey)
    return;

  ssgLevIter   = ssgLevel.find(activeKey);
  numPtsIter   = numCollocPts.find(activeKey);
  anisoWtsIter = anisoLevelWts.find(activeKey);
  //axisLBndsIter = axisLowerBounds.find(activeKey);

  /* So long as we only create new keys and avoid modifying existing ones,
     this deep copy is not needed.
  ActiveKey active_copy; // share 1 deep copy of current active key
  if (ssgLevIter == ssgLevel.end() || numPtsIter == numCollocPts.end() ||
      anisoWtsIter == anisoLevelWts.end())
    active_copy = activeKey.copy();
  */

  if (ssgLevIter == ssgLevel.end()) {
    std::pair<ActiveKey, unsigned short> us_pair(activeKey/*active_copy*/, 0);
    ssgLevIter = ssgLevel.insert(us_pair).first;
  }
  if (numPtsIter == numCollocPts.end()) {
    // use special value for grid_size() (instead of 1 pt for ssgLev 0)
    std::pair<ActiveKey, int> ui_pair(activeKey/*active_copy*/, 0);
    numPtsIter = numCollocPts.insert(ui_pair).first;
  }
  if (anisoWtsIter == anisoLevelWts.end()) {
    std::pair<ActiveKey, RealVector>
      urv_pair(activeKey/*active_copy*/, RealVector());
    anisoWtsIter = anisoLevelWts.insert(urv_pair).first;
  }
  //if (axisLBndsIter == axisLowerBounds.end()) {
  //  std::pair<ActiveKey, RealVector>
  //    urv_pair(activeKey/*active_copy*/, RealVector());
  //  axisLBndsIter = axisLowerBounds.insert(urv_pair).first;
  //}
}


inline void SparseGridDriver::clear_keys()
{
  activeKey.clear();

  ssgLevel.clear();          ssgLevIter    =        ssgLevel.end();
  numCollocPts.clear();      numPtsIter    =    numCollocPts.end();
  anisoLevelWts.clear();     anisoWtsIter  =   anisoLevelWts.end();
  axisLowerBounds.clear(); //axisLBndsIter = axisLowerBounds.end();

  oldMultiIndex.clear();  activeMultiIndex.clear();  poppedTrialSets.clear();

  // this database is shared among all keys
  clear_1d_collocation_points_weights();
}


inline int SparseGridDriver::collocation_points() const
{ return numPtsIter->second; }


inline void SparseGridDriver::clear_size()
{ numPtsIter->second = 0; } // special value indicating a grid update is reqd


inline void SparseGridDriver::reset()
{ IntegrationDriver::reset(); clear_size(); }


inline size_t SparseGridDriver::push_index() const
{ return push_index(activeKey); }


inline unsigned short SparseGridDriver::level() const
{ return ssgLevIter->second; }


inline void SparseGridDriver::level(unsigned short ssg_level)
{
  if (ssgLevIter->second != ssg_level) {
    ssgLevIter->second = ssg_level;
    clear_size();
  }
}


inline const RealVector& SparseGridDriver::anisotropic_weights() const
{ return anisoWtsIter->second; }


inline bool SparseGridDriver::isotropic() const
{ return anisoWtsIter->second.empty(); }


//inline short SparseGridDriver::refinement_type() const
//{ return refineType; }


inline void SparseGridDriver::refinement_control(short cntl)
{ refineControl = cntl; }


inline short SparseGridDriver::refinement_control() const
{ return refineControl; }


inline void SparseGridDriver::growth_rate(short growth_rate)
{ growthRate = growth_rate; }


inline short SparseGridDriver::growth_rate() const
{ return growthRate; }


inline const UShortArraySet& SparseGridDriver::
active_multi_index(const ActiveKey& key) const
{
  std::map<ActiveKey, UShortArraySet>::const_iterator cit
    = activeMultiIndex.find(key);
  if (cit == activeMultiIndex.end()) {
    PCerr << "Error: active key not found in SparseGridDriver::"
	  << "active_multi_index()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const UShortArraySet& SparseGridDriver::active_multi_index() const
{ return active_multi_index(activeKey); }


inline void SparseGridDriver::
level_to_order(size_t i, unsigned short level, unsigned short& order)
{
  //int ilevel = level, iorder;
  //webbur::level_growth_to_order(1, &ilevel, &apiIntegrationRules[i],
  //				  &apiGrowthRules[i], &iorder);
  //order = iorder;

  // if INTERPOLATION_MODE, use mappings that synchronize on the number of
  // interpolated points.  For INTEGRATION_MODE or DEFAULT_MODE, use mappings
  // that synchronize on integrand precision.
  switch (collocRules[i]) {
  case GAUSS_PATTERSON:
    order = (driverMode == INTERPOLATION_MODE) ?
      level_to_order_exp_open_interp(level, growthRate) :
      webbur::level_to_order_exp_gp(level, growthRate);          break;
  case GENZ_KEISTER:
    order = (driverMode == INTERPOLATION_MODE) ?
      level_to_order_exp_hgk_interp(level, growthRate) :
      webbur::level_to_order_exp_hgk(level, growthRate);         break;
  case CLENSHAW_CURTIS: case NEWTON_COTES:
    order = (driverMode == INTERPOLATION_MODE) ?
      level_to_order_exp_closed_interp(level, growthRate) :
      webbur::level_to_order_exp_cc(level, growthRate);          break;
  case FEJER2:
    order = (driverMode == INTERPOLATION_MODE) ?
      level_to_order_exp_open_interp(level, growthRate) :
      webbur::level_to_order_exp_f2(level, growthRate);          break;
  case GAUSS_HERMITE: case GAUSS_LEGENDRE: // weakly nested Gaussian
    order = webbur::level_to_order_linear_wn(level, growthRate); break;
  default:                                 // non-nested Gaussian
    order = webbur::level_to_order_linear_nn(level, growthRate); break;
  }
}


inline void SparseGridDriver::
level_to_order(const UShortArray& levels, UShortArray& orders)
{
  size_t i, num_lev = levels.size();
  if (orders.size() != num_lev)
    orders.resize(num_lev);
  for (i=0; i<num_lev; ++i)
    level_to_order(i, levels[i], orders[i]);
}


inline void SparseGridDriver::
print_index_set(std::ostream& s, const UShortArray& mi) const
{
  size_t j, num_mi = mi.size();
  for (j=0; j<num_mi; ++j)
    s << std::setw(5) << mi[j];
  s << '\n';
}

} // namespace Pecos

#endif
