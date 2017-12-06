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
  //- Heading: Virtual functions
  //

  /// initializes old/active/evaluation sets for use within the 
  /// generalized sparse grid procedure
  virtual void initialize_sets() = 0;
  /// update smolyakMultiIndex with a new trial set for use within the
  /// generalized sparse grid procedure
  virtual void push_trial_set(const UShortArray& set) = 0;
  /// update collocKey, collocIndices, and uniqueIndexMapping based on
  /// restoration of previous trial to smolyakMultiIndex
  virtual void restore_set();
  /// computes the tensor grid for the index set from push_trial_set()
  virtual void compute_trial_grid(RealMatrix& var_sets);
  /// remove the previously pushed trial set from smolyakMultiIndex
  /// during the course of the generalized sparse grid procedure
  virtual void pop_trial_set() = 0;
  /// merge reference sets with trial set and update reference set
  virtual void merge_set();
  /// accept all remaining trial sets within the generalized sparse
  /// grid procedure
  virtual void finalize_sets(bool output_sets, bool converged_within_tol);

  /// update derived reference data, if required
  virtual void update_reference();

  /// return the trial index set from push_trial_set()
  virtual const UShortArray& trial_set() const = 0;
  /// return the number of unique collocation points in the trial index set
  virtual int unique_trial_points() const;

  /// computes tensor grids for new index sets due to an isotropic/anisotropic
  /// refinement
  virtual void compute_grid_increment(RealMatrix& var_sets);

  /// print smolyakMultiIndex
  virtual void print_smolyak_multi_index() const = 0;

  //
  //- Heading: Member functions
  //

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(unsigned short ssg_level,
    const RealVector& dim_pref, const ShortArray& u_types,
    const ExpansionConfigOptions& ec_options, BasisConfigOptions& bc_options,
    short growth_rate = MODERATE_RESTRICTED_GROWTH);

  /// precompute quadrature rules to the maximum current order for each basis
  /// polynomial (efficiency optimization when rules are expensive to compute)
  void precompute_rules();

  /// initialize collocPts1D and type{1,2}CollocWts1D
  void assign_1d_collocation_points_weights();

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

  /// set ssgLevel
  void level(unsigned short ssg_level);
  /// return ssgLevel
  unsigned short level() const;

  /// set anisoLevelWts
  void dimension_preference(const RealVector& dim_pref);
  /// set anisoLevelWts
  void anisotropic_weights(const RealVector& aniso_wts);
  /// return anisoLevelWts
  const RealVector& anisotropic_weights() const;
  /// return dimIsotropic
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

  /// return activeMultiIndex
  const UShortArraySet& active_multi_index() const;
  /// return computedTrialSets
  const UShortArraySet& computed_trial_sets() const;

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

  /// the Smolyak sparse grid level
  unsigned short ssgLevel;

  /// flag indicating a dimension isotropic grid
  bool dimIsotropic;
  // vector of dimension preference levels for dimension anisotropic grids
  //RealVector dimPref;
  /// weighting vector for dimension anisotropic grids
  RealVector anisoLevelWts;

  /// enumeration for rate of exponential growth in nested rules
  short growthRate;

  // type of expansion refinement
  //short refineType;
  /// algorithm control governing expansion refinement
  short refineControl;

  /// the current number of unique points in the grid
  int numCollocPts;
  /// flag indicating when numCollocPts needs to be recomputed due to an
  /// update to the sparse grid settings
  bool updateGridSize;

  /// old reference index sets for generalized sparse grids
  UShortArraySet oldMultiIndex; // or UShort2DArray
  /// active index sets under current consideration for inclusion in a
  /// generalized sparse grid
  UShortArraySet activeMultiIndex; // or UShort2DArray
  /// subset of active set that have been evaluated as trial sets
  /// (incremented in compute_trial_grid() and decremented in update_sets())
  UShortArraySet computedTrialSets; // or UShort2DArray

private:

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Data
  //

  /// refinement constraints that ensure that level/anisotropic weight updates
  /// contain all previous multi-index sets
  RealVector axisLowerBounds;
};


inline SparseGridDriver::SparseGridDriver():
  IntegrationDriver(BaseConstructor()), ssgLevel(0), dimIsotropic(true),
  growthRate(MODERATE_RESTRICTED_GROWTH), numCollocPts(0), updateGridSize(true),
  refineControl(NO_CONTROL)//refineType(NO_REFINEMENT)
{ }


inline SparseGridDriver::
SparseGridDriver(unsigned short ssg_level, const RealVector& dim_pref,
		 short growth_rate, short refine_control):
  IntegrationDriver(BaseConstructor()), ssgLevel(ssg_level),
  growthRate(growth_rate), numCollocPts(0), updateGridSize(true),
  refineControl(refine_control) //refineType(NO_REFINEMENT)
{
  if (dim_pref.empty())
    dimIsotropic = true;
  else {
    numVars = dim_pref.length(); // unit length option not supported
    dimension_preference(dim_pref);
  }
}


inline SparseGridDriver::~SparseGridDriver()
{ }


inline unsigned short SparseGridDriver::level() const
{ return ssgLevel; }


inline void SparseGridDriver::level(unsigned short ssg_level)
{
  if (ssgLevel != ssg_level)
    { ssgLevel  = ssg_level; updateGridSize = true; }
}


inline const RealVector& SparseGridDriver::anisotropic_weights() const
{ return anisoLevelWts; }


inline bool SparseGridDriver::isotropic() const
{ return dimIsotropic; }


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


inline const UShortArraySet& SparseGridDriver::active_multi_index() const
{ return activeMultiIndex; }


inline const UShortArraySet& SparseGridDriver::computed_trial_sets() const
{ return computedTrialSets; }


inline void SparseGridDriver::update_reference()
{ /* default implementation is no-op */ }


inline void SparseGridDriver::push_trial_set(const UShortArray& set)
{ /* default implementation is no-op */ }


inline void SparseGridDriver::restore_set()
{ /* default implementation is no-op */ }


inline void SparseGridDriver::pop_trial_set()
{ /* default implementation is no-op */ }


inline void SparseGridDriver::merge_set()
{ /* default implementation is no-op */ }


inline void SparseGridDriver::
finalize_sets(bool output_sets, bool converged_within_tol)
{ /* default implementation is no-op */ }


inline void SparseGridDriver::compute_trial_grid(RealMatrix& var_sets)
{ /* default implementation is no-op */ }


inline void SparseGridDriver::compute_grid_increment(RealMatrix& var_sets)
{ /* default implementation is no-op */ }


inline int SparseGridDriver::unique_trial_points() const
{ return 0; /* default implementation */ }


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
