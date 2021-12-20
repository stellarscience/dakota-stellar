/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedPolyApproxData
//- Description:  Base Class for Orthogonal/Interpolation Polynomial
//-               Approximations
//-               
//- Owner:        Mike Eldred

#ifndef SHARED_POLY_APPROX_DATA_HPP
#define SHARED_POLY_APPROX_DATA_HPP

#include "SharedBasisApproxData.hpp"
#include "SparseGridDriver.hpp"

namespace Pecos {


/// Container class for various expansion configuration options

/** The ExpansionConfigOptions class provides a simple container class
    for expansion configuration options related to data modes, verbosity,
    refinement, and VBD controls. */

class ExpansionConfigOptions
{
  //
  //- Heading: Friends
  //

  //friend class PolynomialApproximation;

public:

  /// default constructor
  ExpansionConfigOptions();
  /// constructor
  ExpansionConfigOptions(short exp_soln_approach, short exp_basis_type,
			 short combine_type, short discrep_type,
			 short output_level, bool vbd_flag,
			 unsigned short vbd_order, //short refine_type,
			 short refine_cntl, short refine_metric,
			 short refine_stats, size_t max_refine_iter,
			 size_t max_solver_iter, Real conv_tol,
			 unsigned short sc_limit);
  /// copy constructor
  ExpansionConfigOptions(const ExpansionConfigOptions& ec_options);
  /// destructor
  ~ExpansionConfigOptions();

//private:

  /// identifies the approach taken in compute_coefficients(): QUADRATURE,
  /// CUBATURE, *_SPARSE_GRID, *_REGRESSION, or SAMPLING
  short expCoeffsSolnApproach;

  /// identifies the type of basis for the expansion: DEFAULT_BASIS or
  /// {NODAL,HIERARCHICAL}_INTERPOLANT for SC or
  /// {TENSOR_PRODUCT,TOTAL_ORDER,ADAPTED}_BASIS for PCE regression
  short expBasisType;

  /// type of expansion combination: additive, multiplicative, or both
  short combineType;
  /// type of multilevel discrepancy emulation: none, distinct, recursive
  short discrepReduction;

  /// output verbosity level: {SILENT,QUIET,NORMAL,VERBOSE,DEBUG}_OUTPUT
  short outputLevel;

  /// flag indicated use of variance-based decomposition for computing
  /// Sobol' indices
  bool vbdFlag;
  /// limit for order of interactions computed in variance-based decomposition
  unsigned short vbdOrderLimit;

  // type of refinement: {NO,P,H}_REFINEMENT
  //short refineType;
  /// approach for control of refinement: {NO,UNIFORM,LOCAL_ADAPTIVE}_CONTROL
  /// or DIMENSION_ADAPTIVE_CONTROL_{SOBOL,DECAY,GENERALIZED}
  short refineControl;
  /// metric employed for controlling adaptive refinement: COVARIANCE_METRIC,
  /// MIXED_STATS_METRIC, or LEVEL_STATS_METRIC
  short refineMetric;
  /// type of expansion statistics used in computing refinement metric:
  /// ACTIVE_EXPANSION_STATS or COMBINED_EXPANSION_STATS
  short refineStatsType;

  /// control for limiting the maximum number of refinement iterations
  /// in adapted approximation algorithms
  size_t maxRefineIterations;
  /// control for limiting the maximum number of solver iterations in iterative
  /// solver-based approximation algorithms (e.g., regularized regression)
  size_t maxSolverIterations;
  /// convergence tolerance for adapted or iterated approximation algorithms
  Real convergenceTol;
  /// number of consecutive cycles for which convergence criterion
  /// must be met prior to termination
  unsigned short softConvLimit;
};


inline ExpansionConfigOptions::ExpansionConfigOptions():
  expCoeffsSolnApproach(QUADRATURE), expBasisType(DEFAULT_BASIS),
  combineType(NO_COMBINE), discrepReduction(NO_DISCREPANCY),
  outputLevel(NORMAL_OUTPUT), vbdFlag(false), vbdOrderLimit(0),
  /*refineType(NO_REFINEMENT),*/ refineControl(NO_CONTROL),
  refineMetric(NO_METRIC), refineStatsType(NO_EXPANSION_STATS),
  maxRefineIterations(100), maxSolverIterations(100),
  convergenceTol(1.e-4), softConvLimit(3)
{ }


inline ExpansionConfigOptions::
ExpansionConfigOptions(short exp_soln_approach, short exp_basis_type,
		       short combine_type, short discrep_type,
		       short output_level, bool vbd_flag,
		       unsigned short vbd_order, //short refine_type,
		       short refine_cntl, short refine_metric,
		       short refine_stats, size_t max_refine_iter,
		       size_t max_solver_iter, Real conv_tol,
		       unsigned short sc_limit):
  expCoeffsSolnApproach(exp_soln_approach), expBasisType(exp_basis_type),
  combineType(combine_type), discrepReduction(discrep_type),
  outputLevel(output_level), vbdFlag(vbd_flag), vbdOrderLimit(vbd_order),
  /*refineType(refine_type),*/ refineControl(refine_cntl),
  refineMetric(refine_metric), refineStatsType(refine_stats),
  maxRefineIterations(max_refine_iter), maxSolverIterations(max_solver_iter),
  convergenceTol(conv_tol), softConvLimit(sc_limit)
{ }


inline ExpansionConfigOptions::
ExpansionConfigOptions(const ExpansionConfigOptions& ec_options):
  expCoeffsSolnApproach(ec_options.expCoeffsSolnApproach),
  expBasisType(ec_options.expBasisType), combineType(ec_options.combineType),
  discrepReduction(ec_options.discrepReduction),
  outputLevel(ec_options.outputLevel), vbdFlag(ec_options.vbdFlag),
  vbdOrderLimit(ec_options.vbdOrderLimit),
  //refineType(ec_options.refineType),
  refineControl(ec_options.refineControl),
  refineMetric(ec_options.refineMetric),
  refineStatsType(ec_options.refineStatsType),
  maxRefineIterations(ec_options.maxRefineIterations),
  maxSolverIterations(ec_options.maxSolverIterations),
  convergenceTol(ec_options.convergenceTol),
  softConvLimit(ec_options.softConvLimit)
{ }


inline ExpansionConfigOptions::~ExpansionConfigOptions()
{ }


/// Container class for various basis configuration options

/** The BasisConfigOptions class provides a simple container class
    for basis configuration options related to rule nesting, 
    piecewise basis polynomials, and derivative enhancement. */

class BasisConfigOptions
{
  //
  //- Heading: Friends
  //

  //friend class PolynomialApproximation;

public:

  /// default constructor
  BasisConfigOptions();
  /// constructor
  BasisConfigOptions(bool nested_rules, bool piecewise_basis,
		     bool equidistant_rules, bool use_derivs);
  /// copy constructor
  BasisConfigOptions(const BasisConfigOptions& bc_options);
  /// destructor
  ~BasisConfigOptions();

//private:

  /// flag for use of piecewise basis polynomials
  bool piecewiseBasis;
  /// flag for utilizing derivatives during formation/calculation of expansions
  bool useDerivs;

  /// flag for use of nested integration rules
  bool nestedRules;
  /// flag for use of equidistant points for forming piecewise basis polynomials
  bool equidistantRules;
  /// override interpolation rules to employ Gaussian quadrature
  bool gaussRuleOverride;
  /// override interpolation rules to employ open rules (e.g., Fejer2)
  bool openRuleOverride;
};


inline BasisConfigOptions::BasisConfigOptions():
  piecewiseBasis(false), useDerivs(false), nestedRules(true),
  equidistantRules(true), gaussRuleOverride(true)/*(false)*/, // TO DO: replace
  openRuleOverride(false)
{ }


inline BasisConfigOptions::
BasisConfigOptions(bool nested_rules, bool piecewise_basis,
		   bool equidistant_rules, bool use_derivs):
  piecewiseBasis(piecewise_basis), useDerivs(use_derivs),
  nestedRules(nested_rules), equidistantRules(equidistant_rules),
  gaussRuleOverride(true)/*(false)*/, // TO DO: replace
  openRuleOverride(false)
{ }


inline BasisConfigOptions::
BasisConfigOptions(const BasisConfigOptions& bc_options):
  piecewiseBasis(bc_options.piecewiseBasis), useDerivs(bc_options.useDerivs),
  nestedRules(bc_options.nestedRules),
  equidistantRules(bc_options.equidistantRules),
  gaussRuleOverride(true)/*(false)*/, // TO DO: replace
  openRuleOverride(false)
{ }


inline BasisConfigOptions::~BasisConfigOptions()
{ }


class MultivariateDistribution; // fwd declaration


/// Derived approximation class for global basis polynomials.

/** The SharedPolyApproxData class provides a global approximation
    based on basis polynomials.  This includes orthogonal polynomials
    used for polynomial chaos expansions and interpolation polynomials
    used for stochastic collocation. */

class SharedPolyApproxData: public SharedBasisApproxData
{
  //
  //- Heading: Friends
  //

  friend class PolynomialApproximation;
  friend class InterpPolyApproximation;
  friend class NodalInterpPolyApproximation;
  friend class HierarchInterpPolyApproximation;
  friend class OrthogPolyApproximation;
  friend class ProjectOrthogPolyApproximation;
  friend class RegressOrthogPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  SharedPolyApproxData();
  /// standard constructor
  SharedPolyApproxData(short basis_type, size_t num_vars);
  /// alternate constructor with full configuration options
  SharedPolyApproxData(short basis_type, size_t num_vars,
		       const ExpansionConfigOptions& ec_options,
		       const BasisConfigOptions&     bc_options);
  /// destructorboth
  ~SharedPolyApproxData();

  //
  //- Heading: Virtual member functions
  //

  /// set activeKey
  virtual void active_key(const ActiveKey& key);
  /// remove all keys
  virtual void clear_keys();

  /// allocate the shared data prior to building the set of approximations
  virtual void allocate_data() = 0;

  /// update the shared data prior to rebuilding the set of approximations
  virtual void increment_data();
  /// decrement the previous increment and store its shared data for
  /// later retrieval
  virtual void decrement_data();

  /// test for whether current trial set requires a new approximation
  /// increment or can be restored from a previous trial
  virtual bool push_available();

  /// restores previously popped approximation data
  virtual void pre_push_data();
  /// restores previously popped approximation data
  virtual void post_push_data();

  /// finalizes the shared approximation data following a set of increments
  virtual void pre_finalize_data();
  /// finalizes the shared approximation data following a set of increments
  virtual void post_finalize_data();

  /// tests for whether further advancement in refinement candidates is
  /// available or if refinement opportunities have "saturated"
  virtual bool advancement_available();

  /// clear inactive approximation data
  virtual void clear_inactive_data();

  /// combines all approximation data levels into multilevel approximation data
  virtual void pre_combine_data();
  /// cleans up stored approximation data following combination
  virtual void post_combine_data();
  /// promotes combined approximation data into active data
  virtual void combined_to_active(bool clear_combined = true);

  /// construct polynomial basis, for approximation and/or Gauss points/weights
  virtual void
    construct_basis(const Pecos::MultivariateDistribution& mv_dist) = 0;
  /// update distribution parameters for orthogonal polynomial basis
  virtual void update_basis_distribution_parameters(
    const Pecos::MultivariateDistribution& mv_dist) = 0;

  /// set derived polynomialBasis
  virtual void
    polynomial_basis(const std::vector<Pecos::BasisPolynomial>& poly_basis) = 0;
  /// get derived polynomialBasis
  virtual const std::vector<Pecos::BasisPolynomial>&
    polynomial_basis() const = 0;
  /// get derived polynomialBasis
  virtual std::vector<Pecos::BasisPolynomial>& polynomial_basis() = 0;

  //
  //- Heading: Member functions
  //

  /// return activeKey
  const ActiveKey& active_key() const;

  /// allocate orthogonal polynomial basis types and integration rules
  /// based on u_types and rule options
  static void initialize_orthogonal_basis_types_rules(
    const MultivariateDistribution& u_dist,
    const BasisConfigOptions& bc_options, ShortArray& basis_types,
    ShortArray& colloc_rules);
  /// assign orthogonal polynomial basis type and integration rule based
  /// on u_type and basis configuration options
  static void initialize_orthogonal_basis_type_rule(short u_type,
    const BasisConfigOptions& bc_options, short& basis_type,
    short& colloc_rule);
  /// allocate poly_basis based on basis_types and colloc_rules
  static void initialize_polynomial_basis(const ShortArray& basis_types,
    const ShortArray& colloc_rules, std::vector<BasisPolynomial>& poly_basis);
  /// pass distribution parameters from u_dist to poly_basis
  static void update_basis_distribution_parameters(
    const MultivariateDistribution& u_dist,
    std::vector<BasisPolynomial>& poly_basis);

  // returns index of the data set to be restored from within popped
  // bookkeeping (entry in poppedLevMultiIndex corresponding to key)
  //size_t push_index(const ActiveKey& key, const UShortArray& tr_set);
  /// returns index of the data set to be restored from within popped
  /// bookkeeping (entry in poppedLevMultiIndex corresponding to key)
  size_t push_index(const ActiveKey& key);
  /// returns index of the data set to be restored from within popped
  /// bookkeeping (entry in poppedLevMultiIndex corresponding to activeKey)
  size_t push_index();
  /// maps from push_index (flattened or hierarchical index for internal use)
  /// to a consistent (flattened) index for external use
  size_t restore_index(const ActiveKey& key);
  /// maps from push_index (flattened or hierarchical index for internal use)
  /// to a consistent (flattened) index for external use
  size_t restore_index();
  /// returns index of the i-th data set to be restored from within popped
  /// bookkeeping (entry in poppedLevMultiIndex corresponding to key)
  /// during finalization
  size_t finalize_index(size_t i, const ActiveKey& key);
  /// returns index of the i-th data set to be restored from within popped
  /// bookkeeping (entry in poppedLevMultiIndex corresponding to activeKey)
  /// during finalization
  size_t finalize_index(size_t i);

  /// set expConfigOptions as a group (instead of per option below)
  void configuration_options(const ExpansionConfigOptions& ec_options);
  /// set basisConfigOptions (instead of per option below)
  void configuration_options(const BasisConfigOptions& bc_options);
  /*
  /// set ExpansionConfigOptions::expCoeffsSolnApproach
  void solution_approach(short soln_approach);
  /// get ExpansionConfigOptions::expCoeffsSolnApproach
  short solution_approach() const;

  /// set ExpansionConfigOptions::vbdFlag
  void vbd_flag(bool flag);
  /// get ExpansionConfigOptions::vbdFlag
  bool vbd_flag() const;

  /// set ExpansionConfigOptions::vbdOrderLimit
  void vbd_order_limit(unsigned short vbd_order);
  /// get ExpansionConfigOptions::vbdOrderLimit
  unsigned short vbd_order_limit() const;

  // set ExpansionConfigOptions::refineType
  //void refinement_type(short refine_type);
  // get ExpansionConfigOptions::refineType
  //short refinement_type() const;

  /// set ExpansionConfigOptions::refineControl
  void refinement_control(short refine_cntl);
  /// get ExpansionConfigOptions::refineControl
  short refinement_control() const;

  /// set ExpansionConfigOptions::maxIterations
  void maximum_iterations(size_t max_iter);
  /// get ExpansionConfigOptions::maxIterations
  size_t maximum_iterations() const;

  /// set ExpansionConfigOptions::convergenceTol
  void convergence_tolerance(Real conv_tol);
  /// get ExpansionConfigOptions::convergenceTol
  Real convergence_tolerance() const;
  */
  /// get ExpansionConfigOptions::discrepReduction
  short discrepancy_reduction() const;
  /// set ExpansionConfigOptions::refineStatsType
  void refinement_statistics_type(short stats_type);

  /// return sobolIndexMap
  const BitArrayULongMap& sobol_index_map() const;

  /// set driverRep
  void integration_driver_rep(std::shared_ptr<IntegrationDriver> driver_rep);

  /// set randomVarsKey
  void random_variables_key(const BitArray& random_vars_key);

  /// return the number of expansion terms for a tensor-product
  /// expansion with the provided (anisotropic) quadrature orders
  /// (include_upper_bound = false) or expansion orders (default)
  static size_t tensor_product_terms(const UShortArray& orders,
				     bool include_upper_bound = true);
  /// return the number of expansion terms for a total order expansion of
  /// given (isotropic/anisotropic) expansion orders
  static size_t total_order_terms(const UShortArray& orders,
				  short lower_bound_offset = -1);
  /// return the number of expansion terms for an isotropic total
  /// order expansion of given order
  static size_t total_order_terms(unsigned short order, size_t num_vars,
				  short lower_bound_offset = -1);
  /// return the number of expansion terms for an anisotropic total
  /// order expansion of given max_order and dimension preference
  static size_t total_order_terms(unsigned short max_order,
				  const RealVector& dim_pref,
				  short lower_bound_offset = -1);

  /// initialize expansion multi_index using a tensor-product expansion
  static void tensor_product_multi_index(const UShortArray& orders,
					 UShort2DArray& multi_index,
					 bool include_upper_bound = true);
  /// initialize multi_index using a hierarchical tensor-product expansion
  static void hierarchical_tensor_product_multi_index(
    const UShort2DArray& delta_quad, UShort2DArray& multi_index);
  /// initialize expansion multi_index using an (isotropic/anisotropic)
  /// total-order expansion of given orders 
  static void total_order_multi_index(const UShortArray& orders,
    UShort2DArray& multi_index, short lower_bound_offset = -1,
    size_t max_terms = _NPOS);
  /// initialize expansion multi_index using an isotropic total-order
  /// expansion of given order
  static void total_order_multi_index(unsigned short max_order,
    size_t num_vars, UShort2DArray& multi_index,
    short lower_bound_offset = -1, size_t max_terms = _NPOS);
  /// initialize expansion multi_index using an anisotropic total-order
  /// expansion of given order and dimension preference
  static void total_order_multi_index(unsigned short max_order,
    const RealVector& dim_pref, UShort2DArray& multi_index,
    short lower_bound_offset = -1, size_t max_terms = _NPOS);
  /// define a portion of the multi_index for a total-order expansion
  /// corresponding to a single scalar level
  static void total_order_multi_index_by_level(unsigned short level,
					       size_t num_vars, 
					       UShort2DArray& multi_index);

  /// utility function for incrementing a set of multidimensional indices
  static void increment_indices(UShortArray& indices, const UShortArray& limits,
				bool include_upper_bound);
  /// utility function for incrementing a set of multidimensional terms
  static void increment_terms(UShortArray& terms, size_t& last_index,
			      size_t& prev_index, size_t  term_limit,
			      bool& order_complete);//,bool unique_terms=false);

protected:

  //
  //- Heading: Member functions
  //

  /// return true if matching key values within random variable subset
  bool match_random_key(const UShortArray& key_1, 
			const UShortArray& key_2) const;
  /// return true if matching variable values within nonrandom variable subset
  bool match_nonrandom_vars(const RealVector& vars_1,
			    const RealVector& vars_2) const;

  /// allocate sobolIndices and sobolIndexMap for main effects only
  void allocate_main_sobol();
  // allocate sobolIndices and sobolIndexMap for main and interaction
  // effects including m-way interactions for m <= max_order
  //void allocate_main_interaction_sobol(unsigned short max_order);

  /// Define the sobolIndexMap (which defines the set of Sobol' indices)
  /// from the incoming multi-index.  sobolIndexMap values are initialized
  /// to interaction orders, prior to updating with multi-index increments
  /// in assign_sobol_index_map_values().
  void multi_index_to_sobol_index_map(const UShort2DArray& mi);
  /// return the sobolIndexMap values to interaction orders, prior to updating
  /// with multi-index increments in assign_sobol_index_map_values().
  void reset_sobol_index_map_values();
  /// Define the mapping from sobolIndexMap into sobolIndices
  void assign_sobol_index_map_values();

  //
  //- Heading: Data
  //

  /// pointer to integration driver instance
  std::shared_ptr<IntegrationDriver> driverRep;

  /// an encapsulation of expansion configuration options
  ExpansionConfigOptions expConfigOptions;
  /// an encapsulation of basis configuration options
  BasisConfigOptions basisConfigOptions;

  /// previous quadrature order;
  /// used for tracking need for expansion form updates
  UShortArray quadOrderPrev;
  /// previous Smolyak sparse grid level;
  /// used for tracking need for expansion form updates
  unsigned short ssgLevelPrev;
  /// previous anisotropic weighting;
  /// used for tracking need for expansion form updates
  RealVector anisoWtsPrev;

  /// array of bits identifying the random variable subset within the
  /// active variables (used in all_variables mode)
  BitArray randomVarsKey;
  /// list of indices identifying the random variable subset within the active
  /// variables (used in all_variables mode; defined from randomVarsKey)
  SizetList randomIndices;
  /// list of indices identifying the non-random variable subset within the
  /// active variables (used in all_variables mode; defined from randomVarsKey)
  SizetList nonRandomIndices;

  /// database key indicating the currently active polynomial expansion;
  /// the key is a multi-index managing multiple modeling dimensions such
  /// as model form, discretization level, etc.
  ActiveKey activeKey;

  /// mapping to manage different global sensitivity index options
  /// (e.g. univariate/main effects only vs all effects)
  BitArrayULongMap sobolIndexMap;

private:

  //
  //- Heading: Member functions
  //

  static Real dot(const UShortArray& mi, const RealVector& wts);

  static size_t sum(const UShortArray& mi);

  //
  //- Heading: Data
  //

};


inline SharedPolyApproxData::SharedPolyApproxData():
  ssgLevelPrev(USHRT_MAX)
{ }


inline SharedPolyApproxData::
SharedPolyApproxData(short basis_type, size_t num_vars):
  SharedBasisApproxData(BaseConstructor(), basis_type, num_vars),
  ssgLevelPrev(USHRT_MAX)
{ }


inline SharedPolyApproxData::
SharedPolyApproxData(short basis_type, size_t num_vars,
		     const ExpansionConfigOptions& ec_options,
		     const BasisConfigOptions& bc_options):
  SharedBasisApproxData(BaseConstructor(), basis_type, num_vars),
  expConfigOptions(ec_options), basisConfigOptions(bc_options),
  ssgLevelPrev(USHRT_MAX)
{ }


inline SharedPolyApproxData::~SharedPolyApproxData()
{ }


inline const ActiveKey& SharedPolyApproxData::active_key() const
{ return activeKey; }


inline void SharedPolyApproxData::
configuration_options(const ExpansionConfigOptions& ec_options)
{
  expConfigOptions = ec_options;
  // Note: additional settings that are dependent on configuration options
  // (such as SharedInterpPolyApproxData::barycentricFlag and
  // SharedNodalInterpPolyApproxData::momentInterpType) are updated at run
  // time in allocate_data().
}


inline void SharedPolyApproxData::
configuration_options(const BasisConfigOptions& bc_options)
{
  basisConfigOptions = bc_options;
  // Note: additional settings that are dependent on configuration options
  // (such as SharedInterpPolyApproxData::barycentricFlag and
  // SharedNodalInterpPolyApproxData::momentInterpType) are updated at run
  // time in allocate_data().
}


/*
inline void SharedPolyApproxData::solution_approach(short soln_approach)
{ expConfigOptions.expCoeffsSolnApproach = soln_approach; }


inline short SharedPolyApproxData::solution_approach() const
{ return expConfigOptions.expCoeffsSolnApproach; }


inline void SharedPolyApproxData::vbd_flag(bool flag)
{ expConfigOptions.vbdFlag = flag; }


inline bool SharedPolyApproxData::vbd_flag() const
{ return expConfigOptions.vbdFlag; }


inline void SharedPolyApproxData::vbd_order_limit(unsigned short vbd_order)
{ expConfigOptions.vbdOrderLimit = vbd_order; }


inline unsigned short SharedPolyApproxData::vbd_order_limit() const
{ return expConfigOptions.vbdOrderLimit; }


//inline void SharedPolyApproxData::refinement_type(short refine_type)
//{ expConfigOptions.refineType = refine_type; }


//inline short SharedPolyApproxData::refinement_type() const
//{ return expConfigOptions.refineType; }


inline void SharedPolyApproxData::refinement_control(short refine_cntl)
{ expConfigOptions.refineControl = refine_cntl; }


inline short SharedPolyApproxData::refinement_control() const
{ return expConfigOptions.refineControl; }


inline void SharedPolyApproxData::maximum_iterations(size_t max_iter)
{ expConfigOptions.maxIterations = max_iter; }


inline size_t SharedPolyApproxData::maximum_iterations() const
{ return expConfigOptions.maxIterations; }


inline void SharedPolyApproxData::convergence_tolerance(Real conv_tol)
{ expConfigOptions.convergenceTol = conv_tol; }


inline Real SharedPolyApproxData::convergence_tolerance() const
{ return expConfigOptions.convergenceTol; }
*/


inline short SharedPolyApproxData::discrepancy_reduction() const
{ return expConfigOptions.discrepReduction; }


inline void SharedPolyApproxData::refinement_statistics_type(short stats_type)
{ expConfigOptions.refineStatsType = stats_type; }


inline const BitArrayULongMap& SharedPolyApproxData::sobol_index_map() const
{ return sobolIndexMap; }


inline void SharedPolyApproxData::
integration_driver_rep(std::shared_ptr<IntegrationDriver> driver_rep)
{ driverRep = driver_rep; }


inline void SharedPolyApproxData::
random_variables_key(const BitArray& random_vars_key)
{
  randomVarsKey = random_vars_key;
  randomIndices.clear();
  nonRandomIndices.clear();
  for (size_t i=0; i<numVars; i++)
    if (random_vars_key[i]) randomIndices.push_back(i);
    else                 nonRandomIndices.push_back(i);
}


inline void SharedPolyApproxData::
increment_indices(UShortArray& indices, const UShortArray& limits,
		  bool include_upper_bound)
{
  // perform increment
  size_t n = indices.size(), increment_index = 0;
  ++indices[increment_index];
  // if limit exceeded (including or excluding upper bound value within
  // range of indices), push to next index
  while ( increment_index < n &&
	  ( ( !include_upper_bound && 
	       indices[increment_index] >= limits[increment_index] ) ||
	    (  include_upper_bound && 
	       indices[increment_index] >  limits[increment_index] ) ) ) {
    indices[increment_index] = 0;
    ++increment_index;
    if (increment_index < n)
      ++indices[increment_index];
  }
}


inline void SharedPolyApproxData::
increment_terms(UShortArray& terms, size_t& last_index, size_t& prev_index,
		size_t term_limit, bool& order_complete)//, bool unique_terms)
{
  bool increment_complete = false;
  while (!increment_complete) {
    terms[last_index] = 1; // [1,limit] (not [0,limit-1])
    ++terms[prev_index];
    if (prev_index == 0) {
      increment_complete = true;
      if (terms[prev_index] > term_limit)
	order_complete = true;
    }
    else {
      last_index = prev_index;
      --prev_index;
      if ( //( unique_terms && terms[last_index] <  terms[prev_index]) ||
	   //(!unique_terms &&
	      terms[last_index] <= terms[prev_index]) //)
	increment_complete = true;
    }
  }
}


/*
inline size_t SharedPolyApproxData::
push_index(const ActiveKey& key, const UShortArray& tr_set)
{
  switch (expConfigOptions.refineControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
    SparseGridDriver* sg_driver = (SparseGridDriver*)driverRep;
    return sg_driver->push_trial_index(key, tr_set);// lookup for incoming data
    break;
  }
  default:
    PCerr << "Error: trial set not supported in SharedPolyApproxData::"
	  << "push_index()." << std::endl;
    abort_handler(-1);
    return 0;                                             break;
  }
}
*/


inline size_t SharedPolyApproxData::push_index(const ActiveKey& key)
{
  switch (expConfigOptions.refineControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
    std::shared_ptr<SparseGridDriver> sg_driver =
      std::static_pointer_cast<SparseGridDriver>(driverRep);
    size_t p_index = sg_driver->push_index(key); // retrieve value (no lookup)
    return (p_index == _NPOS) ?
      sg_driver->push_trial_index(key) : // not pre-computed -> perform lookup
      p_index;                           //     pre-computed -> return
    break;
  }
  default: // other refinement types support a single candidate with index 0
    return 0;  break;
  }
}


inline size_t SharedPolyApproxData::push_index()
{
  //return push_index(activeKey); // induces unnecessary lookup(s)

  switch (expConfigOptions.refineControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
    std::shared_ptr<SparseGridDriver> sg_driver =
      std::static_pointer_cast<SparseGridDriver>(driverRep);
    size_t p_index = sg_driver->push_index(); // retrieve pushIndex (no lookup)
    return (p_index == _NPOS) ?
      sg_driver->push_trial_index() : // not pre-computed -> perform lookup
      p_index;                        //     pre-computed -> return
    break;
  }
  default: // other refinement types support a single candidate with index 0
    return 0;  break;
  }
}


inline size_t SharedPolyApproxData::restore_index(const ActiveKey& key)
{
  switch (expConfigOptions.refineControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
    std::shared_ptr<SparseGridDriver> sg_driver =
      std::static_pointer_cast<SparseGridDriver>(driverRep);
    // Option 1: default to be overridden by SharedHierarchPAD:
    //return sg_driver->push_index(key);
    // Option 2: add/use virtual sg_driver->restore_index(key) [more consistent with finalize_index()]
    return sg_driver->restore_index(key);
    break;
  }
  default: // not currently used, but a valid operation
    return 0;  break;
  }
}


inline size_t SharedPolyApproxData::restore_index()
{ return restore_index(activeKey); }


inline size_t SharedPolyApproxData::
finalize_index(size_t i, const ActiveKey& key)
{
  switch (expConfigOptions.refineControl) {
  case DIMENSION_ADAPTIVE_CONTROL_GENERALIZED: {
    std::shared_ptr<SparseGridDriver> sg_driver =
      std::static_pointer_cast<SparseGridDriver>(driverRep);
    return sg_driver->finalize_index(i, key);  break;

    /*
    // Note: popped sets are not explicitly added in computed_trial_sets()
    //       order as in {Incr,Hierarch}SparseGridDriver::finalize_sets().
    //       However, poppedLevMultiIndex et al. become ordered due to
    //       enumeration of ordered active_multi_index().  Rather than
    //       incurring additional overhead by mapping indices, the
    //       compile-time verification below can provide spot checking.

    const UShortArraySet& comp_tr_sets = sg_driver->computed_trial_sets(key);
    UShortArraySet::const_iterator cit = comp_tr_sets.begin();
    std::advance(cit, i); // no operator+ for std::set iterator advancement
    // {Incremental,Hierarch}SparseGridDriver::finalize_sets() updates the
    // grid data with remaining computed trial sets.  This function maps from
    // the order of SparseGridDriver::computedTrialSets to the order of
    // SharedPolyApproxData::poppedLevMultiIndex:
    // > computedTrialSets are in sorted order (std::set)
    // > poppedLevMultiIndex is a std::deque updated by push_back(), but when
    //   generated from increment_sets(), it will reflect the sorted order of
    //   activeMultiIndex (minus one candidate following its selection).
    size_t candidate = push_index(key, *cit);
#ifdef DEBUG
    if (candidate != i) { // activate to test need for this mapping
      PCerr << "Error: SharedPolyApproxData::finalize_index() found index "
	    << "mismatch (" << candidate << ", " << i << ")." << std::endl;
      abort_handler(-1);
    }
#endif // DEBUG
    return candidate;  break;
    */
  }
  default: // not currently used, but a valid operation
    return 0;  break;
  }
}


inline size_t SharedPolyApproxData::finalize_index(size_t i)
{ return finalize_index(i, activeKey); }


inline bool SharedPolyApproxData::
match_random_key(const UShortArray& key_1, const UShortArray& key_2) const
{
  SizetList::const_iterator cit;
  for (cit=randomIndices.begin(); cit!=randomIndices.end(); ++cit)
    if (key_1[*cit] != key_2[*cit])
      return false;
  return true;
}


inline bool SharedPolyApproxData::
match_nonrandom_vars(const RealVector& vars_1, const RealVector& vars_2) const
{
  SizetList::const_iterator cit;
  for (cit=nonRandomIndices.begin(); cit!=nonRandomIndices.end(); ++cit)
    if (vars_1[*cit] != vars_2[*cit]) // double match w/o tolerance
      return false;
  return true;
}


inline size_t SharedPolyApproxData::sum(const UShortArray& mi)
{
  size_t i, num_v = mi.size(), sum = 0;
  for (i=0; i<num_v; ++i)
    sum += mi[i];
  return sum;
}


inline Real SharedPolyApproxData::
dot(const UShortArray& mi, const RealVector& wts)
{
  size_t i, num_v = mi.size();
  Real inner_prod = 0.;
  if (wts.empty())
    for (i=0; i<num_v; ++i)
      inner_prod += mi[i];
  else
    for (i=0; i<num_v; ++i)
      inner_prod += mi[i] * wts[i];
  return inner_prod;
}

} // namespace Pecos

#endif
