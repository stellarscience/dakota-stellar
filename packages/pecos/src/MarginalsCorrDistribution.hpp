/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef MARGINALS_CORR_DISTRIBUTION_HPP
#define MARGINALS_CORR_DISTRIBUTION_HPP

#include "MultivariateDistribution.hpp"
#include "RandomVariable.hpp"


namespace Pecos {


/// Class for multivariate distribution based on 1D marginals + correlations

/** The MarginalsCorrDistribution models a multivariate random variable
    distribution using 1D marginals plus a correlation matrix. */

class MarginalsCorrDistribution: public MultivariateDistribution
{
public:

  //
  //- Heading: Constructors and destructor
  //

  MarginalsCorrDistribution();  ///< constructor
  ~MarginalsCorrDistribution(); ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  /// return randomVars
  const std::vector<RandomVariable>& random_variables() const;
  /// return randomVars
  std::vector<RandomVariable>& random_variables();

  /// return randomVars[i]
  const RandomVariable& random_variable(size_t i) const;
  /// return randomVars[i]
  RandomVariable& random_variable(size_t i);

  /// return ranVarTypes
  const ShortArray& random_variable_types() const;
  /// set ranVarTypes
  void random_variable_types(const ShortArray& rv_types);

  /// return ranVarTypes[i]
  short random_variable_type(size_t i) const;
  /// set ranVarTypes[i]
  void random_variable_type(short rv_type, size_t i);

  // BMA TODO: Review why these don't follow typical letter/envelope
  // pattern. Left them for now...

  /// pull non-standardized distribution parameters from pull_mvd
  void pull_distribution_parameters(const MultivariateDistribution& pull_mvd);
  /// pull non-standardized distribution parameters from pull_mvd_rep
  void pull_distribution_parameters
  (const std::shared_ptr<MultivariateDistribution> pull_mvd_rep);

  /// pull non-standardized distribution parameters from pull_mvd, aligning
  /// variables based on label lookups
  void pull_distribution_parameters(const MultivariateDistribution& pull_mvd,
    const StringArray& pull_labels, const StringArray& push_labels);
  /// pull non-standardized distribution parameters from pull_mvd_rep, aligning
  /// variables based on label lookups
  void pull_distribution_parameters
  (const std::shared_ptr<MultivariateDistribution> pull_mvd_rep,
   const StringArray& pull_labels, const StringArray& push_labels);

  /// pull non-standardized distribution parameters for a particular
  /// random variable from pull_mvd
  void pull_distribution_parameters(const MultivariateDistribution& pull_mvd,
				    size_t pull_index, size_t push_index);
  /// pull non-standardized distribution parameters for a particular
  /// random variable from pull_mvd_rep
  void pull_distribution_parameters
  (const std::shared_ptr<MultivariateDistribution> pull_mvd_rep,
   size_t pull_index, size_t push_index);

  /// return activeVars
  const BitArray& active_variables() const;
  /// set activeVars
  void active_variables(const BitArray& );
  /// return activeCorr
  const BitArray& active_correlations() const;
  /// set activeCorr
  void active_correlations(const BitArray& );

  /// return corrMatrix
  const RealSymMatrix& correlation_matrix() const;
  /// set corrMatrix
  void correlation_matrix(const RealSymMatrix& corr);

  /// assemble means and standard deviations from RandomVariable::moments()
  RealRealPairArray moments() const;
  /// assemble means from RandomVariable::mean()
  RealVector means() const;
  /// assemble standard deviations from RandomVariable::standard_deviation()
  RealVector std_deviations() const;
  /// assemble variances from RandomVariable::variance()
  RealVector variances() const;

  /// assemble distribution lower and upper bounds from RandomVariable::bounds()
  RealRealPairArray distribution_bounds() const;
  /// assemble distribution lower bounds from RandomVariable::bounds()
  RealVector distribution_lower_bounds() const;
  /// assemble distribution upper bounds from RandomVariable::bounds()
  RealVector distribution_upper_bounds() const;

  // these bounds are distinct from distribution bounds: for Dakota, they are
  // "global" bounds; in Pecos, they include RangeVariable bounds

  bool global_bounds() const;

  /// check incoming vec for length equal to number of active random variables
  template <typename OrdinalType, typename ScalarType>
  void check_active_length(
    const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& vec,
    const BitArray& mask) const;

  /// set real lower bounds for variable ranges
  void lower_bounds(const RealVector& l_bnds,
		    const BitArray& mask = BitArray());
  /// set integer lower bounds for variable ranges
  void lower_bounds(const IntVector& l_bnds,
		    const BitArray& mask = BitArray());
  /// set real lower bound for a variable range
  void lower_bound(Real l_bnd, size_t rv_index);
  /// set int lower bound for a variable range
  void lower_bound(int  l_bnd, size_t rv_index);

  /// set real upper bounds for variable ranges
  void upper_bounds(const RealVector& u_bnds,
		    const BitArray& mask = BitArray());
  /// set int  upper bounds for variable ranges
  void upper_bounds(const IntVector& u_bnds,
		    const BitArray& mask = BitArray());
  /// set real upper bound for a variable range
  void upper_bound(Real u_bnd, size_t rv_index);
  /// set int  upper bound for a variable range
  void upper_bound(int  u_bnd, size_t rv_index);

  /* This requires clients to bypass the base class and cast a rep pointer
  template <typename OrdinalType, typename ScalarType>
  void lower_bounds(
    const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& l_bnds,
    const BitArray& mask = BitArray());
  template <typename ScalarType>
  void lower_bound(ScalarType l_bnd, size_t rv_index);

  template <typename OrdinalType, typename ScalarType>
  void upper_bounds(
    const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& u_bnds,
    const BitArray& mask = BitArray());
  template <typename ScalarType>
  void upper_bound(ScalarType u_bnd, size_t rv_index);
  */

  //
  //- Heading: Member function definitions
  //

  /// set ranVarTypes and initialize randomVars
  void initialize_types(const ShortArray& rv_types,
			const BitArray& active_vars = BitArray());
  /// assigns corrMatrix and activeCorr; invokes initialize_correlations()
  void initialize_correlations(const RealSymMatrix& corr,
			       const BitArray& active_corr = BitArray());
  /// initializes correlationFlag and performs sanity checks
  void initialize_correlations();

  /// update a scalar distribution parameter within randomVars[v]
  template <typename ValueType>
  void push_parameter(size_t v, short dist_param, ValueType value);
  /// update values for one distribution parameter across a sequence
  /// of random variables
  template <typename OrdinalType, typename ScalarType>
  void push_parameters(size_t start_v, size_t num_v, short dist_param,
    const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values);
  /// update values for one distribution parameter across the set
  /// of random variables with matching RV type
  template <typename OrdinalType, typename ScalarType>
  void push_parameters(short rv_type, short dist_param,
    const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values);
  /// update values for one distribution parameter across a sequence
  /// of random variables
  template <typename ValueType>
  void push_parameters(size_t start_v, size_t num_v, short dist_param,
		       const std::vector<ValueType>& values);
  /// update values for one distribution parameter across the set
  /// of random variables with matching RV type
  template <typename ValueType>
  void push_parameters(short rv_type, short dist_param,
		       const std::vector<ValueType>& values);

  /// update val from a scalar distribution parameter from randomVars[v]
  template <typename ValueType>
  void pull_parameter(size_t v, short dist_param, ValueType& val) const;
  /// define array of values using the identified distribution parameter
  /// across the specified range of random variables
  template <typename ValueType>
  void pull_parameters(size_t start_v, size_t num_v, short dist_param,
		       std::vector<ValueType>& values) const;
  /// define array of values using the identified distribution parameter
  /// across the specified range of random variables
  template <typename ValueType>
  void pull_parameters(short rv_type, short dist_param,
		       std::vector<ValueType>& values) const;
  /// update values for one distribution parameter across the set
  /// of random variables with matching RV type
  template <typename OrdinalType, typename ScalarType>
  void pull_parameters(short rv_type, short dist_param,
    Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values) const;

  /// return a distribution parameter from randomVars[v]
  template <typename ValueType>
  ValueType pull_parameter(size_t v, short dist_param) const;
  /// return array of values for one distribution parameter across the set
  /// of random variables with matching RV type
  template <typename ValueType>
  std::vector<ValueType> pull_parameters(short rv_type, short dist_param) const;

  /// return the size of (non-scalar) distribution data from randomVars[v]
  template <typename ValueType>
  size_t pull_parameter_size(size_t v, short dist_param) const;

  /*
  /// expand corrMatrix from probabilistic variables to combined variables
  void expand_correlation_matrix(size_t num_lead_v, size_t num_prob_v,
				  size_t num_trail_v);
  /// contract corrMatrix from combined variables to probabilistic variables
  void contract_correlation_matrix(size_t num_lead_v, size_t num_prob_v,
				  size_t num_trail_v);
  */

  /// set globalBndsFlag based on ranVarTypes
  void check_global_bounds();
  /// verify that randomVars[i].type() equals rv_type
  void check_random_variable_type(size_t i, short rv_type) const;

  /// draw a sample from the i-th RandomVariable
  template <typename Engine> 
  Real draw_sample(size_t i, Engine& rng) const;
  /// draw a sample from the i-th standardized RandomVariable
  template <typename Engine> 
  Real draw_standard_sample(size_t i, Engine& rng) const;

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  Real pdf(const RealVector& pt) const;
  Real log_pdf(const RealVector& pt) const;

  Real pdf(Real val, size_t rv_index) const;
  Real pdf_gradient(Real val, size_t rv_index) const;
  Real pdf_hessian(Real val, size_t rv_index) const;
  Real log_pdf(Real val, size_t rv_index) const;
  Real log_pdf_gradient(Real val, size_t rv_index) const;
  Real log_pdf_hessian(Real val, size_t rv_index) const;

  /// copy marginals + correlation data between representations
  void copy_rep(std::shared_ptr<MultivariateDistribution> mvd_rep);

  //
  //- Heading: Data
  //

  /// vector of types of each u-space standardized uncertain variable to
  /// which each x-space variable is transformed
  ShortArray ranVarTypes;
  /// vector of random variables encapsulating distribution parameters and
  /// statistical functions (pdf, cdf, etc.)
  std::vector<RandomVariable> randomVars;
  /// subset of randomVars that are currently active (if empty, then
  /// no subset: all variables are active)
  BitArray activeVars;

  /// matrix of random variable correlation coefficients
  RealSymMatrix corrMatrix;
  /// subset of randomVars to which the corrMatrix refers (if empty,
  /// then no subset: correlations are provided for all variables)
  BitArray activeCorr;

  /// precompute the presence of random variables that support global bounds,
  /// for fto accelerate run-time operations
  bool globalBndsFlag;

private:

  //
  //- Heading: Data
  //

};


inline MarginalsCorrDistribution::MarginalsCorrDistribution():
  MultivariateDistribution(BaseConstructor()), globalBndsFlag(false)
{ }


inline MarginalsCorrDistribution::~MarginalsCorrDistribution()
{ }


inline const std::vector<RandomVariable>& MarginalsCorrDistribution::
random_variables() const
{ return randomVars; }


inline std::vector<RandomVariable>& MarginalsCorrDistribution::
random_variables()
{ return randomVars; }


inline const RandomVariable& MarginalsCorrDistribution::
random_variable(size_t rv_index) const
{ return randomVars[rv_index]; }


inline RandomVariable& MarginalsCorrDistribution::
random_variable(size_t rv_index)
{ return randomVars[rv_index]; }


inline void MarginalsCorrDistribution::check_global_bounds()
{
  // As a first cut, identify range variables (design, state) as having
  // "global" bounds, which need to synchronize from global parameter space
  // updates, separate from distribution parameter updates
  globalBndsFlag = false;
  size_t i, num_rv = ranVarTypes.size();
  for (i=0; i<num_rv; ++i)
    if (ranVarTypes[i] == CONTINUOUS_RANGE || ranVarTypes[i] == DISCRETE_RANGE)
      { globalBndsFlag = true; break; }
}


inline void MarginalsCorrDistribution::
check_random_variable_type(size_t rv_index, short rv_type) const
{
  if (randomVars[rv_index].type() != rv_type) {
    PCerr << "Error: inconsistent random variable type in MarginalsCorr"
	  << "Distribution::check_random_variable_type()." << std::endl;
    abort_handler(-1);
  }
}


inline const ShortArray& MarginalsCorrDistribution::
random_variable_types() const
{ return ranVarTypes; }


inline void MarginalsCorrDistribution::
random_variable_types(const ShortArray& rv_types)
{
  ranVarTypes = rv_types;
  check_global_bounds();
}


inline short MarginalsCorrDistribution::
random_variable_type(size_t rv_index) const
{ return ranVarTypes[rv_index]; }


inline void MarginalsCorrDistribution::
random_variable_type(short rv_type, size_t rv_index)
{
  bool new_rv_global = (rv_type == CONTINUOUS_RANGE ||
			rv_type ==   DISCRETE_RANGE);
  if (globalBndsFlag) {
    short& curr_rv_type   = ranVarTypes[rv_index];
    bool   curr_rv_global = (curr_rv_type == CONTINUOUS_RANGE ||
			     curr_rv_type ==   DISCRETE_RANGE);
    curr_rv_type = rv_type;
    if (curr_rv_global && !new_rv_global) // previous global type removed
      check_global_bounds();              // --> check for flag update
  }
  else { // no current types define global bnds --> update globalBndsFlag
         // to true if incoming rv_type defines global bnds
    ranVarTypes[rv_index] = rv_type;
    globalBndsFlag = new_rv_global;
  }
}


inline const BitArray& MarginalsCorrDistribution::active_variables() const
{ return activeVars; }


inline void MarginalsCorrDistribution::
active_variables(const BitArray& active_vars)
{ activeVars = active_vars; }


inline const BitArray& MarginalsCorrDistribution::active_correlations() const
{ return activeCorr; }


inline void MarginalsCorrDistribution::
active_correlations(const BitArray& active_corr)
{ activeCorr = active_corr; }


inline const RealSymMatrix& MarginalsCorrDistribution::
correlation_matrix() const
{ return corrMatrix; }


inline void MarginalsCorrDistribution::
correlation_matrix(const RealSymMatrix& corr)
{ corrMatrix = corr; }


//inline const RealMatrix& MarginalsCorrDistribution::correlation_factor() const
//{ return corrCholeskyFactor; }


/** For consistent random variable ordering. */
inline void MarginalsCorrDistribution::
pull_distribution_parameters
(const std::shared_ptr<MultivariateDistribution> pull_mvd_rep)
{
  size_t v, num_rv = ranVarTypes.size();
  for (v=0; v<num_rv; ++v)
    pull_distribution_parameters(pull_mvd_rep, v, v);
}


inline void MarginalsCorrDistribution::
pull_distribution_parameters(const MultivariateDistribution& pull_mvd)
{ pull_distribution_parameters(pull_mvd.multivar_dist_rep()); }


/** For potentially inconsistent random variable ordering that requires
    a lookup. */
inline void MarginalsCorrDistribution::
pull_distribution_parameters
(const std::shared_ptr<MultivariateDistribution> pull_mvd_rep,
 const StringArray& pull_labels, const StringArray& push_labels)
{
  size_t v, num_rv = ranVarTypes.size();
  for (v=0; v<num_rv; ++v) {
    size_t push_index = find_index(push_labels, pull_labels[v]);
    if (push_index != _NPOS)
      pull_distribution_parameters(pull_mvd_rep, v, push_index);
  }
}


inline void MarginalsCorrDistribution::
pull_distribution_parameters(const MultivariateDistribution& pull_mvd,
			     const StringArray& pull_labels,
			     const StringArray& push_labels)
{
  pull_distribution_parameters(pull_mvd.multivar_dist_rep(),
			       pull_labels, push_labels);
}


inline void MarginalsCorrDistribution::
pull_distribution_parameters(const MultivariateDistribution& pull_mvd,
			     size_t pull_index, size_t push_index)
{
  pull_distribution_parameters(pull_mvd.multivar_dist_rep(),
			       pull_index, push_index);
}


template <typename ValueType>
void MarginalsCorrDistribution::
push_parameter(size_t v, short dist_param, ValueType value)
{ randomVars[v].push_parameter(dist_param, value); }


template <typename OrdinalType, typename ScalarType>
void MarginalsCorrDistribution::
push_parameters(size_t start_v, size_t num_v, short dist_param,
  const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values)
{
  // set one type of distribution parameter for a range of random variables
  // TO DO: would like to retire this version if Dakota migrates from Teuchos
  //        to std::vector for dist params

  size_t i, v, num_updates = std::min((size_t)values.length(), num_v);
  for (i=0, v=start_v; i<num_updates; ++i, ++v)
    randomVars[v].push_parameter(dist_param, values[i]);
}


template <typename ValueType>
void MarginalsCorrDistribution::
push_parameters(size_t start_v, size_t num_v, short dist_param,
		const std::vector<ValueType>& values)
{
  // set one distribution parameter type for a range of random variables

  size_t i, v, num_updates = std::min(values.size(), num_v);
  for (i=0, v=start_v; i<num_updates; ++i, ++v)
    randomVars[v].push_parameter(dist_param, values[i]);
}


template <typename OrdinalType, typename ScalarType>
void MarginalsCorrDistribution::
push_parameters(short rv_type, short dist_param,
  const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values)
{
  // rv_type eliminates need to check for dist_param support
  // TO DO: would like to retire this version if Dakota migrates from Teuchos
  //        to std::vector for dist params

  size_t v, num_v = ranVarTypes.size(), cntr = 0, num_vals = values.length();
  for (v=0; v < num_v && cntr < num_vals; ++v)
    if (ranVarTypes[v] == rv_type)
      randomVars[v].push_parameter(dist_param, values[cntr++]);
}


template <typename ValueType>
void MarginalsCorrDistribution::
push_parameters(short rv_type, short dist_param,
		const std::vector<ValueType>& values)
{
  // rv_type eliminates need to check for dist_param support

  size_t v, num_v = ranVarTypes.size(), cntr = 0, num_vals = values.size();
  for (v=0; v < num_v && cntr < num_vals; ++v)
    if (ranVarTypes[v] == rv_type)
      randomVars[v].push_parameter(dist_param, values[cntr++]);
}


template <typename ValueType>
void MarginalsCorrDistribution::
pull_parameter(size_t v, short dist_param, ValueType& val) const
{ randomVars[v].pull_parameter(dist_param, val); }


template <typename ValueType>
void MarginalsCorrDistribution::
pull_parameters(size_t start_v, size_t num_v, short dist_param,
		std::vector<ValueType>& values) const
{
  // set one distribution parameter type for a range of random variables

  size_t i, v;
  values.resize(num_v);
  for (i=0, v=start_v; i<num_v; ++i, ++v)
    randomVars[v].pull_parameter(dist_param, values[i]);
}


template <typename ValueType>
void MarginalsCorrDistribution::
pull_parameters(short rv_type, short dist_param,
		std::vector<ValueType>& values) const
{
  // rv_type eliminates need to check for dist_param support

  values.resize(std::count(ranVarTypes.begin(), ranVarTypes.end(), rv_type));

  size_t v, num_v = ranVarTypes.size(), cntr = 0;
  for (v=0; v<num_v; ++v)
    if (ranVarTypes[v] == rv_type)
      randomVars[v].pull_parameter(dist_param, values[cntr++]);
}


template <typename OrdinalType, typename ScalarType>
void MarginalsCorrDistribution::
pull_parameters(short rv_type, short dist_param,
  Teuchos::SerialDenseVector<OrdinalType, ScalarType>& values) const
{
  // rv_type eliminates need to check for dist_param support

  size_t v, num_v = ranVarTypes.size(), cntr = 0;
  values.sizeUninitialized(
    std::count(ranVarTypes.begin(), ranVarTypes.end(), rv_type));

  for (v=0; v<num_v; ++v)
    if (ranVarTypes[v] == rv_type)
      randomVars[v].pull_parameter(dist_param, values[cntr++]);
}


template <typename ValueType>
ValueType MarginalsCorrDistribution::
pull_parameter(size_t v, short dist_param) const
{
  ValueType val;
  randomVars[v].pull_parameter(dist_param, val);
  return val;
}


template <typename ValueType>
std::vector<ValueType> MarginalsCorrDistribution::
pull_parameters(short rv_type, short dist_param) const
{
  std::vector<ValueType> vals;
  pull_parameters<ValueType>(rv_type, dist_param, vals);
  return vals;
}


template <typename ValueType>
size_t MarginalsCorrDistribution::
pull_parameter_size(size_t v, short dist_param) const
{
  // overhead of one array copy instead of two or more
  ValueType val;
  randomVars[v].pull_parameter(dist_param, val);
  return val.size();
}


/* These APIs are not currently needed but could be useful in the future:
template <typename VectorType>
void MarginalsCorrDistribution::
push_parameters(size_t v, const ShortArray& dist_params,
                const VectorType& values)
{
  // set multiple distribution parameters for a single variable

  RandomVariable& random_var = randomVars[v];
  size_t i, num_params = std::min(dist_params.size(), values.length());
  for (i=0; i<num_params; ++i)
    random_var.push_parameter(dist_params[i], values[i]);
}


template <typename VectorType>
void MarginalsCorrDistribution::
push_parameters(short dist_param, const VectorType& values)
{
  // Without rv_type, query RV for dist_param support to match values to RV

  size_t v, num_v = randomVars.size(), cntr = 0, num_vals = values.length();
  for (v=0; v < num_v && cntr < num_vals; ++v)
    if (randomVars[v].supports(dist_param))
      randomVars[v].push_parameter(dist_param, values[cntr++]);
}
*/


inline RealRealPairArray MarginalsCorrDistribution::moments() const
{
  size_t v, num_v = randomVars.size();
  RealRealPairArray mom;
  if (activeVars.empty()) {
    mom.resize(num_v);
    for (v=0; v<num_v; ++v)
      mom[v] = randomVars[v].moments();
  }
  else {
    mom.resize(activeVars.count());
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (activeVars[v])
	mom[av_cntr++] = randomVars[v].moments();
  }
  return mom;
}


inline RealVector MarginalsCorrDistribution::means() const
{
  size_t v, num_v = randomVars.size();
  RealVector means;
  if (activeVars.empty()) {
    means.sizeUninitialized(num_v);
    for (v=0; v<num_v; ++v)
      means[v] = randomVars[v].mean();
  }
  else {
    means.sizeUninitialized(activeVars.count());
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (activeVars[v])
	means[av_cntr++] = randomVars[v].mean();
  }
  return means;
}


inline RealVector MarginalsCorrDistribution::std_deviations() const
{
  size_t v, num_v = randomVars.size();
  RealVector std_devs;
  if (activeVars.empty()) {
    std_devs.sizeUninitialized(num_v);
    for (v=0; v<num_v; ++v)
      std_devs[v] = randomVars[v].standard_deviation();
  }
  else {
    std_devs.sizeUninitialized(activeVars.count());
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (activeVars[v])
	std_devs[av_cntr++] = randomVars[v].standard_deviation();
  }
  return std_devs;
}


inline RealVector MarginalsCorrDistribution::variances() const
{
  size_t v, num_v = randomVars.size();
  RealVector vars;
  if (activeVars.empty()) {
    vars.sizeUninitialized(num_v);
    for (v=0; v<num_v; ++v)
      vars[v] = randomVars[v].variance();
  }
  else {
    vars.sizeUninitialized(activeVars.count());
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (activeVars[v])
	vars[av_cntr++] = randomVars[v].variance();
  }
  return vars;
}


inline RealRealPairArray MarginalsCorrDistribution::distribution_bounds() const
{
  size_t v, num_v = randomVars.size();
  RealRealPairArray bnds;
  if (activeVars.empty()) {
    bnds.resize(num_v);
    for (v=0; v<num_v; ++v)
      bnds[v] = randomVars[v].distribution_bounds();
  }
  else {
    bnds.resize(activeVars.count());
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (activeVars[v])
	bnds[av_cntr++] = randomVars[v].distribution_bounds();
  }
  return bnds;
}


inline RealVector MarginalsCorrDistribution::distribution_lower_bounds() const
{
  size_t v, num_v = randomVars.size();
  RealVector lwr_bnds;
  if (activeVars.empty()) {
    lwr_bnds.sizeUninitialized(num_v);
    for (v=0; v<num_v; ++v)
      lwr_bnds[v] = randomVars[v].distribution_bounds().first;
  }
  else {
    lwr_bnds.sizeUninitialized(activeVars.count());
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (activeVars[v])
	lwr_bnds[av_cntr++] = randomVars[v].distribution_bounds().first;
  }
  return lwr_bnds;
}


inline RealVector MarginalsCorrDistribution::distribution_upper_bounds() const
{
  size_t v, num_v = randomVars.size();
  RealVector upr_bnds;
  if (activeVars.empty()) {
    upr_bnds.sizeUninitialized(num_v);
    for (v=0; v<num_v; ++v)
      upr_bnds[v] = randomVars[v].distribution_bounds().second;
  }
  else {
    upr_bnds.sizeUninitialized(activeVars.count());
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (activeVars[v])
	upr_bnds[av_cntr++] = randomVars[v].distribution_bounds().second;
  }
  return upr_bnds;
}


inline bool MarginalsCorrDistribution::global_bounds() const
{ return globalBndsFlag; }


template <typename OrdinalType, typename ScalarType> 
void MarginalsCorrDistribution::check_active_length(
  const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& vec,
  const BitArray& mask) const
{
  size_t vec_len = vec.length(),
      expect_len = (mask.empty()) ? randomVars.size() : mask.count();
  if (vec_len != expect_len) {
    PCerr << "Error: bad active vector length (" << vec_len << "); "
	  << expect_len << "expected." << std::endl;
    abort_handler(-1);
  }
}


inline void MarginalsCorrDistribution::
lower_bounds(const RealVector& l_bnds, const BitArray& mask)
{
  check_active_length(l_bnds, mask);
  size_t v, num_v = randomVars.size();
  if (mask.empty())
    for (v=0; v<num_v; ++v)
      randomVars[v].lower_bound(l_bnds[v]);
  else {
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (mask[v])
	randomVars[v].lower_bound(l_bnds[av_cntr++]);
  }
}


inline void MarginalsCorrDistribution::
lower_bounds(const IntVector& l_bnds, const BitArray& mask)
{
  check_active_length(l_bnds, mask);
  size_t v, num_v = randomVars.size();
  if (mask.empty())
    for (v=0; v<num_v; ++v)
      randomVars[v].lower_bound(l_bnds[v]);
  else {
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (mask[v])
	randomVars[v].lower_bound(l_bnds[av_cntr++]);
  }
}


inline void MarginalsCorrDistribution::lower_bound(Real l_bnd, size_t rv_index)
{
  if (rv_index < randomVars.size())
    randomVars[rv_index].lower_bound(l_bnd);
  else {
    PCerr << "Error: rv_index (" << rv_index << ") out of range in Marginals"
	  << "CorrDistribution::lower_bound(Real, size_t)" << std::endl;
    abort_handler(-1);
  }
}


inline void MarginalsCorrDistribution::lower_bound(int l_bnd, size_t rv_index)
{
  if (rv_index < randomVars.size())
    randomVars[rv_index].lower_bound(l_bnd);
  else {
    PCerr << "Error: rv_index (" << rv_index << ") out of range in Marginals"
	  << "CorrDistribution::lower_bound(int, size_t)" << std::endl;
    abort_handler(-1);
  }
}


inline void MarginalsCorrDistribution::
upper_bounds(const RealVector& u_bnds, const BitArray& mask)
{
  check_active_length(u_bnds, mask);
  size_t v, num_v = randomVars.size();
  if (mask.empty())
    for (v=0; v<num_v; ++v)
      randomVars[v].upper_bound(u_bnds[v]);
  else {
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (mask[v])
	randomVars[v].upper_bound(u_bnds[av_cntr++]);
  }
}


inline void MarginalsCorrDistribution::
upper_bounds(const IntVector& u_bnds, const BitArray& mask)
{
  check_active_length(u_bnds, mask);
  size_t v, num_v = randomVars.size();
  if (mask.empty())
    for (v=0; v<num_v; ++v)
      randomVars[v].upper_bound(u_bnds[v]);
  else {
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (mask[v])
	randomVars[v].upper_bound(u_bnds[av_cntr++]);
  }
}


inline void MarginalsCorrDistribution::upper_bound(Real u_bnd, size_t rv_index)
{
  if (rv_index < randomVars.size())
    randomVars[rv_index].upper_bound(u_bnd);
  else {
    PCerr << "Error: rv_index (" << rv_index << ") out of range in Marginals"
	  << "CorrDistribution::upper_bound(Real, size_t)" << std::endl;
    abort_handler(-1);
  }
}


inline void MarginalsCorrDistribution::upper_bound(int u_bnd, size_t rv_index)
{
  if (rv_index < randomVars.size())
    randomVars[rv_index].upper_bound(u_bnd);
  else {
    PCerr << "Error: rv_index (" << rv_index << ") out of range in Marginals"
	  << "CorrDistribution::upper_bound(int, size_t)" << std::endl;
    abort_handler(-1);
  }
}


inline Real MarginalsCorrDistribution::pdf(Real val, size_t rv_index) const
{ return randomVars[rv_index].pdf(val); }


inline Real MarginalsCorrDistribution::
pdf_gradient(Real val, size_t rv_index) const
{ return randomVars[rv_index].pdf_gradient(val); }


inline Real MarginalsCorrDistribution::
pdf_hessian(Real val, size_t rv_index) const
{ return randomVars[rv_index].pdf_hessian(val); }


inline Real MarginalsCorrDistribution::log_pdf(Real val, size_t rv_index) const
{ return randomVars[rv_index].log_pdf(val); }


inline Real MarginalsCorrDistribution::
log_pdf_gradient(Real val, size_t rv_index) const
{ return randomVars[rv_index].log_pdf_gradient(val); }


inline Real MarginalsCorrDistribution::
log_pdf_hessian(Real val, size_t rv_index) const
{ return randomVars[rv_index].log_pdf_hessian(val); }


inline Real MarginalsCorrDistribution::pdf(const RealVector& pt) const
{
  // correlated density handled via upstream transform to standardized space
  if (correlationFlag) {
    PCerr << "Error: MarginalsCorrDistribution::pdf() currently uses a "
	  << "product of marginal densities\n       and can only be used for "
	  << "independent random variables." << std::endl;
    abort_handler(-1);
  }

  check_active_length(pt, activeVars);
  size_t v, num_v = randomVars.size();
  Real density = 1.;
  if (activeVars.empty())
    for (v=0; v<num_v; ++v)
      density *= pdf(pt[v], v);
  else {
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (activeVars[v])
	density *= pdf(pt[av_cntr++], v);
  }
  return density;
}


inline Real MarginalsCorrDistribution::log_pdf(const RealVector& pt) const
{
  // correlated density handled via upstream transform to standardized space
  if (correlationFlag) {
    PCerr << "Error: MarginalsCorrDistribution::log_pdf() currently uses a "
	  << "sum of log marginal densities\n       and can only be used for "
	  << "independent random variables." << std::endl;
    abort_handler(-1);
  }

  check_active_length(pt, activeVars);
  size_t v, num_v = randomVars.size();
  Real log_density = 0.;
  if (activeVars.empty())
    for (v=0; v<num_v; ++v)
      log_density += log_pdf(pt[v], v);
  else {
    size_t av_cntr = 0;
    for (v=0; v<num_v; ++v)
      if (activeVars[v])
	log_density += log_pdf(pt[av_cntr++], v);
  }
  return log_density;
}


template <typename Engine> 
Real MarginalsCorrDistribution::draw_sample(size_t rv_index, Engine& rng) const
{ return randomVars[rv_index].draw_sample(rng); }


template <typename Engine> 
Real MarginalsCorrDistribution::
draw_standard_sample(size_t rv_index, Engine& rng) const
{ return randomVars[rv_index].draw_standard_sample(rng); }

} // namespace Pecos

#endif
