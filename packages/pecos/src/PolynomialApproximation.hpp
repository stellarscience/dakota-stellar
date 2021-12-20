/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        PolynomialApproximation
//- Description:  Base Class for Orthogonal/Interpolation Polynomial
//-               Approximations
//-               
//- Owner:        Mike Eldred

#ifndef POLYNOMIAL_APPROXIMATION_HPP
#define POLYNOMIAL_APPROXIMATION_HPP

#include "BasisApproximation.hpp"
#include "SurrogateData.hpp"
#include "SharedPolyApproxData.hpp"

namespace Pecos {


/// Derived approximation class for global basis polynomials.

/** The PolynomialApproximation class provides a global approximation
    based on basis polynomials.  This includes orthogonal polynomials
    used for polynomial chaos expansions and interpolation polynomials
    used for stochastic collocation. */

class PolynomialApproximation: public BasisApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// standard constructor
  PolynomialApproximation(const SharedBasisApproxData& shared_data);
  /// destructorboth
  ~PolynomialApproximation();

  //
  //- Heading: Virtual member functions
  //

  /// update *Iter assignments for new active key
  virtual bool update_active_iterators(const ActiveKey& key);

  /// size derived class data attributes
  virtual void allocate_arrays() = 0;

  /// initialize product (covariance) accumulators with pointers to other QoI
  virtual void initialize_covariance(PolynomialApproximation* poly_approx_2);
  /// clear pointers to other QoI
  virtual void clear_covariance_pointers();
  /// (re)initialize and update product accumulators
  virtual void initialize_products();
  /// query whether product interpolants are defined (non-empty)
  virtual bool product_interpolants();
  /// clear bits indicating previously computed moment data
  virtual void clear_computed_bits();

  /// Computes sensitivity indices according to VBD specification
  virtual void compute_component_sobol() = 0;
  /// Computes total sensitivity indices according to VBD specification
  /// and existing computations from compute_component_sobol()
  virtual void compute_total_sobol() = 0;

  /// return RegressOrthogPolyApproximation::sparseSobolIndexMap
  virtual ULongULongMap sparse_sobol_index_map() const;
  /// return the number of non-zero expansion coefficients for this QoI
  virtual size_t expansion_terms() const;

  /// retrieve the gradient for a response expansion with respect to
  /// all variables included in the polynomial bases using the given
  /// parameter vector and default DVV
  virtual const RealVector& gradient_basis_variables(const RealVector& x) = 0;
  /// retrieve the gradient for a response expansion with respect to
  /// variables included in the polynomial basis for a given parameter
  /// vector and a given DVV subset
  virtual const RealVector& gradient_basis_variables(const RealVector& x,
						     const SizetArray& dvv) = 0;
  /// retrieve the gradient for a response expansion with respect to
  /// all variables not included in the polynomial bases
  /// (nonprobabilistic variables such as design or epistemic when not
  /// in "all" mode) using the given parameter vector and default DVV
  virtual const RealVector&
    gradient_nonbasis_variables(const RealVector& x) = 0;
  /// retrieve the Hessian of the response expansion with respect to all
  /// variables included in the polynomial basis (e.g., probabilistic
  /// variables) for a given parameter vector
  virtual const RealSymMatrix& hessian_basis_variables(const RealVector& x) = 0;

  /// retrieve the response value for a stored expansion using the
  /// given parameter vector
  virtual Real stored_value(const RealVector& x, const ActiveKey& key) = 0;
  /// retrieve the response gradient for a stored expansion with
  /// respect to all variables included in the polynomial bases;
  /// evaluate for the given parameter vector.
  virtual const RealVector& stored_gradient_basis_variables(const RealVector& x,
    const ActiveKey& key) = 0;
  /// retrieve the gradient for a stored expansion with respect to
  /// variables included in the polynomial basis for a given parameter
  /// vector and a given DVV subset
  virtual const RealVector& stored_gradient_basis_variables(const RealVector& x,
    const SizetArray& dvv, const ActiveKey& key) = 0;
  /// retrieve the response gradient for a stored expansion with
  /// respect to all variables not included in the polynomial bases;
  /// evaluate for the given parameter vector.
  virtual const RealVector& stored_gradient_nonbasis_variables(
    const RealVector& x, const ActiveKey& key) = 0;
  /// retrieve the Hessian for a stored expansion with respect to all
  /// variables included in the polynomial basis (e.g., probabilistic
  /// variables) for a given parameter vector
  virtual const RealSymMatrix& stored_hessian_basis_variables(
    const RealVector& x, const ActiveKey& key) = 0;

  /// return the mean of the expansion, treating all variables as random
  virtual Real mean() = 0;
  /// return the mean of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual Real mean(const RealVector& x) = 0;
  /// return the covariance between two response expansions, treating
  /// all variables as random
  virtual Real covariance(PolynomialApproximation* poly_approx_2) = 0;
  /// return the covariance between two response expansions for a given
  /// parameter vector, treating a subset of the variables as random
  virtual Real covariance(const RealVector& x,
			  PolynomialApproximation* poly_approx_2) = 0;

  /// return the gradient of the expansion mean for a given parameter
  /// vector, treating all variables as random
  virtual const RealVector& mean_gradient() = 0;
  /// return the gradient of the expansion mean for a given parameter vector
  /// and given DVV, treating a subset of the variables as random
  virtual const RealVector& mean_gradient(const RealVector& x,
					  const SizetArray& dvv) = 0;
  /// return the gradient of the expansion variance for a given parameter
  /// vector, treating all variables as random
  virtual const RealVector& variance_gradient() = 0;
  /// return the gradient of the expansion variance for a given parameter
  /// vector and given DVV, treating a subset of the variables as random
  virtual const RealVector& variance_gradient(const RealVector& x,
					      const SizetArray& dvv) = 0;

  /// return the mean of the combined expansion, treating all variables
  /// as random
  virtual Real combined_mean();
  /// return the mean of the combined expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual Real combined_mean(const RealVector& x);
  /// return the covariance between two combined response expansions,
  /// treating all variables as random
  virtual Real combined_covariance(PolynomialApproximation* poly_approx_2);
  /// return the covariance between two combined response expansions for a
  /// given parameter vector, treating a subset of the variables as random
  virtual Real combined_covariance(const RealVector& x,
				   PolynomialApproximation* poly_approx_2);

  /// return the reliability index (mapped from z_bar), where all active
  /// variables are random
  virtual Real beta(bool cdf_flag, Real z_bar);
  /// return the reliability index (mapped from z_bar), treating a subset of
  /// variables as random
  virtual Real beta(const RealVector& x, bool cdf_flag, Real z_bar);
  /// return the reliability index (mapped from z_bar), where all active
  /// variables are random
  virtual Real combined_beta(bool cdf_flag, Real z_bar);
  /// return the reliability index (mapped from z_bar), treating a subset of
  /// variables as random
  virtual Real combined_beta(const RealVector& x, bool cdf_flag, Real z_bar);

  /// return the change in mean resulting from expansion refinement,
  /// treating all variables as random
  virtual Real delta_mean();
  /// return the change in mean resulting from expansion refinement,
  /// treating a subset of the variables as random
  virtual Real delta_mean(const RealVector& x);
  /// return the change in mean resulting from combined expansion refinement,
  /// treating all variables as random
  virtual Real delta_combined_mean();
  /// return the change in mean resulting from combined expansion refinement,
  /// treating a subset of the variables as random
  virtual Real delta_combined_mean(const RealVector& x);

  /// return the change in standard deviation resulting from expansion
  /// refinement, treating all variables as random
  virtual Real delta_std_deviation();
  /// return the change in standard deviation resulting from expansion
  /// refinement, treating a subset of the variables as random
  virtual Real delta_std_deviation(const RealVector& x);
  /// return the change in standard deviation resulting from combined
  /// expansion refinement, treating all variables as random
  virtual Real delta_combined_std_deviation();
  /// return the change in standard deviation resulting from combined
  /// expansion refinement, treating a subset of the variables as random
  virtual Real delta_combined_std_deviation(const RealVector& x);

  /// return the change in covariance between two response expansions,
  /// treating all variables as random
  virtual Real delta_covariance(PolynomialApproximation* poly_approx_2);
  /// return the change in covariance between two response expansions for a
  /// given parameter vector, treating a subset of the variables as random
  virtual Real delta_covariance(const RealVector& x,
				PolynomialApproximation* poly_approx_2);
  /// return the change in covariance between two combined response expansions,
  /// treating all variables as random
  virtual Real
    delta_combined_covariance(PolynomialApproximation* poly_approx_2);
  /// return the change in covariance between two combined response expansions
  /// for given parameter vector, treating a subset of the variables as random
  virtual Real
    delta_combined_covariance(const RealVector& x,
			      PolynomialApproximation* poly_approx_2);

  /// return the change in reliability index (mapped from z_bar) resulting
  /// from expansion refinement, treating all variables as random
  virtual Real delta_beta(bool cdf_flag, Real z_bar);
  /// return the change in reliability index (mapped from z_bar) resulting
  /// from expansion refinement, treating a subset of the variables as random
  virtual Real delta_beta(const RealVector& x, bool cdf_flag, Real z_bar);
  /// return the change in reliability index (mapped from z_bar) resulting
  /// from expansion refinement, treating all variables as random
  virtual Real delta_combined_beta(bool cdf_flag, Real z_bar);
  /// return the change in reliability index (mapped from z_bar) resulting
  /// from expansion refinement, treating a subset of the variables as random
  virtual Real delta_combined_beta(const RealVector& x, bool cdf_flag,
				   Real z_bar);

  /// return the change in response level (mapped from beta_bar) resulting
  /// from expansion refinement, treating all variables as random
  virtual Real delta_z(bool cdf_flag, Real beta_bar);
  /// return the change in response level (mapped from beta_bar) resulting
  /// from expansion refinement, treating a subset of the variables as random
  virtual Real delta_z(const RealVector& x, bool cdf_flag, Real beta_bar);
  /// return the change in response level (mapped from beta_bar) resulting
  /// from expansion refinement, treating all variables as random
  virtual Real delta_combined_z(bool cdf_flag, Real beta_bar);
  /// return the change in response level (mapped from beta_bar) resulting
  /// from expansion refinement, treating a subset of the variables as random
  virtual Real delta_combined_z(const RealVector& x, bool cdf_flag,
				Real beta_bar);

  /// compute central response moments using some combination of expansion
  /// post-processing and numerical integration
  virtual void compute_moments(bool full_stats = true,
			       bool combined_stats = false);
  /// compute central response moments in all-variables mode using some
  /// combination of expansion post-processing and numerical integration
  virtual void compute_moments(const RealVector& x, bool full_stats = true,
			       bool combined_stats = false);

  /// return moments computed analytically from the expansion
  /// (corresponding to active key if primary)
  virtual const RealVector& expansion_moments() const = 0;
  /// return moments computed numerically by (sparse) quadrature
  /// (corresponding to active key if primary)
  virtual const RealVector& numerical_integration_moments() const = 0;

  //
  //- Heading: Member functions
  //

  /// return preferred response moments (either expansion or numerical
  /// integration, depending on approximation type)
  const RealVector& moments() const;
  /// set preferred response moments (either expansion or numerical
  /// integration, depending on approximation type); this is generally
  /// used to restore previous values when popping an adaptation, without
  /// having to recompute them
  void moments(const RealVector& mom);
  /// return preferred response moments (either expansion or numerical
  /// integration, depending on approximation type)
  Real moment(size_t i) const;
  /// set preferred response moments (either expansion or numerical
  /// integration, depending on approximation type); this is generally
  /// used to restore previous values when popping an adaptation, without
  /// having to recompute them
  void moment(Real mom, size_t i);

  /// return combined response moments
  const RealVector& combined_moments() const;
  /// set combined response moments; this is generally used to restore previous
  /// values when popping an adaptation, without having to recompute
  void combined_moments(const RealVector& mom);
  /// return i-th combined response moment
  Real combined_moment(size_t i) const;
  /// set i-th combined response moment; this is generally used to restore
  /// previous values when popping an adaptation, without having to recompute
  void combined_moment(Real mom, size_t i);

  /// standardize central moments 2-n and eliminate excess kurtosis
  static void standardize_moments(const RealVector& central_moments,
				  RealVector& std_moments);

  /// return the variance of the expansion, treating all variables as random
  Real variance();
  /// return the variance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  Real variance(const RealVector& x);
  /// return the variance of the combined expansion, treating all
  /// variables as random
  Real combined_variance();
  /// return the variance of the combined expansion for a given parameter
  /// vector, treating a subset of the variables as random
  Real combined_variance(const RealVector& x);

  /// return the change in the variance of the expansion, treating all
  /// variables as random
  Real delta_variance();
  /// return the change in the variance of the expansion for a given
  /// parameter vector, treating a subset of the variables as random
  Real delta_variance(const RealVector& x);
  /// return the change in the variance of the combined expansion,
  /// treating all variables as random
  Real delta_combined_variance();
  /// return the change in the variance of the combined expansion for a
  /// given parameter vector, treating a subset of the variables as random
  Real delta_combined_variance(const RealVector& x);

  // number of data points to remove in a decrement
  //size_t pop_count();

  /// set ExpansionConfigOptions::expansionCoeffFlag
  void expansion_coefficient_flag(bool coeff_flag);
  /// get ExpansionConfigOptions::expansionCoeffFlag
  bool expansion_coefficient_flag() const;

  /// set ExpansionConfigOptions::expansionCoeffGradFlag
  void expansion_coefficient_gradient_flag(bool grad_flag);
  /// get ExpansionConfigOptions::expansionCoeffGradFlag
  bool expansion_coefficient_gradient_flag() const;

  /// clear Sobol' indices when inactive
  void clear_component_sobol();

  /// return sobolIndices
  const RealVector& sobol_indices() const;
  /// return totalSobolIndices
  const RealVector& total_sobol_indices() const;

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void surrogate_data(const SurrogateData& data);
  const SurrogateData& surrogate_data() const;
  SurrogateData& surrogate_data();

  void compute_coefficients();
  void combined_to_active(bool clear_combined = true);

  /// generic base class function mapped to gradient_basis_variables(x)
  const RealVector& gradient(const RealVector& x);
  /// generic base class function mapped to hessian_basis_variables(x)
  const RealSymMatrix& hessian(const RealVector& x);

  //
  //- Heading: Member functions
  //

  /// update surrData to define aggregated data from raw data, when indicated
  /// by an active aggregated key
  void synchronize_surrogate_data();
  /// generate synthetic data for the surrogate QoI prediction corresponding
  /// to the level key preceding active key; for use in surplus estimation
  /// for new level data relative to a previous level's surrogate prediction
  void generate_synthetic_data(SurrogateData& surr_data,
			       const ActiveKey& active_key,
			       short combine_type);

  /// compute central moments of response using type1 numerical integration
  void integrate_moments(const RealVector& coeffs, const RealVector& t1_wts,
			 RealVector& moments);
  /// compute central moments of response using type1/2 numerical integration
  void integrate_moments(const RealVector& t1_coeffs,
			 const RealMatrix& t2_coeffs, const RealVector& t1_wts,
			 const RealMatrix& t2_wts, RealVector& moments);

  /// size component Sobol arrays
  void allocate_component_sobol();
  /// size total Sobol arrays
  void allocate_total_sobol();

  /// zero out bit trackers for active moments
  void clear_active_bits();
  /// zero out bit trackers for combined moments
  void clear_combined_bits();

  //
  //- Heading: Data
  //

  /// SurrogateData instance containing the variables (shared) and response
  /// (unique) data arrays for constructing a surrogate of a single response
  /// function; this includes the original unmodified data set for one or more
  /// level keys as well as combinations for aggregated keys (e.g., additive
  /// or multiplicative discrepancies)
  SurrogateData surrData;

  /// flag for calculation of expansion coefficients from response values
  bool expansionCoeffFlag;
  /// flag for calculation of gradients of expansion coefficients from
  /// response gradients
  bool expansionCoeffGradFlag;

  /// gradient of the polynomial approximation returned by gradient()
  RealVector    approxGradient;
  /// Hessian of the polynomial approximation returned by hessian()
  RealSymMatrix approxHessian;

  /// mean and central moments 2/3/4 computed from either the expansion form
  /// (OrthogPolyApproximation) or via numerical integration of the response
  /// (InterpPolyApproximation).  Conversions to standardized moments
  /// (std deviation, skewness, kurtosis) are performed elsewhere as needed.
  std::map<ActiveKey, RealVector> primaryMoments;
  /// iterator to active entry in primaryMoments
  std::map<ActiveKey, RealVector>::iterator primaryMomIter;
  /// track computation of mean and mean gradient to avoid unnecessary
  /// recomputation
  std::map<ActiveKey, short> primaryMeanBits;
  /// iterator to active entry in primaryMeanBits
  std::map<ActiveKey, short>::iterator primaryMeanIter;
  /// track computation of variance and variance gradient to avoid
  /// unnecessary recomputation
  std::map<ActiveKey, short> primaryVarBits;
  /// iterator to active entry in primaryVarBits
  std::map<ActiveKey, short>::iterator primaryVarIter;
  /// track previous evaluation point for all_variables mean to avoid
  /// unnecessary recomputation
  std::map<ActiveKey, RealVector> xPrevMean;
  /// track previous evaluation point for all_variables variance to
  /// avoid unnecessary recomputation
  std::map<ActiveKey, RealVector> xPrevVar;

  /// gradient of mean/variance/etc. for the primary moments (expansion
  /// for OrthogPoly, numerical for InterpPoly) 
  std::map<ActiveKey, RealVectorArray> primaryMomentGrads;
  /// iterator to active entry in primaryMomentGrads
  std::map<ActiveKey, RealVectorArray>::iterator primaryMomGradsIter;
  /// track previous evaluation point for all_variables mean gradient
  /// to avoid unnecessary recomputation
  std::map<ActiveKey, RealVector> xPrevMeanGrad;
  /// track previous evaluation point for all_variables variance
  /// gradient to avoid unnecessary recomputation
  std::map<ActiveKey, RealVector> xPrevVarGrad;

  /// alternate non-active moments (numerical for OrthogPolyApproximation or
  /// expansion for InterpPolyApproximation).  These are computed in final
  /// post-processing and do not currently require key management.
  RealVector secondaryMoments;

  /// moments resulting from expansion roll-up across model index keys
  RealVector combinedMoments;
  /// track computation of combined mean and combined mean gradient to
  /// avoid unnecessary recomputation
  short combinedMeanBits;
  /// track computation of combined variance and combined variance
  /// gradient to avoid unnecessary recomputation
  short combinedVarBits;
  /// track previous evaluation point for all_variables combined mean to avoid
  /// unnecessary recomputation
  RealVector xPrevCombMean;
  /// track previous evaluation point for all_variables combined variance to
  /// avoid unnecessary recomputation
  RealVector xPrevCombVar;

  /// global sensitivities as given by Sobol'
  RealVector sobolIndices;
  /// total global sensitivities as given by Sobol'
  RealVector totalSobolIndices;

private:

  //
  //- Heading: Data
  //
};


inline PolynomialApproximation::
PolynomialApproximation(const SharedBasisApproxData& shared_data):
  BasisApproximation(BaseConstructor(), shared_data), expansionCoeffFlag(true),
  expansionCoeffGradFlag(false), primaryMomIter(primaryMoments.end()),
  combinedMeanBits(0), combinedVarBits(0)
{ }


inline PolynomialApproximation::~PolynomialApproximation()
{ }


inline bool PolynomialApproximation::
update_active_iterators(const ActiveKey& key)
{
  if (primaryMomIter != primaryMoments.end() && primaryMomIter->first == key)
    return false;

  primaryMomIter      = primaryMoments.find(key);
  primaryMomGradsIter = primaryMomentGrads.find(key);
  primaryMeanIter     = primaryMeanBits.find(key);
  primaryVarIter      = primaryVarBits.find(key);

  // share 1 deep copy of current active key
  ActiveKey key_copy;
  if (primaryMomIter      == primaryMoments.end()     ||
      primaryMomGradsIter == primaryMomentGrads.end() ||
      primaryMeanIter     == primaryMeanBits.end()    ||
      primaryVarIter      == primaryVarBits.end())
    key_copy = key.copy();

  if (primaryMomIter == primaryMoments.end()) {
    std::pair<ActiveKey, RealVector> rv_pair(key_copy, RealVector());
    primaryMomIter = primaryMoments.insert(rv_pair).first;
  }

  //if (expansionCoeffGradFlag) { // or all_vars for combined exp grads?
    if (primaryMomGradsIter == primaryMomentGrads.end()) {
      std::pair<ActiveKey, RealVectorArray>
	rva_pair(key_copy, RealVectorArray(2));
      primaryMomGradsIter = primaryMomentGrads.insert(rva_pair).first;
    }
  //}

  if (primaryMeanIter == primaryMeanBits.end()) {
    std::pair<ActiveKey, short> us_pair(key_copy, 0);
    primaryMeanIter = primaryMeanBits.insert(us_pair).first;
  }
  if (primaryVarIter == primaryVarBits.end()) {
    std::pair<ActiveKey, short> us_pair(key_copy, 0);
    primaryVarIter = primaryVarBits.insert(us_pair).first;
  }

  return true;
}


inline const SurrogateData& PolynomialApproximation::surrogate_data() const
{ return surrData; }


inline SurrogateData& PolynomialApproximation::surrogate_data()
{ return surrData; }


inline void PolynomialApproximation::surrogate_data(const SurrogateData& data)
{ surrData = data; }


inline void PolynomialApproximation::clear_active_bits()
{ primaryMeanIter->second = primaryVarIter->second = 0; }


inline void PolynomialApproximation::clear_combined_bits()
{ combinedMeanBits        = combinedVarBits        = 0; }


inline void PolynomialApproximation::clear_computed_bits()
{
  //SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  //if (data_rep->expConfigOptions.refineStatsType == COMBINED_EXPANSION_STATS)
    clear_combined_bits();
  //else
    clear_active_bits();
}


inline const RealVector& PolynomialApproximation::moments() const
{ return primaryMomIter->second; }


inline void PolynomialApproximation::moments(const RealVector& mom)
{
  primaryMomIter->second = mom;//copy_data_partial(mom, primaryMoments, 0);

  // activate bit trackers for moment set (a truncated array _does_ invalidate
  // other moments as this is assumed to be the available moment set)
  switch (mom.length()) {
  case 0:  primaryMeanIter->second &= ~1; primaryVarIter->second &= ~1; break;
  case 1:  primaryMeanIter->second |=  1; primaryVarIter->second &= ~1; break;
  default: primaryMeanIter->second |=  1; primaryVarIter->second |=  1; break;
  }
}


inline Real PolynomialApproximation::moment(size_t i) const
{
  const RealVector& moments = primaryMomIter->second;
  if (moments.length() <= i) {
    PCerr << "Error: index (" << i << ") out of bounds in Polynomial"
	  << "Approximation::moment()." << std::endl;
    abort_handler(-1);
  }
  return moments[i];
}


inline void PolynomialApproximation::moment(Real mom, size_t i)
{
  // see discussion above regarding future consolidation

  RealVector& moments = primaryMomIter->second;
  if (moments.length() <= i) {
    PCerr << "Error: index (" << i << ") out of bounds in Polynomial"
	  << "Approximation::moment()." << std::endl;
    abort_handler(-1);
  }
  moments[i] = mom;

  // activate bit tracker for i-th moment (does _not_ invalidate other moments)
  switch (i) {
  case 0: primaryMeanIter->second |= 1; break;
  case 1: primaryVarIter->second  |= 1; break;
  //default: other indices do not have trackers
  }
}


inline const RealVector& PolynomialApproximation::combined_moments() const
{ return combinedMoments; }


inline void PolynomialApproximation::combined_moments(const RealVector& mom)
{
  combinedMoments = mom;//copy_data_partial(mom, primaryMoments, 0);

  // activate bit trackers for moment set (a truncated array _does_ invalidate
  // other moments as this is assumed to be the available moment set)
  switch (mom.length()) {
  case 0:  combinedMeanBits &= ~1; combinedVarBits &= ~1; break;
  case 1:  combinedMeanBits |=  1; combinedVarBits &= ~1; break;
  default: combinedMeanBits |=  1; combinedVarBits |=  1; break;//usual
  }
}


inline Real PolynomialApproximation::combined_moment(size_t i) const
{
  if (combinedMoments.length() <= i) {
    PCerr << "Error: index (" << i << ") out of bounds in Polynomial"
	  << "Approximation::combined_moment()." << std::endl;
    abort_handler(-1);
  }
  return combinedMoments[i];
}


inline void PolynomialApproximation::combined_moment(Real mom, size_t i)
{
  // see discussion above regarding future consolidation

  if (combinedMoments.length() <= i) {
    PCerr << "Error: index (" << i << ") out of bounds in Polynomial"
	  << "Approximation::moment()." << std::endl;
    abort_handler(-1);
  }
  combinedMoments[i] = mom;

  // activate bit tracker for i-th moment (does _not_ invalidate other moments)
  switch (i) {
  case 0: combinedMeanBits |= 1; break;
  case 1: combinedVarBits  |= 1; break;
  //default: other indices do not have trackers
  }
}


inline Real PolynomialApproximation::variance()
{ return covariance(this); }


inline Real PolynomialApproximation::variance(const RealVector& x)
{ return covariance(x, this); }


inline Real PolynomialApproximation::combined_variance()
{ return combined_covariance(this); }


inline Real PolynomialApproximation::combined_variance(const RealVector& x)
{ return combined_covariance(x, this); }


inline Real PolynomialApproximation::delta_variance()
{ return delta_covariance(this); }


inline Real PolynomialApproximation::delta_variance(const RealVector& x)
{ return delta_covariance(x, this); }


inline Real PolynomialApproximation::delta_combined_variance()
{ return delta_combined_covariance(this); }


inline Real PolynomialApproximation::
delta_combined_variance(const RealVector& x)
{ return delta_combined_covariance(x, this); }


//inline size_t PolynomialApproximation::pop_count()
//{
//  SharedPolyApproxData* spad_rep = (SharedPolyApproxData*)sharedDataRep;
//  SparseGridDriver* ssg_driver = (SparseGridDriver*)(spad_rep->driverRep);
//  return (size_t)ssg_driver->unique_trial_points();
//}


inline void PolynomialApproximation::expansion_coefficient_flag(bool coeff_flag)
{ expansionCoeffFlag = coeff_flag; }


inline bool PolynomialApproximation::expansion_coefficient_flag() const
{ return expansionCoeffFlag; }


inline void PolynomialApproximation::
expansion_coefficient_gradient_flag(bool grad_flag)
{ expansionCoeffGradFlag = grad_flag; }


inline bool PolynomialApproximation::
expansion_coefficient_gradient_flag() const
{ return expansionCoeffGradFlag; }


inline void PolynomialApproximation::clear_component_sobol()
{
  //if (data_rep->expConfigOptions.vbdFlag)
    sobolIndices = 0.;
  //else
  //  sobolIndices.resize(0);
}


inline const RealVector& PolynomialApproximation::sobol_indices() const
{ return sobolIndices; }


inline const RealVector& PolynomialApproximation::total_sobol_indices() const
{ return totalSobolIndices; }


inline const RealVector& PolynomialApproximation::gradient(const RealVector& x)
{ return gradient_basis_variables(x); }


inline const RealSymMatrix& PolynomialApproximation::
hessian(const RealVector& x)
{ return hessian_basis_variables(x); }

} // namespace Pecos

#endif
