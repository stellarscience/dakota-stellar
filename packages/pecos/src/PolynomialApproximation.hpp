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
  virtual Real stored_value(const RealVector& x, const UShortArray& key) = 0;
  /// retrieve the response gradient for a stored expansion with
  /// respect to all variables included in the polynomial bases;
  /// evaluate for the given parameter vector.
  virtual const RealVector& stored_gradient_basis_variables(const RealVector& x,
    const UShortArray& key) = 0;
  /// retrieve the gradient for a stored expansion with respect to
  /// variables included in the polynomial basis for a given parameter
  /// vector and a given DVV subset
  virtual const RealVector& stored_gradient_basis_variables(const RealVector& x,
    const SizetArray& dvv, const UShortArray& key) = 0;
  /// retrieve the response gradient for a stored expansion with
  /// respect to all variables not included in the polynomial bases;
  /// evaluate for the given parameter vector.
  virtual const RealVector& stored_gradient_nonbasis_variables(
    const RealVector& x, const UShortArray& key) = 0;
  /// retrieve the Hessian for a stored expansion with respect to all
  /// variables included in the polynomial basis (e.g., probabilistic
  /// variables) for a given parameter vector
  virtual const RealSymMatrix& stored_hessian_basis_variables(
    const RealVector& x, const UShortArray& key) = 0;

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

  /// return the mean of the expansion, treating all variables as random
  virtual Real combined_mean();
  /// return the mean of the expansion for a given parameter vector,
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
			       bool combined_stats = false) = 0;
  /// compute central response moments in all-variables mode using some
  /// combination of expansion post-processing and numerical integration
  virtual void compute_moments(const RealVector& x, bool full_stats = true,
			       bool combined_stats = false) = 0;

  //
  //- Heading: Member functions
  //

  /// return active expansionMoments
  const RealVector& expansion_moments() const;
  /// return active numericalMoments
  const RealVector& numerical_integration_moments() const;
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

  void modified_surrogate_data(const SurrogateData& data);
  const SurrogateData& modified_surrogate_data() const;
  SurrogateData& modified_surrogate_data();

  void compute_coefficients();

  /// generic base class function mapped to gradient_basis_variables(x)
  const RealVector& gradient(const RealVector& x);
  /// generic base class function mapped to hessian_basis_variables(x)
  const RealSymMatrix& hessian(const RealVector& x);

  //
  //- Heading: New virtual functions
  //

  /// update *Iter for new activeKey from sharedDataRep
  virtual bool update_active_iterators(const UShortArray& key);

  //
  //- Heading: Member functions
  //

  /// update modSurrData from surrData based on deep or shallow copy
  void synchronize_surrogate_data();
  /// compute hierarchical surpluses from the surrData active key
  /// and store in same key within modSurrData
  void response_data_to_surplus_data();
  /// compute discrepancy data using surrData keys (HF and LF pairs)
  /// and store in modSurrData
  void response_data_to_discrepancy_data();

  /// define a LF key corresponding to incoming HF key
  void paired_lf_key(const UShortArray& hf_key, UShortArray& lf_key) const;

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

  //
  //- Heading: Data
  //

  /// SurrogateData instance containing the variables (shared) and response
  /// (unique) data arrays for constructing a surrogate of a single response
  /// function; this is the original unmodified data set for one or more level
  /// keys, prior to any potential manipulations by the approximation classes.
  SurrogateData surrData;
  /// SurrogateData instance used in current approximation builds, potentially
  /// reflecting data modifications (e.g., calculation of hierarchical surplus
  /// of surrData level relative to previous surrogate level or model
  /// discrepancy between two consecutive surrData level keys)
  SurrogateData modSurrData;

  /// flag for calculation of expansion coefficients from response values
  bool expansionCoeffFlag;
  /// flag for calculation of gradients of expansion coefficients from
  /// response gradients
  bool expansionCoeffGradFlag;

  /// mean and central moments 2/3/4 computed from the stochastic expansion
  /// form.  For OrthogPolyApproximation, these are primary, and for
  /// InterpPolyApproximation, they are secondary.  Conversions to standardized
  /// moments (std deviation, skewness, kurtosis) are performed elsewhere.
  std::map<UShortArray, RealVector> expansionMoments;
  /// iterator to active entry in expansionMoments
  std::map<UShortArray, RealVector>::iterator expMomentsIter;
  /// mean and central moments 2/3/4 computed via numerical integration of the
  /// response.  For OrthogPolyApproximation, these are secondary, and for
  /// InterpPolyApproximation, they are primary.  Conversions to standardized
  /// moments (std deviation, skewness, kurtosis) are performed elsewhere.
  std::map<UShortArray, RealVector> numericalMoments;
  /// iterator to active entry in numericalMoments
  std::map<UShortArray, RealVector>::iterator numMomentsIter;

  /// gradient of the polynomial approximation returned by gradient()
  RealVector approxGradient;
  /// Hessian of the polynomial approximation returned by hessian()
  RealSymMatrix approxHessian;
  /// gradient of mean/variance/etc. for the "primary" moments (expansion
  /// for OrthogPoly, numerical for InterpPoly) 
  std::map<UShortArray, RealVectorArray> momentGradients;
  /// iterator to active entry in momentGradients
  std::map<UShortArray, RealVectorArray>::iterator momentGradsIter;

  /// track computation of mean and mean gradient to avoid unnecessary
  /// recomputation
  std::map<UShortArray, short> computedMean;
  /// iterator to active entry in computedMean
  std::map<UShortArray, short>::iterator compMeanIter;
  /// track computation of variance and variance gradient to avoid
  /// unnecessary recomputation
  std::map<UShortArray, short> computedVariance;
  /// iterator to active entry in computedVariance
  std::map<UShortArray, short>::iterator compVarIter;
  /// track previous evaluation point for all_variables mean to avoid
  /// unnecessary recomputation
  std::map<UShortArray, RealVector> xPrevMean;
  /// track previous evaluation point for all_variables mean gradient
  /// to avoid unnecessary recomputation
  std::map<UShortArray, RealVector> xPrevMeanGrad;
  /// track previous evaluation point for all_variables variance to
  /// avoid unnecessary recomputation
  std::map<UShortArray, RealVector> xPrevVar;
  /// track previous evaluation point for all_variables variance
  /// gradient to avoid unnecessary recomputation
  std::map<UShortArray, RealVector> xPrevVarGrad;

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
  BasisApproximation(BaseConstructor(), shared_data),
  expansionCoeffFlag(true), expansionCoeffGradFlag(false),
  expMomentsIter(expansionMoments.end()), numMomentsIter(numericalMoments.end())
{ }


inline PolynomialApproximation::~PolynomialApproximation()
{ }


inline bool PolynomialApproximation::
update_active_iterators(const UShortArray& key)
{
  if (expMomentsIter != expansionMoments.end() && expMomentsIter->first == key)
    return false;

  expMomentsIter = expansionMoments.find(key);
  if (expMomentsIter == expansionMoments.end()) {
    std::pair<UShortArray, RealVector> rv_pair(key, RealVector());
    expMomentsIter = expansionMoments.insert(rv_pair).first;
  }
  numMomentsIter = numericalMoments.find(key);
  if (numMomentsIter == numericalMoments.end()) {
    std::pair<UShortArray, RealVector> rv_pair(key, RealVector());
    numMomentsIter = numericalMoments.insert(rv_pair).first;
  }

  //if (expansionCoeffGradFlag) { // or all_vars for combined exp grads?
    momentGradsIter = momentGradients.find(key);
    if (momentGradsIter == momentGradients.end()) {
      std::pair<UShortArray, RealVectorArray> rva_pair(key, RealVectorArray(2));
      momentGradsIter = momentGradients.insert(rva_pair).first;
    }
  //}

  compMeanIter = computedMean.find(key);
  if (compMeanIter == computedMean.end()) {
    std::pair<UShortArray, short> us_pair(key, 0);
    compMeanIter = computedMean.insert(us_pair).first;
  }
  compVarIter = computedVariance.find(key);
  if (compVarIter == computedVariance.end()) {
    std::pair<UShortArray, short> us_pair(key, 0);
    compVarIter = computedVariance.insert(us_pair).first;
  }

  return true;
}


inline const SurrogateData& PolynomialApproximation::surrogate_data() const
{ return surrData; }


inline SurrogateData& PolynomialApproximation::surrogate_data()
{ return surrData; }


inline void PolynomialApproximation::surrogate_data(const SurrogateData& data)
{ surrData = data; }


inline const SurrogateData& PolynomialApproximation::
modified_surrogate_data() const
{ return modSurrData; }


inline SurrogateData& PolynomialApproximation::modified_surrogate_data()
{ return modSurrData; }


inline void PolynomialApproximation::
modified_surrogate_data(const SurrogateData& data)
{ modSurrData = data; } // shared rep


inline void PolynomialApproximation::
paired_lf_key(const UShortArray& hf_key, UShortArray& lf_key) const
{
  if (hf_key.back() > 0) {
    // decrement trailing index
    lf_key = hf_key; --lf_key.back();
    // append the HF key in order to tag a particular (discrepancy) pairing
    lf_key.insert(lf_key.end(), hf_key.begin(), hf_key.end());
  }
  else
    lf_key.clear();
}


inline void PolynomialApproximation::clear_computed_bits()
{ compMeanIter->second = compVarIter->second = 0; }


inline const RealVector& PolynomialApproximation::expansion_moments() const
{ return expMomentsIter->second; }


inline const RealVector& PolynomialApproximation::
numerical_integration_moments() const
{ return numMomentsIter->second; }


/** All current cases should prefer expansionMoments (with the distinction
    drawn between interpolation and integration rules, we prefer Gauss
    integration rules on interpolant expansions). */
inline const RealVector& PolynomialApproximation::moments() const
{
  // TO DO: return to this (will require activation of exp moments for all view)
  //return expMomentsIter->second;

  // TO DO: remove temporary approach (avoids derived class specialization)
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expBasisType) {
  case NODAL_INTERPOLANT: case HIERARCHICAL_INTERPOLANT:
    return numMomentsIter->second; break;
  default:
    return expMomentsIter->second; break;
  }
}


inline void PolynomialApproximation::moments(const RealVector& mom)
{
  // see discussion above regarding future consolidation

  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expBasisType) {
  case NODAL_INTERPOLANT: case HIERARCHICAL_INTERPOLANT:
    numMomentsIter->second = mom;//copy_data_partial(mom, numericalMoments, 0);
    break;
  default:
    expMomentsIter->second = mom;//copy_data_partial(mom, expansionMoments, 0);
    break;
  }

  // activate bit trackers for moment set (a truncated array _does_ invalidate
  // other moments as this is assumed to be the available moment set)
  size_t num_mom = mom.length();
  switch (num_mom) {
  case 0:  compMeanIter->second &= ~1;  compVarIter->second &= ~1; break;
  case 1:  compMeanIter->second |=  1;  compVarIter->second &= ~1; break;
  default: compMeanIter->second |=  1;  compVarIter->second |=  1; break;//usual
  }
}


inline Real PolynomialApproximation::moment(size_t i) const
{
  // see discussion above regarding future consolidation

  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  short exp_type = data_rep->expConfigOptions.expBasisType;
  const RealVector& moments
    = (exp_type == NODAL_INTERPOLANT || exp_type == HIERARCHICAL_INTERPOLANT)
    ? numMomentsIter->second : expMomentsIter->second;
  if (moments.length() <= i) {
    PCerr << "Error: index (" << i << ") out of bounds in "
	  << "PolynomialApproximation::moment()." << std::endl;
    abort_handler(-1);
  }
  return moments[i];
}


inline void PolynomialApproximation::moment(Real mom, size_t i)
{
  // see discussion above regarding future consolidation

  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  short exp_type = data_rep->expConfigOptions.expBasisType;
  RealVector& moments
    = (exp_type == NODAL_INTERPOLANT || exp_type == HIERARCHICAL_INTERPOLANT)
    ? numMomentsIter->second : expMomentsIter->second;
  if (moments.length() <= i) {
    PCerr << "Error: index (" << i << ") out of bounds in "
	  << "PolynomialApproximation::moment()." << std::endl;
    abort_handler(-1);
  }
  moments[i] = mom;

  // activate bit tracker for i-th moment (does _not_ invalidate other moments)
  switch (i) {
  case 0: compMeanIter->second |= 1; break;
  case 1: compVarIter->second  |= 1; break;
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
