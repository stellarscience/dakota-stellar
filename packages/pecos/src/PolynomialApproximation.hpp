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

class AleatoryDistParams;


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
  /// size component Sobol arrays
  virtual void allocate_component_sobol();

  /// Computes sensitivity indices according to VBD specification
  virtual void compute_component_sobol() = 0;
  /// Computes total sensitivity indices according to VBD specification
  /// and existing computations from compute_component_sobol()
  virtual void compute_total_sobol() = 0;

  /// return RegressOrthogPolyApproximation::sparseSobolIndexMap
  virtual ULongULongMap sparse_sobol_index_map() const;

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
  virtual Real stored_value(const RealVector& x, size_t index) = 0;
  /// retrieve the response gradient for a stored expansion with
  /// respect to all variables included in the polynomial bases;
  /// evaluate for the given parameter vector.
  virtual const RealVector&
    stored_gradient_basis_variables(const RealVector& x, size_t index) = 0;
  /// retrieve the response gradient for a stored expansion with
  /// respect to all variables not included in the polynomial bases;
  /// evaluate for the given parameter vector.
  virtual const RealVector&
    stored_gradient_nonbasis_variables(const RealVector& x, size_t index) = 0;

  /// return the mean of the expansion, treating all variables as random
  virtual Real mean() = 0;
  /// return the mean of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual Real mean(const RealVector& x) = 0;
  /// return the gradient of the expansion mean for a given parameter
  /// vector, treating all variables as random
  virtual const RealVector& mean_gradient() = 0;
  /// return the gradient of the expansion mean for a given parameter vector
  /// and given DVV, treating a subset of the variables as random
  virtual const RealVector& mean_gradient(const RealVector& x,
					  const SizetArray& dvv) = 0;

  /// return the variance of the expansion, treating all variables as random
  virtual Real variance() = 0;
  /// return the variance of the expansion for a given parameter vector,
  /// treating a subset of the variables as random
  virtual Real variance(const RealVector& x) = 0;
  /// return the gradient of the expansion variance for a given parameter
  /// vector, treating all variables as random
  virtual const RealVector& variance_gradient() = 0;
  /// return the gradient of the expansion variance for a given parameter
  /// vector and given DVV, treating a subset of the variables as random
  virtual const RealVector& variance_gradient(const RealVector& x,
					      const SizetArray& dvv) = 0;

  /// return the covariance between two response expansions, treating
  /// all variables as random
  virtual Real covariance(PolynomialApproximation* poly_approx_2) = 0;
  /// return the covariance between two response expansions for a given
  /// parameter vector, treating a subset of the variables as random
  virtual Real covariance(const RealVector& x,
			  PolynomialApproximation* poly_approx_2) = 0;

  /// return the change in covariance between two response expansions,
  /// treating all variables as random
  virtual Real delta_covariance(PolynomialApproximation* poly_approx_2);
  /// return the change in covariance between two response expansions for a
  /// given parameter vector, treating a subset of the variables as random
  virtual Real delta_covariance(const RealVector& x,
				PolynomialApproximation* poly_approx_2);

  /// return the change in mean between two response expansions,
  /// treating all variables as random
  virtual Real delta_mean();
  /// return the change in mean between two response expansions,
  /// treating a subset of the variables as random
  virtual Real delta_mean(const RealVector& x);
  /// return the change in standard deviation between two response
  /// expansions, treating all variables as random
  virtual Real delta_std_deviation();
  /// return the change in standard deviation between two response
  /// expansions, treating a subset of the variables as random
  virtual Real delta_std_deviation(const RealVector& x);
  /// return the change in reliability index (mapped from z_bar) between
  /// two response expansions, treating all variables as random
  virtual Real delta_beta(bool cdf_flag, Real z_bar);
  /// return the change in reliability index (mapped from z_bar) between
  /// two response expansions, treating a subset of the variables as random
  virtual Real delta_beta(const RealVector& x, bool cdf_flag, Real z_bar);
  /// return the change in response level (mapped from beta_bar) between
  /// two response expansions, treating all variables as random
  virtual Real delta_z(bool cdf_flag, Real beta_bar);
  /// return the change in response level (mapped from beta_bar) between
  /// two response expansions, treating a subset of the variables as random
  virtual Real delta_z(const RealVector& x, bool cdf_flag, Real beta_bar);

  /// compute central response moments using some combination of expansion
  /// post-processing and numerical integration
  virtual void compute_moments() = 0;
  /// compute central response moments in all-variables mode using some
  /// combination of expansion post-processing and numerical integration
  virtual void compute_moments(const RealVector& x) = 0;

  //
  //- Heading: Member functions
  //

  /// set surrData (shared representation)
  void surrogate_data(const SurrogateData& data);
  /// get surrData
  const SurrogateData& surrogate_data() const;

  /// return expansionMoments
  const RealVector& expansion_moments() const;
  /// return numericalMoments
  const RealVector& numerical_integration_moments() const;
  /// return preferred response moments (either expansion or numerical
  /// integration, depending on approximation type)
  const RealVector& moments() const;

  /// standardize central moments 2-n and eliminate excess kurtosis
  void standardize_moments(const RealVector& central_moments,
			   RealVector& std_moments);

  // number of data points to remove in a decrement (implemented at this
  // intermediate level since surrData not defined at base level)
  //size_t pop_count();

  /// set ExpansionConfigOptions::expansionCoeffFlag
  void expansion_coefficient_flag(bool coeff_flag);
  /// get ExpansionConfigOptions::expansionCoeffFlag
  bool expansion_coefficient_flag() const;

  /// set ExpansionConfigOptions::expansionCoeffGradFlag
  void expansion_coefficient_gradient_flag(bool grad_flag);
  /// get ExpansionConfigOptions::expansionCoeffGradFlag
  bool expansion_coefficient_gradient_flag() const;

  /// return sobolIndices
  const RealVector& sobol_indices() const;
  /// return totalSobolIndices
  const RealVector& total_sobol_indices() const;

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// generic base class function mapped to gradient_basis_variables(x)
  const RealVector& gradient(const RealVector& x);
  /// generic base class function mapped to hessian_basis_variables(x)
  const RealSymMatrix& hessian(const RealVector& x);

  //
  //- Heading: Member functions
  //

  /// compute central moments of response using type1 numerical integration
  void integrate_moments(const RealVector& coeffs, const RealVector& t1_wts,
			 RealVector& moments);
  /// compute central moments of response using type1/2 numerical integration
  void integrate_moments(const RealVector& t1_coeffs,
			 const RealMatrix& t2_coeffs, const RealVector& t1_wts,
			 const RealMatrix& t2_wts, RealVector& moments);

  /// size total Sobol arrays
  void allocate_total_sobol();

  //
  //- Heading: Data
  //

  /// instance containing the variables (shared) and response (unique) data
  /// arrays for constructing a surrogate of a single response function
  SurrogateData surrData;

  /// flag for calculation of expansion coefficients from response values
  bool expansionCoeffFlag;
  /// flag for calculation of gradients of expansion coefficients from
  /// response gradients
  bool expansionCoeffGradFlag;

  /// mean and central moments 2/3/4 computed from the stochastic expansion
  /// form.  For OrthogPolyApproximation, these are primary, and for
  /// InterpPolyApproximation, they are secondary.  Conversions to standardized
  /// moments (std deviation, skewness, kurtosis) are performed elsewhere.
  RealVector expansionMoments;
  /// mean and central moments 2/3/4 computed via numerical integration of the
  /// response.  For OrthogPolyApproximation, these are secondary, and for
  /// InterpPolyApproximation, they are primary.  Conversions to standardized
  /// moments (std deviation, skewness, kurtosis) are performed elsewhere.
  RealVector numericalMoments;

  /// gradient of the polynomial approximation returned by gradient()
  RealVector approxGradient;
  /// gradient of the polynomial approximation returned by gradient()
  RealSymMatrix approxHessian;
  /// gradient of the primary mean (expansion mean for OrthogPoly,
  /// numerical integration mean for InterpPoly)
  RealVector meanGradient;
  /// gradient of the primary variance (expansion variance for OrthogPoly,
  /// numerical integration variance for InterpPoly)
  RealVector varianceGradient;

  /// track computation of mean and mean gradient to avoid unnecessary
  /// recomputation
  short computedMean;
  /// track computation of variance and variance gradient to avoid
  /// unnecessary recomputation
  short computedVariance;
  /// track previous evaluation point for all_variables mean to avoid
  /// unnecessary recomputation
  RealVector xPrevMean;
  /// track previous evaluation point for all_variables mean gradient
  /// to avoid unnecessary recomputation
  RealVector xPrevMeanGrad;
  /// track previous evaluation point for all_variables variance to
  /// avoid unnecessary recomputation
  RealVector xPrevVar;
  /// track previous evaluation point for all_variables variance
  /// gradient to avoid unnecessary recomputation
  RealVector xPrevVarGrad;

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
  BasisApproximation(BaseConstructor(), shared_data), computedMean(0),
  computedVariance(0), expansionCoeffFlag(true), expansionCoeffGradFlag(false)
{ }


inline PolynomialApproximation::~PolynomialApproximation()
{ }


inline const SurrogateData& PolynomialApproximation::surrogate_data() const
{ return surrData; }


inline void PolynomialApproximation::surrogate_data(const SurrogateData& data)
{ surrData = data; /* shared representation */ }


inline const RealVector& PolynomialApproximation::expansion_moments() const
{ return expansionMoments; }


inline const RealVector& PolynomialApproximation::
numerical_integration_moments() const
{ return numericalMoments; }


/** All current cases prefer expansionMoments (with the distinction
    drawn between interpolation and integration rules, we prefer Gauss
    integration rules on interpolant expansions). */
inline const RealVector& PolynomialApproximation::moments() const
{
  // TO DO: return to this (will require activation of exp moments for all view)
  //return expansionMoments;

  // TO DO: remove this temp hack
  // (avoids returning to derived class specialization)
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  short exp_type = data_rep->expConfigOptions.expBasisType;
  return (exp_type == NODAL_INTERPOLANT || exp_type == HIERARCHICAL_INTERPOLANT)
    ? numericalMoments : expansionMoments;
}


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
