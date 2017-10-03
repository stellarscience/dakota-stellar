/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PROBABILITY_TRANSFORMATION_HPP
#define PROBABILITY_TRANSFORMATION_HPP

#include "pecos_data_types.hpp"
#include "NormalRandomVariable.hpp"
#include "UniformRandomVariable.hpp"
#include "ExponentialRandomVariable.hpp"

namespace Pecos {

class AleatoryDistParams;
class EpistemicDistParams;


/// Base class for all nonlinear distribution transformations

/** The base class for nonlinear distribution transformations,
    including Nataf, Rosenblatt, et al. */

class ProbabilityTransformation
{
public:

  /// default constructor
  ProbabilityTransformation();
  /// standard constructor for envelope
  ProbabilityTransformation(const String& prob_trans_type);
  /// copy constructor
  ProbabilityTransformation(const ProbabilityTransformation& prob_trans);

  /// destructor
  virtual ~ProbabilityTransformation();

  /// assignment operator
  ProbabilityTransformation
    operator=(const ProbabilityTransformation& prob_trans);

  //
  //- Heading: Virtual functions
  //

  /// Transformation routine from u-space of uncorrelated standard normal
  /// variables to x-space of correlated random variables
  virtual void trans_U_to_X(const RealVector& u_vars, RealVector& x_vars);

  /// Transformation routine from x-space of correlated random variables 
  /// to u-space of uncorrelated standard normal variables
  virtual void trans_X_to_U(const RealVector& x_vars, RealVector& u_vars);

  /// As part of the Nataf distribution model (Der Kiureghian & Liu, 1986),
  /// this procedure modifies the user-specified correlation matrix
  /// (corrMatrixX) to account for correlation warping from the nonlinear
  /// X->Z transformation and performs a Cholesky factorization to create
  /// corrCholeskyFactorZ.
  virtual void transform_correlations();

  /// Transformation routine for gradient vector from x-space to u-space
  virtual void trans_grad_X_to_U(const RealVector& fn_grad_x,
				 RealVector& fn_grad_u,
				 const RealVector& x_vars,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);
  /// Transformation routine for gradient vector from x-space to u-space
  virtual void trans_grad_X_to_U(const RealVector& fn_grad_x,
				 RealVector& fn_grad_u,
				 const RealMatrix& jacobian_xu,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);

  /// Transformation routine from x-space gradient vector to design space
  virtual void trans_grad_X_to_S(const RealVector& fn_grad_x,
				 RealVector& fn_grad_s,
				 const RealVector& x_vars,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids,
				 SizetMultiArrayConstView acv_ids,
				 const SizetArray& acv_map1_indices,
				 const ShortArray& acv_map2_targets);
  /// Transformation routine from x-space gradient vector to design space
  virtual void trans_grad_X_to_S(const RealVector& fn_grad_x,
				 RealVector& fn_grad_s,
				 const RealMatrix& jacobian_xs,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids,
				 SizetMultiArrayConstView acv_ids,
				 const SizetArray& acv_map1_indices,
				 const ShortArray& acv_map2_targets);

  /// Transformation routine for gradient vector from u-space to x-space
  virtual void trans_grad_U_to_X(const RealVector& fn_grad_u,
				 RealVector& fn_grad_x,
				 const RealVector& x_vars,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);
  /// Transformation routine for gradient vector from u-space to x-space
  virtual void trans_grad_U_to_X(const RealVector& fn_grad_u,
				 RealVector& fn_grad_x,
				 const RealMatrix& jacobian_ux,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);

  /// Transformation routine for Hessian matrix from x-space to u-space
  virtual void trans_hess_X_to_U(const RealSymMatrix& fn_hess_x,
				 RealSymMatrix& fn_hess_u,
				 const RealVector& x_vars,
				 const RealVector& fn_grad_x,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);
  /// Transformation routine for Hessian matrix from x-space to u-space
  virtual void trans_hess_X_to_U(const RealSymMatrix& fn_hess_x,
				 RealSymMatrix& fn_hess_u,
				 const RealMatrix& jacobian_xu,
				 const RealSymMatrixArray& hessian_xu,
				 const RealVector& fn_grad_x,
				 const SizetArray& x_dvv,
				 SizetMultiArrayConstView cv_ids);

  /// Jacobian of x(u) mapping obtained from dX/dZ dZ/dU
  virtual void jacobian_dX_dU(const RealVector& x_vars,
			      RealMatrix& jacobian_xu);

  /// Jacobian of u(x) mapping obtained from dU/dZ dZ/dX
  virtual void jacobian_dU_dX(const RealVector& x_vars,
			      RealMatrix& jacobian_ux);

  /// Design Jacobian of x(u,s) mapping obtained from differentiation of
  /// trans_U_to_X() with respect to distribution parameters S
  virtual void jacobian_dX_dS(const RealVector& x_vars, RealMatrix& jacobian_xs,
			      SizetMultiArrayConstView cv_ids,
			      SizetMultiArrayConstView acv_ids,
			      const SizetArray& acv_map1_indices,
			      const ShortArray& acv_map2_targets);

  /// Hessian of x(u) mapping obtained from dZ/dU^T d^2X/dZ^2 dZ/dU
  virtual void hessian_d2X_dU2(const RealVector& x_vars,
			       RealSymMatrixArray& hessian_xu);

  //
  //- Heading: Member functions
  //

  /// perform a deep copy of incoming prob_trans
  void copy(const ProbabilityTransformation& prob_trans);

  /// initializes randomVarsX (no transformation: u-space not needed)
  void initialize_random_variable_types(const ShortArray& x_types);
  /// initializes randomVarsX and ranVarTypesU
  void initialize_random_variable_types(const ShortArray& x_types,
					const ShortArray& u_types);
  /// updates parameters within randomVarsX
  void initialize_random_variable_parameters(const RealVector& cd_l_bnds,
					     const RealVector& cd_u_bnds,
					     const AleatoryDistParams& adp,
					     const EpistemicDistParams& edp,
					     const RealVector& cs_l_bnds,
					     const RealVector& cs_u_bnds);
  /// initializes corrMatrixX and correlationFlagX
  void initialize_random_variable_correlations(const RealSymMatrix& x_corr);

  /// reshape corrMatrixX for an all_variables specification
  void reshape_correlation_matrix(size_t num_leading_vars,
				  size_t num_probabilistic_vars,
				  size_t num_trailing_vars);

  /// return randomVarsX
  const std::vector<RandomVariable>& x_random_variables() const;

  /// assemble means and standard deviations in original space from
  /// RandomVariable::moments()
  RealRealPairArray x_moments() const;
  /// assemble means and standard deviations in transformed space
  /// combining data from ranVarTypesU with RandomVariable::moments()
  RealRealPairArray u_moments() const;
  /// assemble means in original space from RandomVariable::mean()
  RealVector x_means() const;
  /// assemble standard deviations in original space from
  /// RandomVariable::standard_deviation()
  RealVector x_std_deviations() const;

  /// assemble lower and upper bounds from RandomVariable::bounds()
  RealRealPairArray x_bounds() const;
  /// assemble lower and upper bounds in transformed space combining
  /// data from ranVarTypesU with RandomVariable::bounds()
  RealRealPairArray u_bounds() const;
  /// assemble lower bounds from RandomVariable::bounds()
  RealVector x_lower_bounds() const;
  /// assemble upper bounds from RandomVariable::bounds()
  RealVector x_upper_bounds() const;

  /// return the univariate PDF value for an x-space random variable
  Real x_pdf(Real x_val, size_t i) const;
  /// return the univariate log PDF value for an x-space random variable
  Real x_log_pdf(Real x_val, size_t i) const;
  /// return the gradient of the univariate log PDF for an x-space
  /// random variable
  Real x_log_pdf_gradient(Real x_val, size_t i) const;
  /// return the Hessian of the univariate log PDF for an x-space
  /// random variable
  Real x_log_pdf_hessian(Real x_val, size_t i) const;
  /// return the univariate PDF value for a u-space random variable
  Real u_pdf(Real u_val, size_t i) const;
  /// return the univariate log PDF value for a u-space random variable
  Real u_log_pdf(Real u_val, size_t i) const;
  /// return the gradient of the univariate log PDF for a u-space
  /// random variable
  Real u_log_pdf_gradient(Real u_val, size_t i) const;
  /// return the Hessian of the univariate log PDF for a u-space random variable
  Real u_log_pdf_hessian(Real u_val, size_t i) const;

  /// return the multivariate PDF value for x-space random variables
  Real x_pdf(const RealVector& x_pt) const;
  /// return the multivariate log PDF value for x-space random variables
  Real x_log_pdf(const RealVector& x_pt) const;
  /// return the multivariate PDF value for u-space random variables
  Real u_pdf(const RealVector& u_pt) const;
  /// return the multivariate log PDF value for u-space random variables
  Real u_log_pdf(const RealVector& u_pt) const;

  /// draw a sample from an x-space random variable
  template <typename Engine>
  Real draw_x_sample(size_t i, Engine& rng) const;
  /// draw a sample from a u-space random variable
  template <typename Engine>
  Real draw_u_sample(size_t i, Engine& rng) const;

  /// return ranVarTypesU
  const ShortArray& u_types() const;
  /// set ranVarTypesU
  void u_types(const ShortArray& types);
  /// set ranVarTypesU[i]
  void u_type(short type, size_t i);

  /// verify that randomVarsX[i].type() equals x_type
  void check_x_type(size_t i, short x_type) const;

  /// return correlationFlagX
  bool x_correlation() const;
  /// return corrMatrixX
  const RealSymMatrix& x_correlation_matrix() const;
  /// return corrCholeskyFactorZ
  const RealMatrix& z_correlation_factor() const;

  /// function to check modelRep (does this envelope contain a letter)
  bool is_null() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  ProbabilityTransformation(BaseConstructor);

  //
  //- Heading: Member functions
  //

  /// Computes numerical dx/ds and dz/ds Jacobians as requested by xs
  /// and zs booleans
  void numerical_design_jacobian(const RealVector& x_vars,
                                 bool xs, RealMatrix& num_jacobian_xs,
                                 bool zs, RealMatrix& num_jacobian_zs,
				 SizetMultiArrayConstView cv_ids,
				 SizetMultiArrayConstView acv_ids,
				 const SizetArray& acv_map1_indices,
				 const ShortArray& acv_map2_targets);

#ifdef DERIV_DEBUG
  /// routine for verification of transformation Jacobian/Hessian terms
  void verify_trans_jacobian_hessian(const RealVector& v0);

  /// routine for verification of design Jacobian terms
  void verify_design_jacobian(const RealVector& u0);
#endif // DERIV_DEBUG

  //
  //- Heading: Data members
  //

  /// vector of random variables encapsulating distribution parameters and
  /// statistical functions (pdf, cdf, etc.)
  std::vector<RandomVariable> randomVarsX;
  /// vector of types of each u-space standardized uncertain variable to
  /// which each x-space variable is transformed
  ShortArray ranVarTypesU;

  /// flag for indicating if correlation exists among the x-space
  /// uncertain variables
  bool correlationFlagX;
  /// matrix of random variable correlation coefficients
  RealSymMatrix corrMatrixX;
  /// cholesky factor of a modified correlation matrix (#corrMatrixX
  /// is modified in transform_correlations() for use in z-space)
  RealMatrix corrCholeskyFactorZ;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// probTransRep to the appropriate derived type.
  static ProbabilityTransformation*
    get_prob_trans(const String& prob_trans_type);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  ProbabilityTransformation* probTransRep;
  /// number of objects sharing probTransRep
  int referenceCount;
};


inline const std::vector<RandomVariable>& ProbabilityTransformation::
x_random_variables() const
{ return (probTransRep) ? probTransRep->randomVarsX : randomVarsX; }


inline void ProbabilityTransformation::
check_x_type(size_t i, short x_type) const
{
  if (probTransRep) probTransRep->check_x_type(i, x_type);
  else if (randomVarsX[i].type() != x_type) {
    PCerr << "Error: failure in x_type check in ProbabilityTransformation."
	  << std::endl;
    abort_handler(-1);
  }
}


inline RealRealPairArray ProbabilityTransformation::x_moments() const
{
  if (probTransRep) return probTransRep->x_moments();
  else {
    size_t i, num_v = randomVarsX.size();
    RealRealPairArray x_mom(num_v);
    for (i=0; i<num_v; ++i)
      x_mom[i] = randomVarsX[i].moments();
    return x_mom;
  }
}


inline RealVector ProbabilityTransformation::x_means() const
{
  if (probTransRep) return probTransRep->x_means();
  else {
    size_t i, num_v = randomVarsX.size();
    RealVector means(num_v, false);
    for (i=0; i<num_v; ++i)
      means[i] = randomVarsX[i].mean();
    return means;
  }
}


inline RealVector ProbabilityTransformation::x_std_deviations() const
{
  if (probTransRep) return probTransRep->x_std_deviations();
  else {
    size_t i, num_v = randomVarsX.size();
    RealVector std_devs(num_v, false);
    for (i=0; i<num_v; ++i)
      std_devs[i] = randomVarsX[i].standard_deviation();
    return std_devs;
  }
}


inline RealRealPairArray ProbabilityTransformation::x_bounds() const
{
  if (probTransRep) return probTransRep->x_bounds();
  else {
    size_t i, num_v = randomVarsX.size();
    RealRealPairArray x_bnds(num_v);
    for (i=0; i<num_v; ++i)
      x_bnds[i] = randomVarsX[i].bounds();
    return x_bnds;
  }
}


inline RealVector ProbabilityTransformation::x_lower_bounds() const
{
  if (probTransRep) return probTransRep->x_lower_bounds();
  else {
    size_t i, num_v = randomVarsX.size();
    RealVector lwr_bnds(num_v, false);
    for (i=0; i<num_v; ++i)
      lwr_bnds[i] = randomVarsX[i].bounds().first;
    return lwr_bnds;
  }
}


inline RealVector ProbabilityTransformation::x_upper_bounds() const
{
  if (probTransRep) return probTransRep->x_upper_bounds();
  else {
    size_t i, num_v = randomVarsX.size();
    RealVector upr_bnds(num_v, false);
    for (i=0; i<num_v; ++i)
      upr_bnds[i] = randomVarsX[i].bounds().second;
    return upr_bnds;
  }
}


inline Real ProbabilityTransformation::x_pdf(Real x_val, size_t i) const
{
  return (probTransRep) ? probTransRep->randomVarsX[i].pdf(x_val) :
    randomVarsX[i].pdf(x_val);
}


inline Real ProbabilityTransformation::x_log_pdf(Real x_val, size_t i) const
{
  return (probTransRep) ? probTransRep->randomVarsX[i].log_pdf(x_val) :
    randomVarsX[i].log_pdf(x_val);
}


inline Real ProbabilityTransformation::
x_log_pdf_gradient(Real x_val, size_t i) const
{
  return (probTransRep) ? probTransRep->randomVarsX[i].log_pdf_gradient(x_val) :
    randomVarsX[i].log_pdf_gradient(x_val);
}


inline Real ProbabilityTransformation::
x_log_pdf_hessian(Real x_val, size_t i) const
{
  return (probTransRep) ? probTransRep->randomVarsX[i].log_pdf_hessian(x_val) :
    randomVarsX[i].log_pdf_hessian(x_val);
}


inline Real ProbabilityTransformation::x_pdf(const RealVector& x_pt) const
{
  if (probTransRep) return probTransRep->x_pdf(x_pt);
  else {
    // TO DO: add support for evaluation of correlated MVN density
    if (correlationFlagX) {
      PCerr << "Error: ProbabilityTransformation::x_pdf() currently uses a "
	    << "product of marginal densities\n       and can only be used for "
	    << "independent random variables." << std::endl;
      abort_handler(-1);
    }
    size_t i, num_v = randomVarsX.size();
    Real density = 1.;
    for (i=0; i<num_v; ++i)
      density *= x_pdf(x_pt[i], i);
    return density;
  }
}


inline Real ProbabilityTransformation::x_log_pdf(const RealVector& x_pt) const
{
  if (probTransRep) return probTransRep->x_log_pdf(x_pt);
  else {
    // TO DO: add support for evaluation of correlated MVN density
    if (correlationFlagX) {
      PCerr << "Error: ProbabilityTransformation::x_log_pdf() currently uses a "
	    << "sum of log marginal densities\n       and can only be used for "
	    << "independent random variables." << std::endl;
      abort_handler(-1);
    }
    size_t i, num_v = randomVarsX.size();
    Real log_density = 0.;
    for (i=0; i<num_v; ++i)
      log_density += x_log_pdf(x_pt[i], i);
    return log_density;
  }
}


inline Real ProbabilityTransformation::u_pdf(const RealVector& u_pt) const
{
  if (probTransRep) return probTransRep->u_pdf(u_pt);
  else {
    // u-space is independent -> use product of marginals
    size_t i, num_v = randomVarsX.size();
    Real density = 1.;
    for (i=0; i<num_v; ++i)
      density *= u_pdf(u_pt[i], i);
    return density;
  }
}


inline Real ProbabilityTransformation::u_log_pdf(const RealVector& u_pt) const
{
  if (probTransRep) return probTransRep->u_log_pdf(u_pt);
  else {
    // u-space is independent -> use sum of log marginals
    size_t i, num_v = randomVarsX.size();
    Real log_density = 0.;
    for (i=0; i<num_v; ++i)
      log_density += u_log_pdf(u_pt[i], i);
    return log_density;
  }
}


template <typename Engine> 
Real ProbabilityTransformation::draw_x_sample(size_t i, Engine& rng) const
{
  return (probTransRep) ? probTransRep->randomVarsX[i].draw_sample(rng) :
    randomVarsX[i].draw_sample(rng);
}


template <typename Engine> 
Real ProbabilityTransformation::draw_u_sample(size_t i, Engine& rng) const
{
  if (probTransRep) return probTransRep->draw_u_sample(i, rng);
  else {
    // can only use randomVarsX[i].standard_pdf() for cases where u_type is a
    // standardized form of the x_type.  For STD_NORMAL and STD_UNIFORM, many
    // x_types can be mapped to these u_types, so use global utility fns
    // whenever there are no auxilliary parameters to manage.
    switch (ranVarTypesU[i]) {
    // these cases require static fns since U type may not correspond to X type
    case STD_NORMAL:  return  NormalRandomVariable::draw_std_sample(rng); break;
    case STD_UNIFORM: return UniformRandomVariable::draw_std_sample(rng); break;
    case STD_EXPONENTIAL:
      return ExponentialRandomVariable::draw_std_sample(rng);             break;
    // these cases can rely on correspondence between X and U types
    case STD_BETA:
      check_x_type(i, BETA);
      return randomVarsX[i].draw_standard_sample(rng); break;
    case STD_GAMMA:
      check_x_type(i, GAMMA);
      return randomVarsX[i].draw_standard_sample(rng); break;
    default: // no transformation (e.g., PCE with numerically-generated bases)
      check_x_type(i, ranVarTypesU[i]);
      return randomVarsX[i].draw_sample(rng);          break;
    }
  }
}


inline const ShortArray& ProbabilityTransformation::u_types() const
{ return (probTransRep) ? probTransRep->ranVarTypesU : ranVarTypesU; }


inline void ProbabilityTransformation::u_types(const ShortArray& types)
{
  if (probTransRep) probTransRep->ranVarTypesU = types;
  else              ranVarTypesU = types;
}


inline void ProbabilityTransformation::u_type(short type, size_t i)
{
  if (probTransRep) probTransRep->ranVarTypesU[i] = type;
  else              ranVarTypesU[i] = type;
}


inline bool ProbabilityTransformation::x_correlation() const
{ return (probTransRep) ? probTransRep->correlationFlagX : correlationFlagX; }


inline const RealSymMatrix& ProbabilityTransformation::
x_correlation_matrix() const
{ return (probTransRep) ? probTransRep->corrMatrixX : corrMatrixX; }


inline const RealMatrix& ProbabilityTransformation::z_correlation_factor() const
{
  return (probTransRep) ? probTransRep->corrCholeskyFactorZ :
                          corrCholeskyFactorZ;
}


inline bool ProbabilityTransformation::is_null() const
{ return (probTransRep) ? false : true; }

} // namespace Pecos

#endif
