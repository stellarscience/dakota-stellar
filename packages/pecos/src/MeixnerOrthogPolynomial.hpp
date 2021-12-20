/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        MeixnerOrthogPolynomial
//- Description:  Class for Meixner Orthogonal Polynomial
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#ifndef MEIXNER_ORTHOG_POLYNOMIAL_HPP
#define MEIXNER_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Meixner polynomials

/** The MeixnerOrthogPolynomial class evaluates a univariate Meixner
    polynomial M^(c,Beta)_n of a particular order.  These polynomials
    are orthogonal with respect to the weight function 

    (k+r-1 choose k)*p^k*(1-p)^n, where p is the the probability of success
    of a trial, n is the number of total successes and k is the trial number.
    The geometric distribution (n==1) is a special case of the negative 
    binomial distribution.
    This corresponds to the negative binomial probability mass function.
    See appendix in Xiu & Karniadakis, Siam J. Sci. Comp., v24, n2,
    pp. 619-644, 2002 for more details.  */

class MeixnerOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  MeixnerOrthogPolynomial();
  /// destructor
  ~MeixnerOrthogPolynomial();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  Real type1_value(Real x, unsigned short order);

  void pull_parameter(short dist_param, Real& param) const;
  void pull_parameter(short dist_param, unsigned int& param) const;
  void push_parameter(short dist_param, Real  param);
  void push_parameter(short dist_param, unsigned int  param);
  bool parameterized() const;

private:

  //
  //- Heading: Data
  //

  /// the probability of a "success" for each experiment
  Real probPerTrial;
  /// the number of discrete points on which to base the polynomial 
  unsigned int numTrials;
};


inline MeixnerOrthogPolynomial::MeixnerOrthogPolynomial() :
  probPerTrial(0.), // dummy value prior to update
  numTrials(1) // default for Geometric dist (overridden for Negative Binomial)
{ }


inline MeixnerOrthogPolynomial::~MeixnerOrthogPolynomial()
{ }


inline void MeixnerOrthogPolynomial::
pull_parameter(short dist_param, Real& param) const
{
  switch (dist_param) {
  case NBI_P_PER_TRIAL: case GE_P_PER_TRIAL: param = probPerTrial;    break;
  default:
    PCerr << "Error: unsupported distribution parameter in MeixnerOrthog"
	  << "Polynomial::pull_parameter(Real)." << std::endl;
    abort_handler(-1);
  }
}


inline void MeixnerOrthogPolynomial::
pull_parameter(short dist_param, unsigned int& param) const
{
  switch (dist_param) {
  case NBI_TRIALS: param = numTrials; break;
  default:
    PCerr << "Error: unsupported distribution parameter in MeixnerOrthog"
	  << "Polynomial::pull_parameter(unsigned int)." << std::endl;
    abort_handler(-1);
  }
}


inline void MeixnerOrthogPolynomial::
push_parameter(short dist_param, Real param)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPointsMap.empty() || collocWeightsMap.empty()) // first pass
    switch (dist_param) {
    case NBI_P_PER_TRIAL: case GE_P_PER_TRIAL: probPerTrial = param; break;
    }
  else
    switch (dist_param) {
    case NBI_P_PER_TRIAL: case GE_P_PER_TRIAL:
      if (!real_compare(probPerTrial, param))
	{ probPerTrial = param;  reset_gauss(); }
      break;
    }
}


inline void MeixnerOrthogPolynomial::
push_parameter(short dist_param, unsigned int param)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPointsMap.empty() || collocWeightsMap.empty()) // first pass
    switch (dist_param) {
    case NBI_TRIALS: /* case GE_TRIALS: */ numTrials = param; break;
    }
  else
    switch (dist_param) {
    case NBI_TRIALS: /* case GE_TRIALS: */
      if (numTrials != param)  { numTrials = param;  reset_gauss(); }
      break;
    }
}


inline bool MeixnerOrthogPolynomial::parameterized() const
{ return true; }

} // namespace Pecos

#endif
