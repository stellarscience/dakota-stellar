/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        KrawtchoukOrthogPolynomial
//- Description:  Class for Krawtchouk Orthogonal Polynomial
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#ifndef KRAWTCHOUK_ORTHOG_POLYNOMIAL_HPP
#define KRAWTCHOUK_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Krawtchouk polynomials

/** The KrawtchoukOrthogPolynomial class evaluates a univariate Krawtchouk
    polynomial K^(p,N)_n of a particular order.  These polynomials
    are orthogonal with respect to the weight function 
    (N choose k)*p^k*(1-p)^(n-k).
    This corresponds to the binomial probability mass function,  
    which is the probability of k successes from N trials.
    See appendix in Xiu & Karniadakis, Siam J. Sci. Comp., v24, n2,
    pp. 619-644, 2002 for more details.  */

class KrawtchoukOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  KrawtchoukOrthogPolynomial();
  /// destructor
  ~KrawtchoukOrthogPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  //
  //- Heading: Noninherited memeber functions
  //

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

  /// naming convention consistent with Binomial distribution

  /// the probability of a "success" for each experiment
  Real probPerTrial;
  /// the number of discrete points on which to base the polynomial 
  unsigned int numTrials;
};


inline KrawtchoukOrthogPolynomial::KrawtchoukOrthogPolynomial() :
  probPerTrial(0.), numTrials(0) // dummy values prior to update
{ }


inline KrawtchoukOrthogPolynomial::~KrawtchoukOrthogPolynomial()
{ }


inline void KrawtchoukOrthogPolynomial::
pull_parameter(short dist_param, Real& param) const
{
  switch (dist_param) {
  case BI_P_PER_TRIAL: param = probPerTrial;    break;
  default:
    PCerr << "Error: unsupported distribution parameter in KrawtchoukOrthog"
	  << "Polynomial::pull_parameter(Real)." << std::endl;
    abort_handler(-1);
  }
}


inline void KrawtchoukOrthogPolynomial::
pull_parameter(short dist_param, unsigned int& param) const
{
  switch (dist_param) {
  case BI_TRIALS: param = numTrials; break;
  default:
    PCerr << "Error: unsupported distribution parameter in KrawtchoukOrthog"
	  << "Polynomial::pull_parameter(unsigned int)." << std::endl;
    abort_handler(-1);
  }
}


inline void KrawtchoukOrthogPolynomial::
push_parameter(short dist_param, Real param)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPointsMap.empty() || collocWeightsMap.empty()) // first pass
    switch (dist_param) {
    case BI_P_PER_TRIAL: probPerTrial = param;   break;
    }
  else
    switch (dist_param) {
    case BI_P_PER_TRIAL:
      if (!real_compare(probPerTrial, param))
	{ probPerTrial = param;  reset_gauss(); }
      break;
    }
}


inline void KrawtchoukOrthogPolynomial::
push_parameter(short dist_param, unsigned int param)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPointsMap.empty() || collocWeightsMap.empty()) // first pass
    switch (dist_param) {
    case BI_TRIALS: numTrials = param; break;
    }
  else
    switch (dist_param) {
    case BI_TRIALS:
      if (numTrials != param)  { numTrials = param;  reset_gauss(); }
      break;
    }
}


inline bool KrawtchoukOrthogPolynomial::parameterized() const
{ return true; }

} // namespace Pecos

#endif
