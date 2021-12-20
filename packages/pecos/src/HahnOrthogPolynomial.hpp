/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HahnOrthogPolynomial
//- Description:  Class for Hahn Orthogonal Polynomial
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#ifndef HAHN_ORTHOG_POLYNOMIAL_HPP
#define HAHN_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Hahn polynomials

/** The HahnOrthogPolynomial class evaluates a univariate Hahn
    polynomial Q^(alpha,beta,N)_n of a particular order.  These polynomials
    are orthogonal with respect to the weight function 

    (K choose k)(N-K choose n-k)/( N choose n).

    This corresponds to the hypergeometric probability mass function, which
    describes the probability of k successes in n draws, without 
    replacement, from a finite population of size N that contains exactly 
    K successes.
    See appendix in Xiu & Karniadakis, Siam J. Sci. Comp., v24, n2,
    pp. 619-644, 2002 for more details.  */

class HahnOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  HahnOrthogPolynomial();
  /// destructor
  ~HahnOrthogPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //


protected:

  //
  //- Heading: Virtual function redefinitions
  //
  Real type1_value(Real x, unsigned short order);

  void pull_parameter(short dist_param, unsigned int& param) const;
  void push_parameter(short dist_param, unsigned int  param);
  bool parameterized() const;

private:

  //
  //- Heading: Data
  //

  /// number in total population
  unsigned int totalPop;
  /// number in selected population
  unsigned int selectPop;
  /// number drawn from population
  unsigned int numDrawn;
};


inline HahnOrthogPolynomial::HahnOrthogPolynomial() :
  totalPop(0), selectPop(0), numDrawn(0) // dummy values prior to update
{ }


inline HahnOrthogPolynomial::~HahnOrthogPolynomial()
{ }


inline void HahnOrthogPolynomial::
pull_parameter(short dist_param, unsigned int& param) const
{
  switch (dist_param) {
  case HGE_TOT_POP: param = totalPop;  break;
  case HGE_SEL_POP: param = selectPop; break;
  case HGE_DRAWN:   param = numDrawn;  break;
  default:
    PCerr << "Error: unsupported distribution parameter in HahnOrthogPolynomial"
	  << "::parameter()." << std::endl;
    abort_handler(-1);
  }
}


inline void HahnOrthogPolynomial::
push_parameter(short dist_param, unsigned int param)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPointsMap.empty() || collocWeightsMap.empty()) { // first pass
    switch (dist_param) {
    case HGE_TOT_POP:  totalPop = param; break;
    case HGE_SEL_POP: selectPop = param; break;
    case HGE_DRAWN:    numDrawn = param; break;
    }
  }
  else {
    switch (dist_param) {
    case HGE_TOT_POP:
      if (totalPop  != param)  {  totalPop = param;  reset_gauss(); }
      break;
    case HGE_SEL_POP:
      if (selectPop != param)  { selectPop = param;  reset_gauss(); }
      break;
    case HGE_DRAWN:
      if (numDrawn  != param)  {  numDrawn = param;  reset_gauss(); }
      break;
    }
  }
}


inline bool HahnOrthogPolynomial::parameterized() const
{ return true; }

} // namespace Pecos

#endif
