/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HermiteInterpPolynomial
//- Description:  Class for 1-D Hermite Interpolation Polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef HERMITE_INTERP_POLYNOMIAL_HPP
#define HERMITE_INTERP_POLYNOMIAL_HPP

#include "InterpolationPolynomial.hpp"
#include "pecos_data_types.hpp"


namespace Pecos {

/// Derived basis polynomial class for 1-D Hermite interpolation polynomials

/** The HermiteInterpPolynomial class evaluates a univariate Hermite
    interpolation polynomial.  The order of the polynomial is dictated
    by the number of interpolation points (order = N_p - 1).  It enables
    multidimensional interpolants within InterpPolyApproximation. */

class HermiteInterpPolynomial: public InterpolationPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  HermiteInterpPolynomial();
  /// constructor with collocation rule
  HermiteInterpPolynomial(short colloc_rule);
  /// constructor with set of interpolation points
  HermiteInterpPolynomial(const RealArray& interp_pts);
  /// destructor
  ~HermiteInterpPolynomial();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void precompute_data();

  Real type1_value(Real x, unsigned short i);
  Real type2_value(Real x, unsigned short i);

  Real type1_gradient(Real x, unsigned short i);
  Real type2_gradient(Real x, unsigned short i);

  const RealArray& collocation_points(unsigned short order);
  const RealArray& type1_collocation_weights(unsigned short order);
  const RealArray& type2_collocation_weights(unsigned short order);

  void collocation_rule(short rule);
  short collocation_rule() const;

private:

  //
  //- Heading: Data
  //

  /// name of uniform collocation rule: GAUSS_PATTERSON,
  /// CLENSHAW_CURTIS, FEJER2, or GAUSS_LEGENDRE
  short collocRule;

  /// set of 1-D weights for interpolation of values
  RealArray type1InterpWts;
  /// set of 1-D] weights for interpolation of gradients
  RealArray type2InterpWts;

  /// pre-computed divided difference table for input used in type1/2
  /// value calculation
  RealArray xValDiffTab;
  /// pre-computed divided difference table for input used in type1/2
  /// gradient calculation
  RealArray xGradDiffTab;

  /// pre-computed divided difference table for output used in type1
  /// value calculation
  Real2DArray yT1ValDiffTab;
  /// pre-computed divided difference table for output used in type1
  /// gradient calculation
  Real2DArray yT1GradDiffTab;

  /// pre-computed divided difference table for output used in type2
  /// value calculation
  Real2DArray yT2ValDiffTab;
  /// pre-computed divided difference table for output used in type2
  /// gradient calculation
  Real2DArray yT2GradDiffTab;
};


inline HermiteInterpPolynomial::HermiteInterpPolynomial():
  InterpolationPolynomial(), collocRule(GAUSS_LEGENDRE)
{ wtFactor = 0.5; }


inline HermiteInterpPolynomial::HermiteInterpPolynomial(short colloc_rule):
  InterpolationPolynomial(), collocRule(colloc_rule)
{ wtFactor = 0.5; }


inline HermiteInterpPolynomial::
HermiteInterpPolynomial(const RealArray& interp_pts):
  InterpolationPolynomial(interp_pts), collocRule(GAUSS_LEGENDRE)
{ wtFactor = 0.5; }


inline HermiteInterpPolynomial::~HermiteInterpPolynomial()
{ }


inline void HermiteInterpPolynomial::collocation_rule(short rule)
{ collocRule = rule; }


inline short HermiteInterpPolynomial::collocation_rule() const
{ return collocRule; }

} // namespace Pecos

#endif
