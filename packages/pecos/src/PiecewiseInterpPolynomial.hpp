/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        PiecewiseInterpPolynomial
//- Description:  Class for 1-D piecewise Interpolation Polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef PIECEWISE_INTERP_POLYNOMIAL_HPP
#define PIECEWISE_INTERP_POLYNOMIAL_HPP

#include "InterpolationPolynomial.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {


/// Derived basis polynomial class for 1-D piecewise interpolation polynomials

/** The PiecewiseInterpPolynomial class evaluates a univariate
    interpolation polynomial with local support.  The order of the
    polynomial may be linear, based only on interpolated values, or
    cubic, based on interpolated values and gradients.  It enables
    multidimensional interpolants within InterpPolyApproximation. */

class PiecewiseInterpPolynomial: public InterpolationPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  // default constructor
  PiecewiseInterpPolynomial();
  /// constructor with poly_type and rule arguments
  PiecewiseInterpPolynomial(short poly_type, short rule = NEWTON_COTES);
  /// constructor with poly_type, rule, and set of points to interpolate
  PiecewiseInterpPolynomial(const RealArray& interp_pts, short poly_type,
			    short rule = NEWTON_COTES);
  /// destructor
  ~PiecewiseInterpPolynomial();

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

  /// name of closed nested rule: NEWTON_COTES (equidistant) or
  /// CLENSHAW_CURTIS (non-equidistant)
  short collocRule;

  /// the constant interval between points for an equidistant collocRule
  Real interpInterval;

  /// set of 1-D weights for interpolation of values
  RealArray type1InterpWts;
  /// set of 1-D] weights for interpolation of gradients
  RealArray type2InterpWts;
};


inline PiecewiseInterpPolynomial::PiecewiseInterpPolynomial():
  InterpolationPolynomial(), collocRule(NEWTON_COTES)
{ basisPolyType = PIECEWISE_LINEAR_INTERP; }


inline PiecewiseInterpPolynomial::
PiecewiseInterpPolynomial(short poly_type, short rule):
  InterpolationPolynomial(), collocRule(rule)
{ basisPolyType = poly_type; }


inline PiecewiseInterpPolynomial::
PiecewiseInterpPolynomial(const RealArray& interp_pts, short poly_type,
			  short rule):
  InterpolationPolynomial(interp_pts), collocRule(rule)
{ basisPolyType = poly_type; }


inline PiecewiseInterpPolynomial::~PiecewiseInterpPolynomial()
{ }


inline void PiecewiseInterpPolynomial::collocation_rule(short rule)
{ collocRule = rule; }


inline short PiecewiseInterpPolynomial::collocation_rule() const
{ return collocRule; }

} // namespace Pecos

#endif
