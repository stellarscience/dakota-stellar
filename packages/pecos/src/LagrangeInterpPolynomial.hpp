/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LagrangeInterpPolynomial
//- Description:  Class for 1-D Lagrange Interpolation Polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef LAGRANGE_INTERP_POLYNOMIAL_HPP
#define LAGRANGE_INTERP_POLYNOMIAL_HPP

#include "InterpolationPolynomial.hpp"
#include "pecos_data_types.hpp"


namespace Pecos {

/// Derived basis polynomial class for 1-D Lagrange interpolation polynomials

/** The LagrangeInterpPolynomial class evaluates a univariate Lagrange
    interpolation polynomial.  The order of the polynomial is dictated
    by the number of interpolation points (order = N_p - 1).  It
    enables multidimensional interpolants within
    InterpPolyApproximation.  This class supports both the traditional
    characteristic polynomial form of Lagrange interpolation as well
    as barycentric Lagrange interpolation (the second form from Berrut
    and Trefethen, 2004).  The former is used for actual evaluation of
    1D polynomial values (when needed), whereas the latter allows
    alternative interpolant evaluations with additional precomputation
    that improve efficiency from O(n^2) evaluations to O(n). */

class LagrangeInterpPolynomial: public InterpolationPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  LagrangeInterpPolynomial();
  /// standard constructor
  LagrangeInterpPolynomial(const RealArray& interp_pts);
  /// destructor
  ~LagrangeInterpPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  /// retrieve the value of the i_th Lagrange polynomial for a given
  /// parameter x using barycentric formulation
  Real type1_value(unsigned short i);
  /// retrieve the value of the i_th Lagrange polynomial for a given
  /// parameter x using traditional characteristic polynomial formulation
  Real type1_value(Real x, unsigned short i);

  /// retrieve the gradient of the i_th Lagrange polynomial for a given
  /// parameter x using barycentric formulation
  Real type1_gradient(unsigned short i);
  /// retrieve the gradient of the i_th Lagrange polynomial for a given
  /// parameter x using traditional characteristic polynomial formulation
  Real type1_gradient(Real x, unsigned short i);

  void set_new_point(Real x, short request_order);
  void set_new_point(Real x, short request_order, const UShortArray& delta_key);

  size_t exact_index() const;
  size_t exact_delta_index() const;

  //const RealVector& barycentric_weights() const;
  const RealVector& barycentric_value_factors() const;
  const RealVector& barycentric_gradient_factors() const;

  Real barycentric_value_factor(unsigned short i) const;
  Real barycentric_gradient_factor(unsigned short i) const;

  Real barycentric_value_factor_sum() const;
  Real barycentric_difference_product() const;

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void precompute_data();

private:

  //
  //- Heading: Member functions
  //

  /// define compute order from request order and newPoint match
  void init_new_point(Real x, short request_order, short& compute_order);
  /// based on compute order, size barycentric value/gradient factors
  void allocate_factors(short compute_order);

  //
  //- Heading: Data
  //

  /// set of denominator products calculated from interpPts in
  /// precompute_data(); in barycentric formulations, these are the weights
  RealVector bcWeights;

  // additional data for barycentric formulation (second form)

  /// the parameter value for evaluation of the interpolant
  Real newPoint;
  /// order of data that has been precomputed at newPoint
  short newPtOrder;
  /// index of interpolation point that exactly matches the interpolated value x
  size_t exactIndex;
  /// index within a hierarchical increment to the interpolation points that
  /// exactly matches the interpolated value x
  size_t exactDeltaIndex;

  /// product of point differences (x-x_j) for a particular newPoint x
  /// and interpPts x_j
  Real diffProduct;

  /// terms bcWeights[j]/(x-x[j]) from barycentric formulation
  RealVector bcValueFactors;
  /// sum of bcValueFactors used for evaluating barycentric interpolant
  /// denominator term
  Real bcValueFactorSum;

  /// terms for gradients of barycentric interpolants
  // (bcValueFactors[j] * l(x) * sum of diff inverses)
  RealVector bcGradFactors;
};


inline LagrangeInterpPolynomial::LagrangeInterpPolynomial():
  InterpolationPolynomial(), newPoint(DBL_MAX), exactIndex(_NPOS),
  exactDeltaIndex(_NPOS)
{ }


inline LagrangeInterpPolynomial::
LagrangeInterpPolynomial(const RealArray& interp_pts):
  InterpolationPolynomial(interp_pts), newPoint(DBL_MAX), exactIndex(_NPOS),
  exactDeltaIndex(_NPOS)
{ }


inline LagrangeInterpPolynomial::~LagrangeInterpPolynomial()
{ }


/** Shared initialization code. */
inline void LagrangeInterpPolynomial::allocate_factors(short compute_order)
{
  size_t num_interp_pts = interpPts.size();
  if (bcWeights.length() != num_interp_pts) {
    PCerr << "Error: length of precomputed bcWeights (" << bcWeights.length()
	  << ") is inconsistent with number of collocation points ("
	  << num_interp_pts << ")." << std::endl;
    abort_handler(-1);
  }

  if ( (compute_order & 1) && bcValueFactors.length() != num_interp_pts)
    bcValueFactors.sizeUninitialized(num_interp_pts);
  if ( (compute_order & 2) && bcGradFactors.length()  != num_interp_pts)
    bcGradFactors.sizeUninitialized(num_interp_pts);
}


inline size_t LagrangeInterpPolynomial::exact_index() const
{ return exactIndex; }


inline size_t LagrangeInterpPolynomial::exact_delta_index() const
{ return exactDeltaIndex; }


//inline const RealVector& LagrangeInterpPolynomial::
//barycentric_weights() const
//{ return bcWeights; }


inline const RealVector& LagrangeInterpPolynomial::
barycentric_value_factors() const
{ return bcValueFactors; }


inline Real LagrangeInterpPolynomial::
barycentric_value_factor(unsigned short i) const
{ return bcValueFactors[i]; }


inline const RealVector& LagrangeInterpPolynomial::
barycentric_gradient_factors() const
{ return bcGradFactors; }


inline Real LagrangeInterpPolynomial::
barycentric_gradient_factor(unsigned short i) const
{ return bcGradFactors[i]; }


inline Real LagrangeInterpPolynomial::barycentric_value_factor_sum() const
{ return bcValueFactorSum; }


inline Real LagrangeInterpPolynomial::barycentric_difference_product() const
{ return diffProduct; }

} // namespace Pecos

#endif
