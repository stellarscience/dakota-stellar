/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        InterpolationPolynomial
//- Description:  Class for 1-D Interpolation Polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef INTERPOLATION_POLYNOMIAL_HPP
#define INTERPOLATION_POLYNOMIAL_HPP

#include "BasisPolynomial.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {


/// Derived basis polynomial class for 1-D Lagrange interpolation polynomials

/** The InterpolationPolynomial class evaluates a univariate
    interpolation polynomial.  It enables multidimensional
    interpolants within InterpPolyApproximation. */

class InterpolationPolynomial: public BasisPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  InterpolationPolynomial();
  /// standard constructor
  InterpolationPolynomial(const RealArray& interp_pts);
  /// destructor
  virtual ~InterpolationPolynomial();

  //
  //- Heading: New virtual functions
  //

  /// precompute data that is reused repeatedly within interpolation
  virtual void precompute_data();

  //
  //- Heading: Virtual function redefinitions
  //

  /// get size of interpPts
  size_t interpolation_size() const;

  /// set interpPts
  void interpolation_points(const RealArray& interp_pts);
  /// get interpPts
  const RealArray& interpolation_points() const;

protected:

  //
  //- Heading: Data
  //

  /// set of 1-D interpolation points: the i_th interpolation polynomial
  /// evaluated at the j_th interpolation point produces Kronecker delta_ij
  RealArray interpPts;

private:

  //
  //- Heading: Data
  //

};


inline InterpolationPolynomial::InterpolationPolynomial():
  BasisPolynomial(BaseConstructor())
{ }


inline InterpolationPolynomial::
InterpolationPolynomial(const RealArray& interp_pts):
  BasisPolynomial(BaseConstructor()), interpPts(interp_pts)
{ precompute_data(); }


inline InterpolationPolynomial::~InterpolationPolynomial()
{ }


inline size_t InterpolationPolynomial::interpolation_size() const
{ return interpPts.size(); }


inline void InterpolationPolynomial::
interpolation_points(const RealArray& interp_pts)
{
  interpPts = interp_pts;
  precompute_data();
}


inline const RealArray& InterpolationPolynomial::interpolation_points() const
{ return interpPts; }

} // namespace Pecos

#endif
