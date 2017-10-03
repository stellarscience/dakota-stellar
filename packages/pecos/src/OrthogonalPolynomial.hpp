/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        OrthogonalPolynomial
//- Description:  Abstract base class for orthogonal polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef ORTHOGONAL_POLYNOMIAL_HPP
#define ORTHOGONAL_POLYNOMIAL_HPP

#include "BasisPolynomial.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {


/// Base class for the orthogonal polynomial class hierarchy.

/** The OrthogonalPolynomial class is the base class for the
    univariate orthogonal polynomial class hierarchy in PECOS.  One
    instance of an OrthogonalPolynomial is created for each variable
    within a multidimensional orthogonal polynomial basis function (a
    vector of OrthogonalPolynomials is contained in
    OrthogPolyApproximation, which may be mixed and matched in, e.g.,
    the Wiener-Askey scheme for polynomial chaos). */

class OrthogonalPolynomial: public BasisPolynomial
{
public:

  //
  //- Heading: Constructors, destructor, assignment operator
  //

  OrthogonalPolynomial();  /// default constructor
  ~OrthogonalPolynomial(); /// destructor

  //
  //- Heading: Virtual function redefinitions
  //

  /// destroy history of Gauss pts/wts due to change in distribution parameters
  void reset_gauss();

  //
  //- Heading: Member functions
  //

  /// precompute tripleProductMap
  void precompute_triple_products(const UShortMultiSet& max_ijk);
  /// lookup value based on UShortMultiSet key within tripleProductMap;
  /// returns false if not stored
  bool triple_product(const UShortMultiSet& ijk_key, Real& trip_prod) const;
  /// lookup value based on three size_t keys within tripleProductMap;
  /// returns false if not stored
  bool triple_product(size_t i, size_t j, size_t k, Real& trip_prod) const;

  /// perform unit testing on Gauss points/weights
  void gauss_check(unsigned short order);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// set collocRule
  void collocation_rule(short rule);
  /// get collocRule
  short collocation_rule() const;

  //
  //- Heading: Data
  //

  /// collocation points for one-dimensional quadrature
  /// (x parameter values for which Poly_n(x) = 0)
  RealArray collocPoints;
  /// collocation weights for one-dimensional quadrature
  RealArray collocWeights;

  /// the type of integration rule associated with the orthogonal polynomial
  /** In most cases, this is just the corresponding Gauss quadrature
      rule.  However, for Legendre, collocRule manages the option of
      GAUSS_LEGENDRE or GAUSS_PATTERSON, for Chebyshev, it manages the
      option of CLENSHAW_CURTIS or FEJER2, and for Hermite, it manages
      the option of GAUSS_HERMITE or GENZ_KEISTER. */
  short collocRule;

private:

  //
  //- Heading: Data
  //

  /// mapping from an ijk sorted index set into <Psi_i Psi_j Psi_k>.
  /// These are precomputed with precompute_triple_products(order)
  /// and retrieved with triple_product(key)
  UShortMultiSetRealMap tripleProductMap;
  /// tracks precomputations to prevent redundancy
  UShortMultiSet tripleProductOrder;
};


inline OrthogonalPolynomial::OrthogonalPolynomial():
  BasisPolynomial(BaseConstructor())
{ }


inline OrthogonalPolynomial::~OrthogonalPolynomial()
{ }


inline void OrthogonalPolynomial::reset_gauss()
{ collocPoints.clear(); collocWeights.clear(); }


inline bool OrthogonalPolynomial::
triple_product(const UShortMultiSet& ijk_key, Real& trip_prod) const
{
  UShortMultiSetRealMap::const_iterator cit = tripleProductMap.find(ijk_key);
  if (cit == tripleProductMap.end())
    { trip_prod = 0.;          return false; }
  else
    { trip_prod = cit->second; return  true; }
}


inline bool OrthogonalPolynomial::
triple_product(size_t i, size_t j, size_t k, Real& trip_prod) const
{
  UShortMultiSet ijk_key;
  ijk_key.insert(i); ijk_key.insert(j); ijk_key.insert(k);
  return triple_product(ijk_key, trip_prod);
}


inline void OrthogonalPolynomial::collocation_rule(short rule)
{ collocRule = rule; }


inline short OrthogonalPolynomial::collocation_rule() const
{ return collocRule; }

} // namespace Pecos

#endif
