/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_MONOMIAL_HPP
#define PECOS_SURROGATES_MONOMIAL_HPP

#include "PolyApproximation.hpp"

namespace Pecos {
namespace surrogates {

/**
\class Monomial
\brief A multivariate monomial approximation.

This class was generated mainly for unit-testing purposes. Much greater
functionality can be reached by including the PECOS library and utilizing
the polynomial chaos wrappers.
*/
class Monomial : public PolyApproximation {
public:
  Monomial();

  ~Monomial();

  void set_options(const util::OptionsList &opts);

  /** \copydoc PolyApproximation::generate_canonical_basis_matrix() */
  void generate_canonical_basis_matrix(const RealMatrix &samples, RealMatrix &result_0);

}; // class Monomial

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
