/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_POLYNOMIAL_CHAOS_EXPANSION_HPP
#define PECOS_SURROGATES_POLYNOMIAL_CHAOS_EXPANSION_HPP

#include "PolyApproximation.hpp"
#include "OrthogonalPolynomialBasis.hpp"
#include "AleatoryVariableTransformation.hpp"

namespace Pecos {
namespace surrogates {

/**
\class PolynomialChaosExpansion
\brief A multivariate polynomial chaos expansion approximation.

This class was requires the PECOS library.
*/
class PolynomialChaosExpansion : public PolyApproximation {
private:
  OrthogonalPolynomialBasis orthogPolyBasis_;
public:

  //std::shared_ptr<Surrogates::AleatoryVariableTransformation> aleatoryVarTransform_;

  PolynomialChaosExpansion();

  ~PolynomialChaosExpansion();

  void set_options(const util::OptionsList &opts);

  /** \copydoc PolyApproximation::generate_canonical_basis_matrix() */
  void generate_canonical_basis_matrix(const RealMatrix &samples, RealMatrix &result_0);

  /** \copydoc PolyApproximation::set_variable_transformation() */
  void set_variable_transformation(const std::shared_ptr<VariableTransformation> var_transform);

  void initialize_polynomial_basis_from_basis_types(const Pecos::ShortArray &basis_types);

}; // class PolynomialChaosExpansion

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
