/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_POLYNOMIAL_CHAOS_EXPANSION_WRAP_HPP
#define PECOS_SURROGATES_POLYNOMIAL_CHAOS_EXPANSION_WRAP_HPP

#include "PolyApproximation.hpp"
#include "BasisApproximation.hpp"
#include "OptionsList.hpp"
#include "teuchos_data_types.hpp"
#include "pecos_data_types.hpp"
#include "SharedBasisApproxData.hpp"

namespace Pecos {
namespace surrogates {

/**
\class PolynomialChaosExpansionWrap

\brief Polynomial chaos expansion of a \f$L^2\f$ function, i.e a function with finite variance.
*/
class  PolynomialChaosExpansionWrap: public PolyApproximation{
protected:
/// TODO Do I need this???
Pecos::SharedBasisApproxData sharedData_;

/// The Pecos basis approximation which is being wrapped
Pecos::BasisApproximation poly_;

/// Model specific options
util::OptionsList opts_;

/// The number of QoI (outputs) of the vector-valued function
int numQOI_;

public:
  /// Default constructtor
  PolynomialChaosExpansionWrap();

  /// Destructor
  virtual ~PolynomialChaosExpansionWrap();


  /** \copydoc PolyApproximation::generate_basis_matrix() */
  void generate_basis_matrix(const RealMatrix &samples, RealMatrix &result_0);

}; // class PolynomialChaosExpansionWrap

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
