/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_ORTHOGONAL_POLY_BASIS_HPP
#define PECOS_SURROGATES_ORTHOGONAL_POLY_BASIS_HPP

//#include "AleatoryVariableTransformation.hpp"
#include "BasisPolynomial.hpp"
#include "teuchos_data_types.hpp"

namespace Pecos {
namespace surrogates {

/**
   \class OrthogonalPolynomialBasis()
   \brief A multivariate orthogonal polynomial basis

   This class was requires the PECOS library.
*/
class OrthogonalPolynomialBasis{
protected:
/**\brief Get the basis type from a variable type*/
static short get_basis_type(short var_type);

/**\brief Get the quadrature type from a basis type*/
static short get_quadrature_rule(short basis_type, bool nested_rule);

/**\brief Get the set of basis types from a set of variable types*/
static void get_basis_types(const Pecos::ShortArray& var_types, Pecos::ShortArray& basis_types);

public:
//std::shared_ptr<Surrogates::AleatoryVariableTransform> varTransform_;
std::vector<Pecos::BasisPolynomial> polynomialBasis_;
bool nestedRules_;

/// Constructor
OrthogonalPolynomialBasis();

/// Destructor
virtual ~OrthogonalPolynomialBasis();

/**\brief Initialize the polynomial basis from a set of basis types*/
void initialize_polynomial_basis_from_basis_types(const Pecos::ShortArray& basis_types);

/**\brief Initialize the polynomial basis from a set of variable types*/
void initialize_polynomial_basis_from_variable_types(const Pecos::ShortArray& var_types);
#ifdef DEPENDENCIES_NOT_IMPLEMENTED
/**\brief Initialize the polynomials to be orthogonal to the aleatory
   variables
*/
static void initialize_polynomial_basis(const std::shared_ptr<Surrogates::AleatoryVariables> vars);
#endif

/**\brief Evaluate a multivariate polynomial with a given index
*/
Real value(const RealVector& sample, const IntVector& indices);

/**\brief Evaluate a multivariate polynomial with a set of indices.
   TODO replace this function  with one that takes advantage of recursion
  //formulae of orthogonal polynomials
*/
void value(const RealMatrix& sample, const IntMatrix& indices,
           RealMatrix &values);

}; // class OrthogonalPolynomialBasis

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
