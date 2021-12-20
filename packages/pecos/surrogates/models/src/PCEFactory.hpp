/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_PCE_FACTORY_HPP
#define PECOS_SURROGATES_PCE_FACTORY_HPP

//#include "Function.hpp"
#include "Monomial.hpp"
#include "PolynomialChaosExpansion.hpp"

namespace Pecos {
namespace surrogates {

  enum {MONOMIAL, PCE};

std::shared_ptr<PolyApproximation> polynomial_approximation_factory(
	  const std::shared_ptr<VariableTransformation> &var_transform,
	  const util::OptionsList& opts) {
  int poly_type = opts.get<int>("poly_type");
  switch (poly_type){
  case MONOMIAL : {
    std::shared_ptr<Monomial> poly(new Monomial);
    poly->set_variable_transformation(var_transform);
    return poly;
  }
  case PCE : {
    std::shared_ptr<PolynomialChaosExpansion> poly(new PolynomialChaosExpansion);
    IntVector basis_types;
    basis_types = opts.get("basis_types", basis_types);
    if (basis_types.length()>0){
      // OptionsList can only pass in Teuchos Vectors so convert
      // to std::vector here
      Pecos::ShortArray short_basis_types(basis_types.length());
      for (int i=0; i<basis_types.length(); ++i)
        short_basis_types[i]=static_cast<short>(basis_types[i]);
      poly->initialize_polynomial_basis_from_basis_types(short_basis_types);
    }
    // else //basis will be chosen based upon info in variable transformation
    //   {throw(std::runtime_error("setting pce basis from variables not implemented"));}
    poly->set_variable_transformation(var_transform);
    return poly;
  }
  default : {
    throw(std::runtime_error("Incorrect poly_type"));
  }
  }
}

  
}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
