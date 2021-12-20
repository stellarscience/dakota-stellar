/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "OrthogonalPolynomialBasis.hpp"
#include "pecos_global_defs.hpp"
#include "teuchos_data_types.hpp"
#include "math_tools.hpp"

namespace Pecos {
namespace surrogates {

OrthogonalPolynomialBasis::OrthogonalPolynomialBasis() :nestedRules_(true){};
OrthogonalPolynomialBasis::~OrthogonalPolynomialBasis(){};

short OrthogonalPolynomialBasis::get_basis_type(short var_type){
  short basis_type=Pecos::NUM_GEN_ORTHOG;
  switch (var_type){
  case Pecos::STD_NORMAL:
    basis_type  = Pecos::HERMITE_ORTHOG; break;
  case Pecos::STD_UNIFORM:
    basis_type = Pecos::LEGENDRE_ORTHOG; break;
  case Pecos::STD_EXPONENTIAL:
    basis_type = Pecos::LAGUERRE_ORTHOG; break;
  case Pecos::STD_BETA:
    basis_type = Pecos::JACOBI_ORTHOG; break;
  case Pecos::STD_GAMMA:
    basis_type = Pecos::GEN_LAGUERRE_ORTHOG; break;
  case Pecos::POISSON:
    basis_type = Pecos::CHARLIER_DISCRETE; break;
  case Pecos::BINOMIAL:
    basis_type = Pecos::KRAWTCHOUK_DISCRETE; break;
  case Pecos::NEGATIVE_BINOMIAL:
    basis_type = Pecos::MEIXNER_DISCRETE; break;
  default:
    basis_type = Pecos::NUM_GEN_ORTHOG; break;
  }
  return basis_type;
}

short OrthogonalPolynomialBasis::
get_quadrature_rule(short basis_type, bool nested_rule){
  short quad_rule=Pecos::GOLUB_WELSCH;
  switch (basis_type) {
  case Pecos::HERMITE_ORTHOG:
    quad_rule = (nested_rule) ? Pecos::GENZ_KEISTER : Pecos::GAUSS_HERMITE;
    break;
  case Pecos::LEGENDRE_ORTHOG:
    quad_rule = (nested_rule) ? Pecos::GAUSS_PATTERSON : Pecos::GAUSS_LEGENDRE;
    break;
  case Pecos::LAGUERRE_ORTHOG:
    quad_rule = Pecos::GAUSS_LAGUERRE; break;
  case Pecos::JACOBI_ORTHOG:
    quad_rule = Pecos::GAUSS_JACOBI; break;
  case Pecos::GEN_LAGUERRE_ORTHOG:
    quad_rule = Pecos::GEN_GAUSS_LAGUERRE; break;
  case Pecos::CHARLIER_DISCRETE:
  quad_rule = Pecos::GAUSS_CHARLIER; break;
  case Pecos::KRAWTCHOUK_DISCRETE:
    quad_rule = Pecos::GAUSS_KRAWTCHOUK; break;
  case Pecos::MEIXNER_DISCRETE:
    quad_rule = Pecos::GAUSS_MEIXNER; break;
  default:
    quad_rule = Pecos::GOLUB_WELSCH; break;
  }
  return quad_rule;
}

void OrthogonalPolynomialBasis::
get_basis_types(const Pecos::ShortArray& var_types,
                Pecos::ShortArray& basis_types){
  size_t num_vars = var_types.size();
  basis_types.resize(num_vars);
  for (size_t i=0; i<num_vars; ++i){
    basis_types[i] = get_basis_type(var_types[i]);
  }
}

void OrthogonalPolynomialBasis::
initialize_polynomial_basis_from_variable_types(const Pecos::ShortArray& var_types){
  Pecos::ShortArray basis_types;
  get_basis_types(var_types, basis_types);
  initialize_polynomial_basis_from_basis_types(basis_types);
}

void OrthogonalPolynomialBasis::
initialize_polynomial_basis_from_basis_types(const Pecos::ShortArray& basis_types){
  size_t num_vars = basis_types.size();
  polynomialBasis_.resize(num_vars);
  for (size_t i=0; i<num_vars; ++i){
    short quad_rule = get_quadrature_rule(basis_types[i], nestedRules_);
    polynomialBasis_[i] = Pecos::BasisPolynomial(basis_types[i], quad_rule);
  }
}

Real OrthogonalPolynomialBasis::
value(const RealVector& sample, const IntVector& indices){
  if (polynomialBasis_.size()==0)
    throw(std::runtime_error("polynomialBasis_ has not been initialized"));
  Real val = 1.; unsigned short order_1d; int num_vars = sample.length();
  for (size_t i=0; i<num_vars; ++i) {
    order_1d = indices[i];
    if (order_1d)
      val *= polynomialBasis_[i].type1_value(sample[i], order_1d);
  }
  return val;
}

void OrthogonalPolynomialBasis::
value(const RealMatrix& samples, const IntMatrix& indices, RealMatrix &values){
  if (polynomialBasis_.size()==0)
    throw(std::runtime_error("polynomialBasis_ has not been initialized"));
  unsigned short order_1d; int num_vars = samples.numRows();
  int i,j,d, num_exp_terms = indices.numCols(), num_samples=samples.numCols();
  util::resize_if_needed(values, num_samples, num_exp_terms);
  for (j=0; j<num_exp_terms; ++j){
    for (i=0; i<num_samples; ++i){
      values(i,j) = 1.;
      for (d=0; d<num_vars; ++d) {
	order_1d = indices(d,j);
	if (order_1d)
	  values(i,j) *= polynomialBasis_[d].type1_value(samples(d,i),order_1d);
        }
    }
  }
}

#ifdef DEPENDENCIES_NOT_IMPLEMENTED
void OrthogonalPolynomialBasis::
initialize_polynomial_basis(const std::shared_ptr<Surrogates::AleatoryVariables> vars){

  std::shared_ptr<Surrogates::AleatoryVariables> aleatory_vars = Teuchos::rcp_dynamic_cast<Surrogates::AleatoryVariables>(vars,true);
  if (aleatory_vars.is_null())
    throw(std::runtime_error("vars is not an object of type AleatoryVariables"));
  size_t = aleatory_vars.num_vars();
  Pecos::ShortArray var_types(num_vars);
  for (size_t i=0; i<num_vars; ++i)
    var_types[i]=aleatory_vars.type(i);
  initialize_polynomial_basis_from_variable_types(var_types);
}
#endif

}  // namespace surrogates
}  // namespace Pecos
