/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "Monomial.hpp"

namespace Pecos {
namespace surrogates {

Monomial::Monomial(){}

Monomial::~Monomial(){}

void Monomial::set_options(const util::OptionsList &opts){
  PolyApproximation::set_options(opts);
}

void Monomial::
generate_canonical_basis_matrix(const RealMatrix &samples, RealMatrix &basis_matrix){
  int num_vars = varTransform_->num_vars(),
    num_terms=basisIndices_.numCols(), num_samples=samples.numCols();
  Teuchos::BLAS<int, Real> blas;
  int amax = blas.IAMAX( num_samples*num_vars, samples.values(), 1 ) - 1;
  if (samples.values()[amax]>1.){
    std::string msg =
      "Monomial::genereate_canonical_basis() Samples must be in [-1,1]";
    throw(std::runtime_error(msg));
  }
  util::resize_if_needed(basis_matrix, num_samples, num_terms);
  for(int j=0; j<num_terms; ++j){
    for(int i=0; i<num_samples; ++i){
      Real val = 1.;
      for (int d=0; d<num_vars; d++)
	val *= std::pow(samples(d,i), basisIndices_(d,j));
      basis_matrix(i,j) = val;
    }
  }
}

}  // namespace surrogates
}  // namespace Pecos
