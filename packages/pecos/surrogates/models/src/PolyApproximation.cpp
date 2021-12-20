/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "PolyApproximation.hpp"
#include "teuchos_data_types.hpp"

namespace Pecos {
namespace surrogates {

PolyApproximation::PolyApproximation(){};

PolyApproximation::~PolyApproximation(){};

void PolyApproximation::value(const RealMatrix &samples, RealMatrix &approx_vals) {
  if (basisCoeffs_.numRows()!=basisIndices_.numCols())
    throw(std::runtime_error("basis indices and coefficients are inconsistent"));
  util::resize_if_needed(approx_vals,samples.numCols(),basisCoeffs_.numCols());
  RealMatrix basis_matrix;
  generate_basis_matrix(samples, basis_matrix);
  util::multiply(basis_matrix,basisCoeffs_,approx_vals,1.0,0.0);
}

void PolyApproximation::
generate_basis_matrix(const RealMatrix &samples, RealMatrix &basis_matrix){
  // BMA TODO: Wasn't sure if the intent is to have a default no-op
  // variable transformation, or conditional logic here to avoid the
  // copy...
  if (varTransform_) {
    RealMatrix transformed_samples;
    varTransform_->map_samples_from_user_space(samples,transformed_samples);
    generate_canonical_basis_matrix(transformed_samples, basis_matrix);
  }
  else
    generate_canonical_basis_matrix(samples, basis_matrix);
}


void PolyApproximation::set_coefficients(const RealMatrix &coeffs) {
    basisCoeffs_.shapeUninitialized(coeffs.numRows(), coeffs.numCols());
    basisCoeffs_.assign(coeffs);  // TODO: compare dimensions to basisIndices?
}

void PolyApproximation::get_coefficients(RealMatrix &coeffs) const{
  coeffs = basisCoeffs_;
}

void PolyApproximation::set_basis_indices(const IntMatrix &basis_indices){
  if (num_vars()!=basis_indices.numRows())
    throw(std::runtime_error("basis indices must be num_vars x num_indices"));
  util::resize_if_needed(basisIndices_, basis_indices.numRows(),
		   basis_indices.numCols());
  basisIndices_.assign(basis_indices);
};

void PolyApproximation:: get_basis_indices(IntMatrix &basis_indices) const{
  basis_indices = basisIndices_;
}

void PolyApproximation::gradient(const RealMatrix &samples, int qoi, RealMatrix &gradients) {
  throw(std::string("This Function type does not provide gradients"));
}

void PolyApproximation::jacobian(const RealVector &sample, RealMatrix &jacobian) {
  throw(std::string("This Function type does not provide a jacobian"));
}

void PolyApproximation::hessian(const RealMatrix &samples, int qoi, RealMatrixList &hessians) {
  throw(std::string("This Function type does not provide gradients"));
}

int PolyApproximation::num_terms() const{
  return basisIndices_.numCols();
}

}  // namespace surrogates
}  // namespace Pecos
