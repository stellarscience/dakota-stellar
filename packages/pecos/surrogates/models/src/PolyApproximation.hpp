/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_POLY_APPROXIMATION_HPP
#define PECOS_SURROGATES_POLY_APPROXIMATION_HPP

#include "Approximation.hpp"
#include "OptionsList.hpp"
#include "teuchos_data_types.hpp"

namespace Pecos {
namespace surrogates {

/**
\class PolyApproximation

\brief A polynomial approximation \f$p(x)\f$ of a function \f$f(x)\f$.

\f$f(x):\mathbb{R}^d\rightarrow\mathbb{R}^q\f$ is a vector-valued function
parameterized by a vector of variables \f$x=(x_1,\ldots,x_d)\in\mathbb{R}^d\f$.
The polynomial approximation of \f$p(x)\f$ takes the form
\f[f(x) \approx p(x) = \sum_{\boldsymbol{\lambda}\in \boldsymbol{\Lambda}}
 \alpha_{\boldsymbol{\lambda}}\phi_{\boldsymbol{\lambda}}(x) \f]
for some set of \f$\boldsymbol{\Lambda}\f$ of multivariate indices
\f$\boldsymbol{\lambda}=(\lambda_1,\ldots,\lambda_d)\f$

This class provides an interface to the polynomial chaos expansion software
implemented in Pecos
 */
class PolyApproximation : public Approximation {
protected:

  IntMatrix basisIndices_;

  RealMatrix basisCoeffs_;

public:

  /// Default constructor
  PolyApproximation();

  /// Destructor
  virtual ~PolyApproximation();

  /** \copydoc Function::value() */
  virtual void value(const RealMatrix &samples, RealMatrix &result_0);

  /** \copydoc Function::gradient() */
  virtual void gradient(const RealMatrix &samples, int qoi, RealMatrix &result_0);

  /** \copydoc Function::jacobian() */
  virtual void jacobian(const RealVector &sample, RealMatrix &result_0);

  /** \copydoc Function::hessian() */
  virtual void hessian(const RealMatrix &samples, int qoi, RealMatrixList &hessians);

  /**\brief Evaluate each basis function at a set of samples in the 
     canonical space of the polynomial basis

   Build the matrix \f$V\f$, where \f$V_{ij}=\phi_j(x^{(i)})\f$.

   \param[in] samples (num_vars x num_samples) matrix
       The coordindates of the samples x in the canonical space of 
       the polynomial basis.

   \param[out] result_0
       The matrix \f$V\f$, where \f$V_{ij}=\phi_j(x^{(i)})\f$.
   */
  virtual void generate_canonical_basis_matrix(const RealMatrix &samples, RealMatrix &result_0)=0;


/**\brief Evaluate each basis function at a set of samples in the 
     user-defined space.

   Build the matrix \f$V\f$, where \f$V_{ij}=\phi_j(x^{(i)})\f$.

   \param[in] samples (num_vars x num_samples) matrix
       The coordindates of the samples x in the user-defined space.

   \param[out] result_0
       The matrix \f$V\f$, where \f$V_{ij}=\phi_j(x^{(i)})\f$.
   */
  
  void generate_basis_matrix(const RealMatrix &samples, RealMatrix &result_0);

  /**
  \brief Set the coefficients of the polynomial expansion.

  \param[in] coeffs (num_terms x num_qoi) matrix
      The coefficients of the basis for each QoI
  */
  void set_coefficients(const RealMatrix &coeffs);

    /**
  \brief Get the coefficients of the polynomial expansion.

  \param[out] result (num_terms x num_qoi) matrix
      The coefficients of the basis for each QoI
  */
  void get_coefficients(RealMatrix &result) const;

  /*
   \brief Set the multivariate index of the polynomial basis

   \param[in] basis_indices (num_dims x num_indices) matrix
       The multivariate indices
  */
  void set_basis_indices(const IntMatrix &basis_indices);

  /* \brief Get the multivariate index of the polynomial basis

   \param[out] result (num_dims x num_indices) matrix
       The multivariate indices
  */
  void get_basis_indices(IntMatrix &result) const;

  /* \brief Get the number of terms in the polynomial

   \return num_terms
  */
  int num_terms() const;

}; // class PolyApproximation

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
