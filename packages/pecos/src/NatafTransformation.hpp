/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef NATAF_TRANSFORMATION_HPP
#define NATAF_TRANSFORMATION_HPP

#include "ProbabilityTransformation.hpp"


namespace Pecos {


/// Class for Nataf nonlinear distribution transformation.

/** The Nataf transformation occurs in two steps: (1) transformation
    from the original correlated distributions (x-space) to correlated
    standard normals (z-space) using CDF matching and from correlated
    standard normals to uncorrelated standard normals (u-space) using
    the inverse Cholesky factor of a modified correlation matrix. */

class NatafTransformation: public ProbabilityTransformation
{
public:

  //
  //- Heading: Constructors and destructor
  //

  NatafTransformation();  ///< constructor
  ~NatafTransformation(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// Transformation routine from u-space of uncorrelated standard normal
  /// variables to x-space of correlated random variables
  void trans_U_to_X(const RealVector& u_vars, RealVector& x_vars);

  /// Transformation routine from x-space of correlated random variables 
  /// to u-space of uncorrelated standard normal variables
  void trans_X_to_U(const RealVector& x_vars, RealVector& u_vars);

  /// As part of the Nataf distribution model (Der Kiureghian & Liu, 1986),
  /// this procedure modifies the user-specified correlation matrix
  /// (corrMatrixX) to account for correlation warping from the nonlinear
  /// X->Z transformation and performs a Cholesky factorization to create
  /// corrCholeskyFactorZ.
  void transform_correlations();

  /// Transformation routine for gradient vector from x-space to u-space
  void trans_grad_X_to_U(const RealVector& fn_grad_x, RealVector& fn_grad_u,
			 const RealVector& x_vars,
			 const SizetArray& x_dvv,
			 SizetMultiArrayConstView cv_ids);
  /// Transformation routine for gradient vector from x-space to u-space
  void trans_grad_X_to_U(const RealVector& fn_grad_x,   RealVector& fn_grad_u,
			 const RealMatrix& jacobian_xu,
			 const SizetArray& x_dvv,
			 SizetMultiArrayConstView cv_ids);

  /// Transformation routine from x-space gradient vector to design space
  void trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
			 const RealVector& x_vars, const SizetArray& x_dvv,
			 SizetMultiArrayConstView cv_ids,
			 SizetMultiArrayConstView acv_ids,
			 const SizetArray& acv_map1_indices,
			 const ShortArray& acv_map2_targets);
  /// Transformation routine from x-space gradient vector to design space
  void trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
			 const RealMatrix& jacobian_xs, const SizetArray& x_dvv,
			 SizetMultiArrayConstView cv_ids,
			 SizetMultiArrayConstView acv_ids,
			 const SizetArray& acv_map1_indices,
			 const ShortArray& acv_map2_targets);

  /// Transformation routine for gradient vector from u-space to x-space
  void trans_grad_U_to_X(const RealVector& fn_grad_u, RealVector& fn_grad_x,
			 const RealVector& x_vars, const SizetArray& x_dvv,
			 SizetMultiArrayConstView cv_ids);
  /// Transformation routine for gradient vector from u-space to x-space
  void trans_grad_U_to_X(const RealVector& fn_grad_u,   RealVector& fn_grad_x,
			 const RealMatrix& jacobian_ux, const SizetArray& x_dvv,
			 SizetMultiArrayConstView cv_ids);

  /// Transformation routine for Hessian matrix from x-space to u-space
  void trans_hess_X_to_U(const RealSymMatrix& fn_hess_x,
			 RealSymMatrix& fn_hess_u, const RealVector& x_vars,
			 const RealVector& fn_grad_x, const SizetArray& x_dvv,
			 SizetMultiArrayConstView cv_ids);
  /// Transformation routine for Hessian matrix from x-space to u-space
  void trans_hess_X_to_U(const RealSymMatrix& fn_hess_x,
			 RealSymMatrix& fn_hess_u,
			 const RealMatrix& jacobian_xu,
			 const RealSymMatrixArray& hessian_xu,
			 const RealVector& fn_grad_x, const SizetArray& x_dvv,
			 SizetMultiArrayConstView cv_ids);

  /// Jacobian of x(u) mapping obtained from dX/dZ dZ/dU
  void jacobian_dX_dU(const RealVector& x_vars, RealMatrix& jacobian_xu);

  /// Jacobian of u(x) mapping obtained from dU/dZ dZ/dX
  void jacobian_dU_dX(const RealVector& x_vars, RealMatrix& jacobian_ux);

  /// Design Jacobian of x(u,s) mapping obtained from differentiation of
  /// trans_U_to_X() with respect to distribution parameters S
  void jacobian_dX_dS(const RealVector& x_vars, RealMatrix& jacobian_xs,
		      SizetMultiArrayConstView cv_ids,
		      SizetMultiArrayConstView acv_ids,
		      const SizetArray& acv_map1_indices,
		      const ShortArray& acv_map2_targets);

  /// Hessian of x(u) mapping obtained from dZ/dU^T d^2X/dZ^2 dZ/dU
  void hessian_d2X_dU2(const RealVector& x_vars,
		       RealSymMatrixArray& hessian_xu);

private:

  //
  //- Heading: Utility routines
  //

  /// Transformation routine from u-space of uncorrelated standard normal
  /// variables to z-space of correlated standard normal variables
  void trans_U_to_Z(const RealVector& u_vars, RealVector& z_vars);

  /// Transformation routine from z-space of correlated standard normal
  /// variables to x-space of correlated random variables
  void trans_Z_to_X(const RealVector& z_vars, RealVector& x_vars);
  /// Transformation routine from a single z-space variable to a
  /// corresponding x-space variable
  void trans_Z_to_X(Real z, Real& x, size_t i);

  /// Transformation routine from x-space of correlated random variables
  /// to z-space of correlated standard normal variables
  void trans_X_to_Z(const RealVector& x_vars, RealVector& z_vars);
  /// Transformation routine from a single x-space random variable to
  /// a corresponding z-space variable
  void trans_X_to_Z(Real x, Real& z, size_t i);

  /// Transformation routine from z-space of correlated standard normal
  /// variables to u-space of uncorrelated standard normal variables
  void trans_Z_to_U(RealVector& z_vars, RealVector& u_vars);

  /// Jacobian of x(z) mapping obtained from differentiation of trans_Z_to_X()
  void jacobian_dX_dZ(const RealVector& x_vars, RealMatrix& jacobian_xz);

  /// Jacobian of z(x) mapping obtained from differentiation of trans_X_to_Z()
  void jacobian_dZ_dX(const RealVector& x_vars, RealMatrix& jacobian_zx);

  /// Hessian of x(z) mapping obtained from differentiation of jacobian_dX_dZ()
  void hessian_d2X_dZ2(const RealVector& x_vars,
		       RealSymMatrixArray& hessian_xz);
};


inline NatafTransformation::NatafTransformation():
  ProbabilityTransformation(BaseConstructor())
{ }


inline NatafTransformation::~NatafTransformation()
{ }

} // namespace Pecos

#endif
