/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_VARIABLE_TRANSFORMATION_HPP
#define PECOS_SURROGATES_VARIABLE_TRANSFORMATION_HPP

#include <Variables.hpp>
#include <memory>
#include <teuchos_data_types.hpp>

namespace Pecos {
namespace surrogates {

  /**
     \class VariableTransformation
     \brief Defines transformations of random variables from native function
     variable space to the canonical variable space of an approximation.

     This class should be used for any transformation including Nataf, Rosenblatt
     and rotations (such as those coming from active subspace methods)

     We should allow for combinations of correlated variables, data informed
     densities, e.g. KDE, and tensor product densities. E.g variable is a list
     comprising of at least one of these three types, but possibly all and multiple
     different objects of same class.

     In this scenario it is likely this class will just be a container
     that calls tranformations on each of the supported variable types.
  */

  class VariableTransformation {
  protected:
    /// Variable meta data describing user-defined x-space
    /// where evaluations are performed
    std::shared_ptr<Variables> vars_;

  public:
    VariableTransformation();

    virtual ~VariableTransformation();

    /**
       \brief Tranform a set of samples from transformed space (u-space) to the
       original user-defined x-space (where evaluations are performed)
       via application of the Jacobian dx/du.

       \param[in] samples (num_dims x num_samples) matrix
       Coordinates of the samples in the u-space

       \param[out] transformed_samples
       The coordinates of the transformed samples in the user x-space
    */
    virtual void map_samples_to_user_space(const RealMatrix &samples,
					   RealMatrix &transformed_samples) const=0;

    /**
       \brief Tranform a set of samples from the original
       user-defined x-space (where evaluations are performed)
       to another variable space (u-space) via application of the
       Jacobian dx/du.

       \param[in] samples (num_dims x num_samples) matrix
       Coordinates of the samples in the x-space

       \param[out] transformed_samples
       The coordinates of the transformed samples in the canonical u-space
    */
    virtual void map_samples_from_user_space(const RealMatrix &samples, RealMatrix &transformed_samples) const=0;

    /**
       \brief Tranform a derivative vector dg/dx from the original
       user-defined x-space (where evaluations are performed) to another variables
       space (u-space) via application of the Jacobian dx/du.

       \param[in] derivatives (num_derivatives x 1) vector
       Derivatives in the ith (i=dim) dimension, in the x-space

       \param[in] dim
       The direction of the derivatives

       \param[out] transformed_derivatives
       The transformed derivatives in the canonical u-space
    */
    virtual void map_derivatives_to_user_space(const RealVector &derivatives,
					       int dim, RealVector &transformed_derivatives) const{throw( std::runtime_error("Not yet implemented") );};

    /**
       \brief Tranform a derivative vector dg/du from a transformed variable
       space (u-space) to the original user-defined x-space (where evaluations are
       performed)via application of the Jacobian dx/du.

       \param[in] derivatives (num_derivatives x 1) vector
       Derivatives in the ith (i=dim) dimension, in the u-space

       \param[in] dim
       The direction of the derivatives

       \param[out] transformed_derivatives
       The transformed derivatives in the user-defined x-space
    */
    virtual void map_derivatives_from_user_space(const RealVector &derivatives,
						 int dim, RealVector &transformed_derivatives) const{throw( std::runtime_error("Not yet implemented") );};

    /**
       \brief Tranform a Hessian matrix from the original
       user-defined x-space (where evaluations are performed) to
       another variable space (u-space).

       \param[in] hessian (num_vars x num_vars) matrix
       Hessian matrix computed in the x-space

       \param[out] transformed_hessian
       The transformed hessian in the u-space coordinates
    */
    virtual void map_hessian_from_user_space(const RealMatrix &hessian,
					     RealMatrix &transformed_hessian) const{throw( std::runtime_error("Not yet implemented") );};

    /**
       \brief Tranform a Hessian matrix from a transformed variable space
       (u-space) to the original user-defined x-space (where evaluations are
       performed).

       \param[in] hessian (num_derivatives) matrix
       Hessian computed in the u-space

       \param[out] transformed_hessian
       The transformed hessian in the user-defined x-space
    */
    virtual void map_hessian_to_user_space(const RealMatrix &hessian,
					   RealMatrix &transformed_hessian) const{throw( std::runtime_error("Not yet implemented") );};

    /** \brief Set options of the transformation

	For example (in pecos language) set transformation is from x to u space
	of x to z space or u to z etc.
    */
    virtual void set_options(util::OptionsList &opts){}

    /**\brief Set the variables defining the user x-space
     */
    virtual void set_variables(const std::shared_ptr<Variables> &vars);

    /**\copydoc Variables::num_vars() */
    int num_vars();

  };

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
