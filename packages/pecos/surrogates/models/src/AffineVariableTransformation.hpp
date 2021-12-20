/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_AFFINE_VARIABLE_TRANSFORMATION_HPP
#define PECOS_SURROGATES_AFFINE_VARIABLE_TRANSFORMATION_HPP

#include <BoundedVariables.hpp>
#include <teuchos_data_types.hpp>
#include <VariableTransformation.hpp>

namespace Pecos {
namespace surrogates {

  /**
     \class AffineVariableTransformation
     \brief Defines transformations of BOUNDED random variables from
     native function variable space to the canonical variable space of
     an approximation.

     Transforms $d$-dimensional samples, derivatives, and Hessians
     from user defined x-space \f$[a_i,b_i,\ldots,a_d,b_d]\f$
     to standardized u-space \f$[-1,1]^d\f$
  */

  class AffineVariableTransformation : public VariableTransformation{
  protected:
    std::shared_ptr<BoundedVariables> boundedVars_;

  public:
    AffineVariableTransformation();

    virtual ~AffineVariableTransformation();

    /** \copydoc VariableTransformation::set_variables() */
    void set_variables(const std::shared_ptr<Variables> &vars);

    /** \copydoc VariableTransformation::map_samples_to_user_space() */
    virtual void map_samples_to_user_space(const RealMatrix &samples,
					   RealMatrix &result) const;

    /** \copydoc VariableTransformation::map_samples_from_user_space() */
    virtual void map_samples_from_user_space(const RealMatrix &samples, RealMatrix &result) const;

    /** \copydoc VariableTransformation::map_derivatives_to_user_space() */
    virtual void map_derivatives_to_user_space(const RealVector &derivatives,
					       int dim, RealVector &transformed_derivatives) const;

    /** \copydoc VariableTransformation::map_derivatives_from_user_space() */
    virtual void map_derivatives_from_user_space(const RealVector &derivatives,
						 int dim, RealVector &transformed_derivatives) const;

    /** \copydoc VariableTransformation::map_hessian_from_user_space() */
    virtual void map_hessian_from_user_space(const RealMatrix &hessian,
					     RealMatrix &transformed_hessian) const{throw( std::runtime_error("Not yet implemented") );};

    /** \copydoc VariableTransformation::map_hessian_to_user_space() */
    virtual void map_hessian_to_user_space(const RealMatrix &hessian,
					   RealMatrix &transformed_hessian) const{throw( std::runtime_error("Not yet implemented") );};

    /** \copydoc VariableTransformation::set_options() */
    virtual void set_options(util::OptionsList &opts){};

  };

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
