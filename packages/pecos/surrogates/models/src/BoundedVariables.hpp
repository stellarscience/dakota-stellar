/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_BOUNDED_VARIABLES_HPP
#define PECOS_SURROGATES_BOUNDED_VARIABLES_HPP

#include <Variables.hpp>
#include <teuchos_data_types.hpp>

namespace Pecos {
namespace surrogates {

  /**
     \class BoundedVariables
     \brief The object representing the meta data of a vector of bounded
     variables.
  */
  class BoundedVariables : public Variables{
  public:
    BoundedVariables();

    virtual ~BoundedVariables();

    /**\brief return the upper bound of the ith variable

     \param[in] i The ith dimension for which the upper bound is requested
     */
    Real ub(int i) const;

    /**\brief return the lower bound of the ith variable

     \param[in] i The ith dimension for which the lower bound is requested
     */
    Real lb(int i) const;

    /** \brief set the ranges of the variables.

     \param[in] ranges (2*num_vars x 1) vector
          [lb1,ub1,...,lbd,upb] where lb1 and ub1 are the lower and
	  upper bound of the 1st variable and so on.
     */
    void set_ranges(const RealVector &ranges);

    void set_options(const util::OptionsList &opts);

  private:

    /// The ranges of the bounded random variables
    RealVector ranges_;

  }; // class BoundedVariables


void define_homogeneous_ranges(int num_vars, Real lb, Real ub,
			       RealVector &result); 

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
