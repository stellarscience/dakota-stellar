/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_REGRESSION_BUILDER_HPP
#define PECOS_SURROGATES_REGRESSION_BUILDER_HPP

#include "SurrogateBuilder.hpp"
#include "math_tools.hpp"
#include "linear_solvers.hpp"


// \todo: to be removed once generalized sampler is created
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <VariableTransformation.hpp>

namespace Pecos {
namespace surrogates {

/** \brief Solve a regression problem for a given set of
    samples and values and set the coefficients of the
    approximation.
 */
void solve_regression(const RealMatrix &samples,
		      const RealMatrix &values,
		      util::OptionsList &opts,
		      Approximation &approx);

/** \brief helper function for gen_uniform_samples. \odo remove when gen_uniform_samples is removed*/
void get_canonical_uniform_samples(int M, int N, int seed, RealMatrix &samples);

  /** \brief Generate samples from a uniform distribution.
This should be replaced by a generalized sampling class
   */
void generate_uniform_samples(int num_vars, int num_samples, int seed, const VariableTransformation &var_transform, RealMatrix &samples );
  
  
/**
   \class RegressionBuilder
   \brief Class used to build surrogates using regression
*/
class RegressionBuilder : public SurrogateBuilder{
public:
  
  /// Constructor
  RegressionBuilder(){};
  
  /// Destructor
  virtual ~RegressionBuilder(){};
  
  /** \brief Build a surrogate to a set of specifications 
   * using regression.
   *
   * \param[in] opts Spefications used to build the surrogate
   *
   */
  virtual void build(util::OptionsList &opts, Approximation &approx,
                     util::OptionsList &result);
};

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
