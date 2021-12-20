/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_SURROGATE_BUILDER_HPP
#define PECOS_SURROGATES_SURROGATE_BUILDER_HPP

#include "Approximation.hpp"
#include "OptionsList.hpp"

namespace Pecos {
namespace surrogates {

/**
   \class SurrogateBuilder
   \brief Abstract base class for methods used to construct 
   * surrogates
*/
class SurrogateBuilder {
public:

  /// Function that is being emulated
  Function* targetFunction_;
  
  /// Constructor
  SurrogateBuilder(){};
  
  /// Destructor
  virtual ~SurrogateBuilder(){};

  /** \brief Build a surrogate to a set of specifications.
   *
   * \param[in] opts Spefications used to build the surrogate
   *
   */
  virtual void build(util::OptionsList &opts, Approximation &approx,
                     util::OptionsList &result) = 0;

  void set_target_function(Function &target_function){
    targetFunction_ = &target_function;};
};

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
