/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_VARIABLES_HPP
#define PECOS_SURROGATES_VARIABLES_HPP

#include <OptionsList.hpp>

namespace Pecos {
namespace surrogates {

  /**
     \class Variables
     \brief The object representing the meta data of a vector of variables.
  */
  class Variables {
  public:
    Variables();

    virtual ~Variables();

    /**\brief Set options specific to the model
     */
    virtual void set_options(const util::OptionsList &opts);

    /**\brief Get the model specific options
     */
    virtual void get_options(util::OptionsList &opts);

    /**\brief Set the number of variables
     */
    
    void set_num_vars(int num_vars);

    /**\brief Get the number of variables
     */
    int num_vars();

  protected:
    /// The number of variables
    int numVars_;

  }; // class Variables


}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
