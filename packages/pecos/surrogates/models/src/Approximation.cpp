/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "Approximation.hpp"

namespace Pecos {
namespace surrogates {

void Approximation::
set_variable_transformation(const std::shared_ptr<VariableTransformation> &var_transform){
  varTransform_ = var_transform;
}

std::shared_ptr<VariableTransformation> Approximation::
get_variable_transformation(){
  return varTransform_;
}

int Approximation::num_vars(){
  if (!varTransform_)
    throw(std::runtime_error("Variable transform has not been set"));
  return varTransform_->num_vars();
}

}  // namespace surrogates
}  // namespace Pecos
