/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <VariableTransformation.hpp>

namespace Pecos {
namespace surrogates {

VariableTransformation::VariableTransformation(){}

VariableTransformation::~VariableTransformation(){}

void VariableTransformation::
set_variables(const std::shared_ptr<Variables> &vars){
  vars_ = vars;
}

int  VariableTransformation::num_vars(){
  return vars_->num_vars();
}

}  // namespace surrogates
}  // namespace Pecos
