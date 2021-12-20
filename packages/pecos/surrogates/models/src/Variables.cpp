/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <Variables.hpp>

namespace Pecos {
namespace surrogates {

Variables::Variables() : numVars_(0){}

Variables::~Variables(){}

void Variables::set_options(const util::OptionsList &opts){
}

void Variables::get_options(util::OptionsList &opts){
}

void Variables::set_num_vars(int num_vars){
  numVars_ = num_vars;
}

int Variables::num_vars(){return numVars_;}

}  // namespace surrogates
}  // namespace Pecos
