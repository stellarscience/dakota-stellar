/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <Function.hpp>

using Pecos::util::OptionsList;

namespace Pecos {
namespace surrogates {

Function::Function() : numQOI_(0) {}

Function::~Function(){}

void Function::set_options(const OptionsList &opts) {
  opts_=opts;
  numQOI_ = opts_.get<int>("num_qoi");
}

void Function::get_options(OptionsList &opts) {
  opts=opts_;
}

void Function::gradient(const RealMatrix &samples, int qoi, RealMatrix &gradients) {
  throw(std::string("This Function type does not provide gradients"));
}

void Function::jacobian(const RealVector &sample, RealMatrix &jacobian) {
  throw(std::string("This Function type does not provide a jacobian"));
}

void Function::hessian(const RealMatrix &samples, int qoi, RealMatrixList &hessians) {
  throw(std::string("This Function type does not provide gradients"));
}

}  // namespace surrogates
}  // namespace Pecos
