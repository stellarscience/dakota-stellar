/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include<BoundedVariables.hpp>

namespace Pecos {
namespace surrogates {

BoundedVariables::BoundedVariables(){}

BoundedVariables::~BoundedVariables(){}

void BoundedVariables::set_options(const util::OptionsList &opts){

  //set_ranges(ranges);
}

void BoundedVariables::set_ranges(const RealVector &ranges){
  ranges_.sizeUninitialized(ranges.length());
  ranges_.assign(ranges);
  set_num_vars(ranges.length()/2);
}

Real BoundedVariables::lb(int i) const{
  return ranges_[2*i];
}

Real BoundedVariables::ub(int i) const{
  return ranges_[2*i+1];
}

void define_homogeneous_ranges(int num_vars, Real lb, Real ub,
			       RealVector &ranges){
  ranges.sizeUninitialized(2*num_vars);
  for (int i=0; i<num_vars; ++i){
    ranges[2*i] = lb;
    ranges[2*i+1] = ub;
  }
}

}  // namespace surrogates
}  // namespace Pecos
