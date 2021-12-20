/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <CppFunction.hpp>


namespace Pecos {
namespace surrogates {

CppFunction::CppFunction(){}

CppFunction::~CppFunction(){}

void CppFunction::
set_function(void (*function)(const Real* sample, Real* func_vals,
			      const util::OptionsList &opts)){
  targetFunction_ = function;
}

void CppFunction::value(const RealVector &sample, RealVector &values){
  util::resize_if_needed(values,numQOI_);
  targetFunction_(sample.values(),values.values(),opts_);
}

void CppFunction::value(const RealMatrix &samples, RealMatrix &values){
  util::resize_if_needed(values,samples.numCols(),numQOI_);
  RealVector val(numQOI_,false); //zeroOut = false
  for (int i=0; i<samples.numCols(); ++i){
    const RealVector sample = Teuchos::getCol(Teuchos::View,
					      const_cast<RealMatrix&>(samples),i);
    value(sample,val);
    for(int j = 0; j < numQOI_; j++)
      values(i,j) = val(j);
  }
}

}  // namespace surrogates
}  // namespace Pecos
