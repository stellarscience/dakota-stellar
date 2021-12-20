// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_PD_MEANSEMIDEVIATIONFROMTARGET_HPP
#define ROL_PD_MEANSEMIDEVIATIONFROMTARGET_HPP

#include "ROL_PD_RandVarFunctional.hpp"

namespace ROL {

template<class Real>
class PD_MeanSemiDeviationFromTarget : public PD_RandVarFunctional<Real> {
private:
  Real coeff_;
  Real target_;

  Ptr<SampledScalar<Real>> values_;
  Ptr<SampledScalar<Real>> gradvecs_;
  Ptr<SampledVector<Real>> gradients_;
  Ptr<SampledVector<Real>> hessvecs_;

  using RandVarFunctional<Real>::val_;
  using RandVarFunctional<Real>::g_;
  using RandVarFunctional<Real>::gv_;
  using RandVarFunctional<Real>::hv_;
  using RandVarFunctional<Real>::dualVector_;

  using RandVarFunctional<Real>::point_;
  using RandVarFunctional<Real>::weight_;

  using RandVarFunctional<Real>::computeValue;
  using RandVarFunctional<Real>::computeGradient;
  using RandVarFunctional<Real>::computeGradVec;
  using RandVarFunctional<Real>::computeHessVec;

  using PD_RandVarFunctional<Real>::setValue;
  using PD_RandVarFunctional<Real>::getMultiplier;
  using PD_RandVarFunctional<Real>::getPenaltyParameter;
  using PD_RandVarFunctional<Real>::ppf;

  void initializeStorage(void) {
    values_    = makePtr<SampledScalar<Real>>();
    gradvecs_  = makePtr<SampledScalar<Real>>();
    gradients_ = makePtr<SampledVector<Real>>();
    hessvecs_  = makePtr<SampledVector<Real>>();

    RandVarFunctional<Real>::setStorage(values_,gradients_);
    RandVarFunctional<Real>::setHessVecStorage(gradvecs_,hessvecs_);
  }

  void clear(void) {
    gradvecs_->update();
    hessvecs_->update();
  }

  void checkInputs(void) {
    Real zero(0);
    ROL_TEST_FOR_EXCEPTION((coeff_ < zero), std::invalid_argument,
      ">>> ERROR (ROL::PD_MeanSemiDeviation): Element of coefficient array out of range!");
    initializeStorage();
  }

public:
  PD_MeanSemiDeviationFromTarget(const Real coeff, const Real target)
    : PD_RandVarFunctional<Real>(), coeff_(coeff), target_(target) {
    checkInputs();
  }

  void setStorage(const Ptr<SampledScalar<Real>> &value_storage,
                  const Ptr<SampledVector<Real>> &gradient_storage) {
    values_    = value_storage;
    gradients_ = gradient_storage;
    PD_RandVarFunctional<Real>::setStorage(values_,gradients_);
  }

  void setHessVecStorage(const Ptr<SampledScalar<Real>> &gradvec_storage,
                         const Ptr<SampledVector<Real>> &hessvec_storage) {
    gradvecs_ = gradvec_storage;
    hessvecs_ = hessvec_storage;
    PD_RandVarFunctional<Real>::setHessVecStorage(gradvecs_,hessvecs_);
  }
 
  void initialize(const Vector<Real> &x) {
    PD_RandVarFunctional<Real>::initialize(x);
    clear();
  }
  
  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real lam(0);
    getMultiplier(lam, point_);
    Real val = computeValue(obj,x,tol);
    Real arg = coeff_ * (val - target_);
    Real  pf = ppf(arg, lam, getPenaltyParameter(), 0);
    val_ += weight_ * (val + pf);
    setValue(arg, point_);
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real ev(0);
    sampler.sumAll(&val_,&ev,1);
    return ev;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    const Real one(1);
    Real lam(0);
    getMultiplier(lam, point_);
    Real val = computeValue(obj,x,tol);
    Real arg = coeff_ * (val - target_);
    Real  pf = ppf(arg, lam, getPenaltyParameter(), 1);
    computeGradient(*dualVector_,obj,x,tol);
    g_->axpy(weight_ * (one + coeff_ * pf), *dualVector_);
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    sampler.sumAll(*g_,g);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    const Real zero(0), one(1);
    Real lam(0);
    getMultiplier(lam, point_);
    Real val = computeValue(obj,x,tol);
    Real arg = coeff_ * (val - target_);
    Real pf1 = ppf(arg, lam, getPenaltyParameter(), 1);
    Real pf2 = ppf(arg, lam, getPenaltyParameter(), 2);
    computeHessVec(*dualVector_, obj, v, x, tol);
    hv_->axpy(weight_ * (one + pf1 * coeff_), *dualVector_); 
    if ( pf2 > zero ) {
      Real gv = computeGradVec(*dualVector_, obj, v, x, tol);
      hv_->axpy(weight_ * pf2 * coeff_ * coeff_ * gv, *dualVector_);
    }
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    sampler.sumAll(*hv_,hv);
  }
};

}

#endif
