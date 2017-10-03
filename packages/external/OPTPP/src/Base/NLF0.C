//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cfloat>
#include <cstring>
#else
#include <float.h>
#include <string.h>
#endif

#include "NLF.h"

using namespace std;
using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;


namespace OPTPP {

//-------------------------------------------------------------------------
// NLF0 method routines.
//-------------------------------------------------------------------------
void NLF0::reset() // Reset parameter values 
{
  init_flag = false;  
  nfevals   = 0;
#ifdef OPTPP_HAVE_MPI
  SpecFlag = Spec1;
#else
  SpecFlag = NoSpec;
#endif
  application.reset();
}

void NLF0::initFcn() // Initialize Function
{
  if (init_flag == false)  {
    init_fcn(dim, mem_xc);
    init_flag = true;
  }
  else  {
    cerr << "NLF0:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc);
  }
}

double NLF0::evalF() // Evaluate Function
{
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(mem_xc,fvalue)) {
      fcn_v(dim, mem_xc, fvalue, result, vptr);
      application.update(NLPFunction,dim,mem_xc,fvalue);
      nfevals++;
    }
  }
  else {
    SpecFlag = Spec1;
    (void) evalG();
    SpecFlag = Spec2;
  }

  function_time = get_wall_clock_time() - time0;
  return fvalue;
}

double NLF0::evalF(const SerialDenseVector<int,double>& x) // Evaluate Function at x
{
  double fx;
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(x,fx)) {
      fcn_v(dim, x, fx, result, vptr);
      application.update(NLPFunction,dim,x,fx);
      nfevals++;
    }
  }
  else {
    SpecFlag = Spec1;
    (void) evalG(x);
    fx = specF;
    SpecFlag = Spec2;
  }

  function_time = get_wall_clock_time() - time0;
  return fx;
}

SerialDenseVector<int,double> NLF0::evalG() 
{
  SerialDenseVector<int,double> grad(dim);
  SerialDenseVector<int,double> sx(dim);
  sx = 1.0;

  // Since NLF0 objects do not have analytic gradients supply
  // one by using finite differences

  grad = FDGrad(sx, mem_xc, fvalue, partial_grad);
  return grad;
}

SerialDenseVector<int,double> NLF0::evalG(const SerialDenseVector<int,double>& x) 
{
  SerialDenseVector<int,double> gx(dim);
  SerialDenseVector<int,double> sx(dim);
  sx = 1.0;

  // Since NLF0 objects do not have analytic gradients supply
  // one by using finite differences

  if (SpecFlag == NoSpec) {
    int result = 0;
    if (!application.getF(x, specF)) {
      fcn_v(dim, x, specF, result, vptr);
      nfevals++;
    }
  }

  gx = FDGrad(sx, x, specF, partial_grad);
  return gx;
}

SerialSymDenseMatrix<int,double> NLF0::evalH() 
{
// Since NLF0 objects do not have analytic hessians supply
// one by using finite differences
  std::cout<<"NLF0.C"<<std::endl;
  SerialSymDenseMatrix<int,double> hess(dim);

  hess = FD2Hessian(mem_xc);
  return hess;
}

SerialSymDenseMatrix<int,double> NLF0::evalH(SerialDenseVector<int,double>& x) 
{
// Since NLF0 objects do not have analytic hessians supply
// one by using finite differences

  SerialSymDenseMatrix<int,double> hess(dim);

  hess = FD2Hessian(x);
  return hess;
}

void NLF0::eval()
{
  (void) evalF();
}

real NLF0::evalLagrangian(const SerialDenseVector<int,double>& xc , 
                          SerialDenseVector<int,double>& multiplier,
                          const SerialDenseVector<int,double>& type)
{
   real result = evalF(xc);
   if( hasConstraints()){
     SerialDenseVector<int,double> resid(constraint_->evalResidual(xc));
     //     SerialDenseVector<int,double> resid(constraint_->getNumOfCons());
     //     resid = constraint_->evalResidual(xc);
      result  -= resid.dot(multiplier);
   }
   return result;
}

SerialDenseVector<int,double> NLF0::evalLagrangianGradient(const SerialDenseVector<int,double>& xc, 
                                          const SerialDenseVector<int,double>& multiplier,
					  const SerialDenseVector<int,double>& type) 
{
   SerialDenseVector<int,double> grad  = evalG(xc);
   SerialDenseVector<int,double> gradtmp(grad.length());
   //   SerialDenseVector<int,double> gradtmp(xc.length());
   if(hasConstraints()){
      SerialDenseVector<int,double> tmult = multiplier;
      for (int i = 0; i < getNumOfCons(); i++){
         if(type(i) == NLineq || type(i) == Lineq)
            tmult(i)*= -1;
      }
      gradtmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,1.0,constraint_->evalGradient(xc),tmult,0.0);     
      // grad += constraint_->evalGradient(xc)*tmult;
      grad += gradtmp;
      
   }
   return grad;
}


SerialDenseVector<int,double> NLF0::evalCF(const SerialDenseVector<int,double>& x) // Evaluate Nonlinear Constraint at x
{
  SerialDenseVector<int,double> cfx(ncnln);
  int result = 0;

  double time0 = get_wall_clock_time();


  // *** CHANGE *** //
  if (!application.getCF(x,cfx)) {
    confcn(dim, x, cfx, result);
    application.constraint_update(NLPFunction,dim,ncnln,x,cfx);
   // nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;
  
  setConstraintValue(cfx);
  return cfx;
}

SerialDenseMatrix<int,double> NLF0::evalCG(const SerialDenseVector<int,double>& x) // Evaluate Nonlinear Constraint Gradient at x
{
// Since NLF0 objects do not have analytic gradients supply
// one by using finite differences

  SerialDenseMatrix<int,double> grad(dim,ncnln);
  grad = CONFDGrad(x);
  return grad;
}

SerialSymDenseMatrix<int,double> NLF0::evalCH(SerialDenseVector<int,double>& x) 
{
// CPJW - This is a placeholder routine.  NIPS is the only algorithm which
// supports nonlinear constraints and currently this routine is never accessed.
// The true evaluator will be implemented later.
// Since NLF0 objects do not have analytic hessians supply
// one by using finite differences

  SerialSymDenseMatrix<int,double> hess(dim);
  hess = 0.0;
  return hess;
}

OptppArray<SerialSymDenseMatrix<int,double> > NLF0::evalCH(SerialDenseVector<int,double>& x, int darg) 
{
// CPJW - This is a placeholder routine.  NIPS is the only algorithm which
// supports nonlinear constraints and currently this routine is never accessed.
// The true evaluator will be implemented later.
// Since NLF0 objects do not have analytic hessians supply
// one by using finite differences

  OptppArray<SerialSymDenseMatrix<int,double> > hess(1);
  SerialSymDenseMatrix<int,double> hessT(dim);
  hessT = 0.0;
  hess[0] = hessT;
  return hess;
}

void NLF0::evalC(const SerialDenseVector<int,double>& x)
{
  (void) evalCF(x);
}

} // namespace OPTPP
