//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include "NLF.h"
#include "TOLS.h"
#include "cblas.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

using namespace std;
using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

namespace OPTPP {

void NLF2::reset() // Reset parameter values  
{
  init_flag = false;    
  nfevals   = ngevals = nhevals = 0; 
#ifdef OPTPP_HAVE_MPI
  SpecFlag = Spec1;
#else
  SpecFlag = NoSpec;
#endif
  application.reset();
}

void NLF2::initFcn() // Initialize Function
{
  if (init_flag == false) {
    init_fcn(dim, mem_xc);
    init_flag = true;
  }
  else {
    cerr << "NLF2:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc);
  }
}

double NLF2::evalF() // Evaluate Function
{
  int  result = 0;
  SerialDenseVector<int, double> gtmp(dim);
  SerialSymDenseMatrix<int,double> Htmp(dim);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(mem_xc,fvalue)) {
    fcn_v(NLPFunction, dim, mem_xc, fvalue, gtmp, Htmp,result,vptr);
    application.update(result,dim,mem_xc,fvalue,gtmp,Htmp);
    nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_)  cout << "NLF2::evalF()\n" 
                   << "nfevals       = " << nfevals   << "\n"
                   << "fvalue        = " << fvalue << "\n"
                   << "function time = " << function_time << "\n";
  return fvalue;
}

double NLF2::evalF(const SerialDenseVector<int, double>& x) // Evaluate Function at x
{
  int    result = 0;
  double fx;
  SerialDenseVector<int, double> gtmp(dim);
  SerialSymDenseMatrix<int,double> Htmp(dim);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(x,fx)) {
    fcn_v(NLPFunction, dim, x, fx, gtmp, Htmp,result,vptr);
    application.update(result,dim,x,fx,gtmp,Htmp);
    nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_) cout << "NLF2::evalF(x)\n" 
                  << "nfevals       = " << nfevals   << "\n"
                  << "fvalue        = " << fvalue << "\n"
                  << "function time = " << function_time << "\n";
  return fx;
}

SerialDenseVector<int, double> NLF2::evalG() // Evaluate the gradient
{
  int    result = 0;
  double fx;
  SerialSymDenseMatrix<int,double> Htmp(dim);

  // *** CHANGE *** //
  if (!application.getGrad(mem_xc,mem_grad)) {
    fcn_v(NLPGradient, dim, mem_xc, fx, mem_grad, Htmp,result,vptr);
    application.update(result,dim,mem_xc,fx,mem_grad,Htmp);
    ngevals++;
  }
  // *** CHANGE *** //
  return mem_grad;
}

SerialDenseVector<int, double> NLF2::evalG(const SerialDenseVector<int, double>& x) // Evaluate the gradient at x
{
  int    result = 0;
  double fx;
  SerialDenseVector<int, double> gx(dim);
  SerialSymDenseMatrix<int,double> Htmp(dim);

  // *** CHANGE *** //
  if (!application.getGrad(x,gx)) {
    fcn_v(NLPGradient, dim, x, fx, gx, Htmp, result,vptr);
    application.update(result,dim,x,fx,gx,Htmp);
    ngevals++;
  }
  // *** CHANGE *** //
  return gx;
}

SerialSymDenseMatrix<int,double> NLF2::evalH() // Evaluate the Hessian
{
  int    result = 0;
  double fx;
  SerialDenseVector<int, double> gtmp(dim);

  // *** CHANGE *** //
  if (!application.getHess(mem_xc,Hessian)) {
    fcn_v(NLPHessian, dim, mem_xc, fx, gtmp, Hessian, result,vptr);
    application.update(result,dim,mem_xc,fx,gtmp,Hessian);
    nhevals++;
  }
  // *** CHANGE *** //
  return Hessian;
}

SerialSymDenseMatrix<int,double> NLF2::evalH(SerialDenseVector<int, double>& x) // Evaluate the hessian at x
{
  int    result = 0;
  double fx;
  SerialDenseVector<int, double> gx(dim);
  SerialSymDenseMatrix<int,double> Hx(dim);

  // *** CHANGE *** //
  if (!application.getHess(x,Hx)) {
    fcn_v(NLPHessian, dim, x, fx, gx, Hx, result,vptr);
    application.update(result,dim,x,fx,gx,Hx);
    nhevals++;
  }
  // *** CHANGE *** //
  return Hx;
}

void NLF2::eval() // Evaluate Function, Gradient, and Hessian
{
  int mode = NLPFunction | NLPGradient | NLPHessian, result = 0;

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(mem_xc,fvalue) || !application.getGrad(mem_xc,mem_grad) ||
      !application.getHess(mem_xc,Hessian)) { 
    fcn_v(mode, dim, mem_xc, fvalue, mem_grad, Hessian,result,vptr);
    application.update(result,dim,mem_xc,fvalue,mem_grad,Hessian);
    nfevals++; ngevals++; nhevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_)  cout << "NLF2::eval()\n" 
                   << "mode          = " << mode   << "\n"
                   << "nfevals       = " << nfevals   << "\n"
	           << "fvalue        = " << fvalue << "\n"
	           << "function time = " << function_time << "\n";
  
}

double NLF2::evalLagrangian(const SerialDenseVector<int, double>& xc , 
                          SerialDenseVector<int, double>& multiplier,
                          const SerialDenseVector<int, double>& type) 
{
   double result = evalF(xc);
   if( hasConstraints()){
     SerialDenseVector<int, double> resid(constraint_->evalResidual(xc));
     //     SerialDenseVector<int, double> resid(constraint_->getNumOfCons());
     //     resid = constraint_->evalResidual(xc);
      result  -=  resid.dot(multiplier);
   }
   return result;
}

SerialDenseVector<int, double> NLF2::evalLagrangianGradient(const SerialDenseVector<int, double>& xc, 
                                          const SerialDenseVector<int, double>& multiplier,
					  const SerialDenseVector<int, double>& type) 
{
  //  SerialDenseVector<int,double> tmult2(xc.length());
  SerialDenseVector<int, double> grad(evalG(xc));
  SerialDenseVector<int,double> tmult2(grad.length());
  //  SerialDenseVector<int, double> grad(xc.length());
  //  grad  = evalG(xc);
   if(hasConstraints()){
     SerialDenseVector<int, double> tmult(multiplier);
     //     SerialDenseVector<int, double> tmult(multiplier.length());
     //     tmult = multiplier;
    tmult *= -1;
    tmult2.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, constraint_->evalGradient(xc),tmult, 0.0);    
// tmult *= constraint_->evalGradient(xc);
    grad += tmult2;
   }
   return grad;
}

SerialSymDenseMatrix<int,double> NLF2::evalLagrangianHessian(SerialDenseVector<int, double>& xc,
                                            const SerialDenseVector<int, double>& multiplier,
                                            const SerialDenseVector<int, double>& type) 
{
  SerialSymDenseMatrix<int,double> hessian(evalH(xc));
  //  SerialSymDenseMatrix<int,double> hessian(xc.length());
  //  hessian = evalH(xc);
   if(hasConstraints()){
     SerialSymDenseMatrix<int,double> thessian(xc.length());
	 thessian = constraint_->evalHessian(xc);
       Print(thessian);
   }
   return hessian;
}

SerialDenseVector<int, double> NLF2::evalCF(const SerialDenseVector<int, double>& x) // Evaluate Function at x
{
  int    result = 0;
  SerialDenseVector<int, double> cfx(ncnln);
  SerialDenseMatrix<int,double> gtmp(dim,ncnln);
  OptppArray<SerialSymDenseMatrix<int,double> > Htmp(ncnln);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getCF(x,cfx)) {

    if(confcn1 != NULL){   
       confcn1(NLPFunction, dim, x, cfx, gtmp, result);
       application.constraint_update(result,dim,ncnln,x,cfx,gtmp);
    }
    else if(confcn2 != NULL){   
       confcn2(NLPFunction, dim, x, cfx, gtmp, Htmp,result);
       application.constraint_update(result,dim,ncnln,x,cfx,gtmp,Htmp);
    }
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_) cout << "NLF2::evalCF(x)\n" 
                  << "nfevals       = " << nfevals   << "\n"
                  << "fvalue(1)        = " << cfx(1) << "\n"
                  << "function time = " << function_time << "\n";
  return cfx;
}

SerialDenseMatrix<int,double> NLF2::evalCG(const SerialDenseVector<int, double>& x) // Evaluate the gradient at x
{
  int    result = 0;
  SerialDenseVector<int, double> cfx(ncnln);
  SerialDenseMatrix<int,double> cgx(dim,ncnln);
  OptppArray<SerialSymDenseMatrix<int,double> > Htmp(ncnln);

  // *** CHANGE *** //
  if (!application.getCGrad(x,cgx)) {
    if(confcn1 != NULL){
      confcn1(NLPGradient, dim, x, cfx, cgx, result);
      application.constraint_update(result,dim,ncnln,x,cfx,cgx);
    }
    if(confcn2 != NULL){
      confcn2(NLPGradient, dim, x, cfx, cgx, Htmp, result);
      application.constraint_update(result,dim,ncnln,x,cfx,cgx,Htmp);
    }
  }
  // *** CHANGE *** //
  return cgx;
}

SerialSymDenseMatrix<int,double> NLF2::evalCH(SerialDenseVector<int, double>& x) // Evaluate the hessian at x
{
  SerialDenseVector<int, double> cfx(ncnln);
  SerialDenseMatrix<int,double> cgx(dim,ncnln);
  SerialSymDenseMatrix<int,double> cHx(dim);

  cHx = 0;
  return cHx;
}
OptppArray<SerialSymDenseMatrix<int,double> > NLF2::evalCH(SerialDenseVector<int, double>& x, int darg) // Evaluate the hessian at x
{
  int    result = 0;
  SerialDenseVector<int, double> cfx(ncnln);
  SerialDenseMatrix<int,double> cgx(dim,ncnln);
  OptppArray<SerialSymDenseMatrix<int,double> > cHx(ncnln);

  // *** CHANGE *** //
  if (!application.getCHess(x,cHx)) {
    if(confcn2 != NULL){
       confcn2(NLPHessian, dim, x, cfx, cgx, cHx, result);
       application.constraint_update(result,dim,ncnln,x,cfx,cgx,cHx);
       nhevals++;
    }
  }
  // *** CHANGE *** //
  return cHx;
}

void NLF2::evalC(const SerialDenseVector<int, double>& x)
{
  int mode1 = NLPFunction | NLPGradient;
  int mode2 = NLPFunction | NLPGradient | NLPHessian;
  int result = 0;
  SerialDenseVector<int, double> cfx(ncnln);
  SerialDenseMatrix<int,double> cgx(dim,ncnln);
  OptppArray<SerialSymDenseMatrix<int,double> > cHx(ncnln);

  double time0 = get_wall_clock_time();

  // *** CHANGE *** //
  if (!application.getCF(x, cfx) || !application.getCGrad(x, cgx) || !application.getCHess(x, cHx)) {
    if(confcn1 != NULL){
      confcn1(mode1, dim, x, cfx, cgx, result);
      application.constraint_update(result, dim, ncnln, x, cfx, cgx);
    }
    if(confcn2 != NULL){
       confcn2(mode2, dim, x, cfx, cgx, cHx, result);
       application.constraint_update(result, dim, ncnln, x, cfx, cgx, cHx);
       nhevals++;
    }
  }

  function_time = get_wall_clock_time() - time0;
}

} // namespace OPTPP
