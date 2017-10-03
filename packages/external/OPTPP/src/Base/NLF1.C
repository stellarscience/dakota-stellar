//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include "NLF.h"
#include "TOLS.h"
#include "cblas.h"

using namespace std;
using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

namespace OPTPP {

void NLF1::reset() // Reset parameter values  
{
  init_flag = false;    
  nfevals   = ngevals = 0; 
#ifdef OPTPP_HAVE_MPI
  SpecFlag = Spec1;
#else
  SpecFlag = NoSpec;
#endif
  application.reset();
}

void NLF1::initFcn() // Initialize Function
{
  if (init_flag == false)  {
      init_fcn(dim, mem_xc);
      init_flag = true;
  }
  else  {
    cerr << "NLF1:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc);
  }
}

double NLF1::evalF() // Evaluate Function
{
  int result = 0;
  SerialDenseVector<int,double> gtmp(dim);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(mem_xc,fvalue)) {
    fcn_v(NLPFunction, dim, mem_xc, fvalue, gtmp, result, vptr);
    application.update(result,dim,mem_xc,fvalue,gtmp);
    nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_)
  cout << "NLF1::evalF()\n" 
    << "nfevals       = " << nfevals << "\n"
    << "fvalue        = " << fvalue << "\n"
    << "function time = " << function_time << "\n";
  return fvalue;
}

double NLF1::evalF(const SerialDenseVector<int,double>& x) // Evaluate Function at x
{
  int    result = 0;
  double fx;
  SerialDenseVector<int,double> gtmp(dim);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(x,fx)) {
    fcn_v(NLPFunction, dim, x, fx, gtmp, result, vptr);
    application.update(result,dim,x,fx,gtmp);
    nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_)
  cout << "NLF1::evalF(x)\n" 
    << "nfevals       = " << nfevals << "\n"
    << "fvalue        = " << fx << "\n"
    << "function time = " << function_time << "\n";
  return fx;
}

SerialDenseVector<int,double> NLF1::evalG() // Evaluate the gradient
{
  int    result = 0;
  double fx;

  // *** CHANGE *** //
  if (!application.getGrad(mem_xc,mem_grad)) {
    fcn_v(NLPGradient, dim, mem_xc, fx, mem_grad, result, vptr);
    application.update(result,dim,mem_xc,fx,mem_grad);
    ngevals++;
  }
  // *** CHANGE *** //
  return mem_grad;
}

SerialDenseVector<int,double> NLF1::evalG(const SerialDenseVector<int,double>& x) // Evaluate the gradient at x
{
  int    result = 0 ;
  double fx;
  SerialDenseVector<int,double> gx(dim);

  // *** CHANGE *** //
  if (!application.getGrad(x,gx)) {
    fcn_v(NLPGradient, dim, x, fx, gx, result, vptr);
    application.update(result,dim,x,fx,gx);
    ngevals++;
  }
  // *** CHANGE *** //
  return gx;
}

SerialSymDenseMatrix<int,double> NLF1::evalH() // Evaluate the Hessian
{
  SerialDenseVector<int,double> sx(dim);
  SerialSymDenseMatrix<int,double> Hessian(dim);

  sx = 1.0;
  Hessian = FDHessian(sx);
  return Hessian;
}

SerialSymDenseMatrix<int,double> NLF1::evalH(SerialDenseVector<int,double>& x) // Evaluate the Hessian at x
{
  SerialSymDenseMatrix<int,double> Hessian(dim);

  Hessian = FDHessian(x);
  return Hessian;
}

void NLF1::eval() // Evaluate Function and Gradient
{
  int mode = NLPFunction | NLPGradient, result = 0;

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(mem_xc,fvalue) || !application.getGrad(mem_xc,mem_grad)) { 
    fcn_v(mode, dim, mem_xc, fvalue, mem_grad, result, vptr);
    application.update(result,dim,mem_xc,fvalue,mem_grad);
    nfevals++; ngevals++;
  }
  // *** CHANGE *** //

  function_time = get_wall_clock_time() - time0;

  if (debug_)
  cout << "NLF1::eval()\n" 
    << "mode          = " << mode   << "\n"
    << "nfevals       = " << nfevals << "\n"
    << "fvalue        = " << fvalue << "\n"
    << "function time = " << function_time << "\n";
}

double NLF1::evalLagrangian(const SerialDenseVector<int,double>& xc , 
                          SerialDenseVector<int,double>& multiplier,
                          const SerialDenseVector<int,double>& type) 
{
   double result = evalF(xc);
   if( hasConstraints()){
     SerialDenseVector<int,double> resid(constraint_->evalResidual(xc));
     //     SerialDenseVector<int,double> resid(constraint_->getNumOfCons());
     //     resid = constraint_->evalResidual(xc);
      result  -=  resid.dot(multiplier);
   }
   return result;
}

SerialDenseVector<int,double> NLF1::evalLagrangianGradient(const SerialDenseVector<int,double>& xc, 
                                          const SerialDenseVector<int,double>& multiplier,
					  const SerialDenseVector<int,double>& type) 
{
   SerialDenseVector<int,double> grad  = evalG(xc);
   if(hasConstraints()){
     SerialDenseVector<int,double> tmult(multiplier);
     //     tmult = multiplier;
      int alpha = grad.length();
      SerialDenseVector<int,double> tmult2(alpha);
      //tmult *= -1;
      // grad += constraint_->evalGradient(xc)*tmult;
      tmult2.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,-1.0,constraint_->evalGradient(xc),tmult,0.0);
      grad+=tmult2;
   }
  return grad;

}


SerialDenseVector<int,double> NLF1::evalCF(const SerialDenseVector<int,double>& x) // Evaluate Function at x
{
  int    result = 0;
  SerialDenseVector<int,double> cfx(ncnln);
  SerialDenseMatrix<int,double> gtmp(dim,ncnln);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getCF(x,cfx)) {
    confcn(NLPFunction, dim, x, cfx, gtmp, result);
    application.constraint_update(result,dim,ncnln,x,cfx,gtmp);
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_)
    cout << "NLF1::evalCF(x)\n" 
         << "nfevals       = " << nfevals << "\n"
       //  << "fvalue        = " << cfx << "\n"
         << "function time = " << function_time << "\n";
  return cfx;
}

SerialDenseMatrix<int,double> NLF1::evalCG(const SerialDenseVector<int,double>& x) // Evaluate the gradient at x
{
  int    result = 0 ;
  SerialDenseVector<int,double> cfx(ncnln);
  SerialDenseMatrix<int,double> cgx(dim,ncnln);

  // *** CHANGE *** //
  if (!application.getCGrad(x,cgx)) {
    confcn(NLPGradient, dim, x, cfx, cgx, result);
    application.constraint_update(result,dim,ncnln,x,cfx,cgx);
  }
  // *** CHANGE *** //
  return cgx;
}

SerialSymDenseMatrix<int,double> NLF1::evalCH(SerialDenseVector<int,double>& x) // Evaluate the Hessian at x
{
  SerialSymDenseMatrix<int,double> Hessian(dim);

  // PJW This is a dummy routine.  NIPS is the only algorithm which supports
  // nonlinear constraints and currently this routine is not being accessed.
  Hessian = 0.0;
  return Hessian;
}


OptppArray<SerialSymDenseMatrix<int,double> > NLF1::evalCH(SerialDenseVector<int,double>& x, int darg) // Evaluate the Hessian at x
{
  OptppArray<SerialSymDenseMatrix<int,double> > Hessian(ncnln);
  Hessian = CONFDHessian(x);
  return Hessian;
}

void NLF1::evalC(const SerialDenseVector<int,double>& x)
{
  int mode = NLPFunction | NLPGradient, result = 0;
  SerialDenseVector<int,double> cfx(ncnln);
  SerialDenseMatrix<int,double> cgx(dim, ncnln);

  double time0 = get_wall_clock_time();

  if (!application.getCF(x, cfx) || !application.getCGrad(x, cgx)) {
    confcn(mode, dim, x, cfx, cgx, result);
    application.constraint_update(result, dim, ncnln, x, cfx, cgx);
  }

  function_time = get_wall_clock_time() - time0;
}

} // namespace OPTPP
