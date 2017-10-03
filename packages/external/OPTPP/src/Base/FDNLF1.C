//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include "NLF.h"
#include "TOLS.h"
#include "ioformat.h"
#include "cblas.h"

using namespace std;
using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;


namespace OPTPP {

//-------------------------------------------------------------------------
// FDNLF1 method routines.
//-------------------------------------------------------------------------
void FDNLF1::reset() // Reset parameter values  
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

void FDNLF1::initFcn() // Initialize Function
{
  if (init_flag == false) {
      init_fcn(dim, mem_xc);
      init_flag = true;
  }
  else  {
    cerr << "FDNLF1:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc);
  }
}

void FDNLF1::eval() // Evaluate Function and Gradient
{
  (void)evalF();
  (void)evalG();
}

double FDNLF1::evalF() // Evaluate Function
{
  int result = 0;
  double time0 = get_wall_clock_time();


  if (SpecFlag == NoSpec) {
    if (!application.getF(mem_xc, fvalue)) {
      fcn_v(dim, mem_xc, fvalue, result, vptr);
      function_time = get_wall_clock_time() - time0;
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

double FDNLF1::evalF(const SerialDenseVector<int,double>& x) // Evaluate Function at x

{
  double fx;
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(x, fx)) {
      fcn_v(dim, x, fx, result, vptr);
      function_time = get_wall_clock_time() - time0;
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

SerialDenseVector<int,double> FDNLF1::evalG() // Evaluate the gradient
{ 

  SerialDenseVector<int,double> sx(dim);
  sx = 1.0;
  ngevals++;

  if (finitediff == ForwardDiff)
    mem_grad =  FDGrad(sx, mem_xc, fvalue, partial_grad);
  else if (finitediff == BackwardDiff)
    mem_grad =  BDGrad(sx, mem_xc, fvalue, partial_grad);
  else if (finitediff == CentralDiff)
    mem_grad =  CDGrad(sx, mem_xc, fvalue, partial_grad);
  else {
    cout << "FDNLF1::evalG: Unrecognized difference option\n";
    cout << "FDNLF1::evalG: Using forward difference option\n";
    mem_grad =  FDGrad(sx, mem_xc, fvalue, partial_grad);
  }
 return mem_grad;
}

SerialDenseVector<int,double> FDNLF1::evalG(const SerialDenseVector<int,double>& x) // Evaluate the gradient at x
{
  SerialDenseVector<int,double> gx(dim);
  SerialDenseVector<int,double> sx(dim);
  sx = 1.0;
  ngevals++;

  if (SpecFlag == NoSpec) {
    int result = 0;
    if (!application.getF(x, specF)) {
      fcn_v(dim, x, specF, result, vptr);
      nfevals++;
    }
  }

  if (finitediff == ForwardDiff) 
    gx = FDGrad(sx, x, specF, partial_grad);
  else if (finitediff == BackwardDiff)
    gx = BDGrad(sx, x, specF, partial_grad);
  else if (finitediff == CentralDiff)
    gx = CDGrad(sx, x, specF, partial_grad);
  else {
    cout << "FDNLF1::evalG: Unrecognized difference option\n";
    cout << "FDNLF1::evalG: Using forward difference option\n";
    mem_grad =  FDGrad(sx, x, specF, partial_grad);
  }
  return gx;
}

SerialSymDenseMatrix<int,double> FDNLF1::FDHessian(SerialDenseVector<int,double>& sx) // Evaluate the Hessian
{
  SerialSymDenseMatrix<int,double> Hessian(dim);

  sx = 1.0;
  Hessian = FD2Hessian(sx);
  return Hessian;
}

SerialSymDenseMatrix<int,double> FDNLF1::evalH() // Evaluate the Hessian
{
  SerialDenseVector<int,double> sx(dim);
  SerialSymDenseMatrix<int,double> Hessian(dim);

  sx = 1.0;
  Hessian = FD2Hessian(sx);
  return Hessian;
}

SerialSymDenseMatrix<int,double> FDNLF1::evalH(SerialDenseVector<int,double>& x) // Evaluate the Hessian
{
  SerialSymDenseMatrix<int,double> Hessian(dim);

  Hessian = FD2Hessian(x);
  return Hessian;
}

void FDNLF1::printState(const char * s) 
{ // Print out current state: x current, gradient and Function value
  cout <<"\n\n=========  " << s << "  ===========\n\n";
  cout <<"\n   i\t    xc \t\t grad \t\t fcn_accrcy \n";
  for (int i=0; i<dim; i++) 
    cout << d(i,6) << e(mem_xc(i),12,4)<< "\t" << e(mem_grad(i),12,4) 
         << "\t"   << e(mem_fcn_accrcy(i),12,4) << "\n";
  cout <<"\nFunction Value     = " << e(fvalue,12,4) << "\n";
  double gnorm = sqrt(mem_grad.dot(mem_grad));
  cout <<"Norm of gradient   = " << e(gnorm,12,4) << "\n";
  cout <<"Derivative Option  = " << finitediff << "\n\n";
}

void FDNLF1::fPrintState(ostream *nlpout, const char * s) 
{ // Print out current state: x current, gradient and Function value
  (*nlpout) <<"\n\n=========  " << s << "  ===========\n\n";
  (*nlpout) <<"\n   i\t    xc \t\t grad \t\t fcn_accrcy \n";
  for (int i=0; i<dim; i++) 
    (*nlpout) << d(i,6) << e(mem_xc(i),12,4)<< "\t" << e(mem_grad(i),12,4) 
         << "\t"   << e(mem_fcn_accrcy(i),12,4) << "\n";
  (*nlpout) <<"\nFunction Value     = " << e(fvalue,12,4) << "\n";
  double gnorm = sqrt(mem_grad.dot(mem_grad));
  (*nlpout) <<"Norm of gradient   = " << e(gnorm,12,4) << "\n";
//  (*nlpout) <<"Function Accuracy  = " << e(mem_fcn_accrcy,12,4) << "\n";
  (*nlpout) <<"Derivative Option  = " << finitediff << "\n\n";
}


double FDNLF1::evalLagrangian(const SerialDenseVector<int,double>& xc , 
                            SerialDenseVector<int,double>& multiplier,
                            const SerialDenseVector<int,double>& type) 
{
  double result = evalF(xc);
   if( hasConstraints()){
     //      SerialDenseVector<int,double> resid = constraint_->evalResidual(xc);
     SerialDenseVector<int,double> resid(constraint_->evalResidual(xc));
      result  -=  resid.dot(multiplier);
   }
   return result;
}

SerialDenseVector<int,double> FDNLF1::evalLagrangianGradient(const SerialDenseVector<int,double>& xc, 
                                          const SerialDenseVector<int,double>& multiplier,
                                          const SerialDenseVector<int,double>& type) 
{
   mem_grad  = evalG(xc);
   SerialDenseVector<int,double> grad(mem_grad);
   //   SerialDenseVector<int,double> grad(mem_grad.length());
   //   grad = mem_grad;
   SerialDenseVector<int,double> temp(grad.length());
   //   SerialDenseVector<int,double> temp(constraint_->getNumOfCons());
   if(hasConstraints())
     //SerialDenseVector<int,double> temp;
   temp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, constraint_->evalGradient(xc),multiplier, 0.0);
   // grad -= constraint_->evalGradient(xc)*multiplier;
   grad -= temp;   
   return grad;
}

SerialDenseVector<int,double> FDNLF1::evalCF(const SerialDenseVector<int,double>& x) // Evaluate Constraint Fcn at x
{
  int result = 0;
  SerialDenseVector<int,double> cfx(ncnln);
  double time0 = get_wall_clock_time();
  confcn(dim, x, cfx,result);
  function_time = get_wall_clock_time() - time0;

  //nfevals++;
  return cfx;
}

SerialDenseMatrix<int,double> FDNLF1::evalCG(const SerialDenseVector<int,double>& x) // Evaluate the gradient at x
{
  SerialDenseVector<int,double> sx(dim);
  sx = 1.0;
  SerialDenseVector<int,double> xsave(dim);
  SerialDenseMatrix<int,double> gx(dim, ncnln);
  xsave = getXc();
  setX(x);
  if (finitediff == ForwardDiff) 
    gx = CONFDGrad(sx);
  else if (finitediff == BackwardDiff)
    gx = CONBDGrad(sx);
  else if (finitediff == CentralDiff)
    gx = CONCDGrad(sx);
  else {
    cout <<"FDNLF1::evalG: Unrecognized difference option\n";
  }
  
  setX(xsave);
  //ngevals++;
  return gx;
}

SerialSymDenseMatrix<int,double> FDNLF1::evalCH(SerialDenseVector<int,double>& x) // Evaluate the Hessian
{
  SerialSymDenseMatrix<int,double> Hessian(dim);

  Hessian = FD2Hessian(x);
  return Hessian;
}

OptppArray<SerialSymDenseMatrix<int,double> > FDNLF1::evalCH( SerialDenseVector<int,double>& x, int darg) // Evaluate the Hessian
{
  SerialSymDenseMatrix<int,double> Hessian(dim);

  Hessian = FD2Hessian(x);
  OptppArray<SerialSymDenseMatrix<int,double> > H(1);
  H[0] = Hessian;
  return H;
}

void FDNLF1::evalC(const SerialDenseVector<int,double>& x) // Evaluate Function and Gradient
{
  (void) evalCF(x);
  (void) evalCG(x);
}

} // namespace OPTPP
