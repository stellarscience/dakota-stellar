// TEST file for derivatives
//
#include <iostream>
#include <fstream>
#include "NLF.h"
#include "tstfcn.h"
#include "ioformat.h"

using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;

using namespace OPTPP;

int main ()
{
  int i;
  int n = 2;
  ColumnVector x(n), grad(n), gfd(n), gbd(n), gcd(n);
  ColumnVector sx(n);
  SymmetricMatrix H(n);
  //  Create a Nonlinear problem object
  double f, accrcy = 1.0e-9;

  NLF1 nlp(n,rosen,init_rosen);
  
  //  Initialize and evaluate the function at x

  nlp.initFcn();
  for(i = 1; i <=n; i++) nlp.setFcnAccrcy(i,accrcy);
  nlp.eval();
  x = nlp.getXc();
  f = nlp.getF();

  cout <<"*****  Testing Finite-Difference Gradients for NLF1 class *****\n";
  cout <<" i \t analytic grad \t forward \t backward \t central\n";
  sx = 1.0;
  nlp.printState("NLF1");
  grad = nlp.evalG();
  gfd  = nlp.FDGrad(sx,x,f,grad);
  gbd  = nlp.BDGrad(sx,x,f,grad);
  gcd  = nlp.CDGrad(sx,x,f,grad);
  for (i=1; i<=n; i++) {
    cout << i << e(grad(i),12,4)  
              << e(gfd(i),12,4)  
	      << e(gbd(i),12,4) 
	      << e(gcd(i),12,4) << "\n";
  }

  cout << "\n*****  Testing Gradients for FDNLF1 class *****\n";
  grad = 0.0;
  FDNLF1 fdnlp(n,rosen0,init_rosen);
  fdnlp.initFcn();
  fdnlp.eval();
  for(i = 1; i <=n; i++) fdnlp.setFcnAccrcy(i,accrcy);
  fdnlp.printState("FDNLF1");

  grad = fdnlp.evalG();

//  fdnlp.setDerivOption(ForwardDiff);
  gfd = fdnlp.evalG();

  fdnlp.setDerivOption(BackwardDiff);
  gbd = fdnlp.evalG();

  fdnlp.setDerivOption(CentralDiff);
  gcd = fdnlp.evalG();

  for (i=1; i<=n; i++) {
    cout << i << e(grad(i),12,4)  
              << e(gfd(i),12,4)  
	      << e(gbd(i),12,4) 
	      << e(gcd(i),12,4) << "\n";
  }

  
  cout <<"\n*****  Testing Hessians  *****\n";

  cout <<"Finite-Difference Hessian for NLF1 using evalH()\n";
  H = nlp.evalH();
  Print (H);
  cout <<"Finite-Difference Hessian for NLF1 using FDHessian\n";
  H = 0.0;
  H = nlp.FDHessian(sx);
  Print (H);
  cout <<"Finite-Difference Hessian for FDNLF1 using evalH()\n";
  H = 0.0;
  H = fdnlp.evalH();
  Print (H);
  cout <<"Finite-Difference Hessian for FDNLF1 using FDHessian\n";
  H = 0.0;
  H = fdnlp.FDHessian(sx);
  Print (H);
  
}
