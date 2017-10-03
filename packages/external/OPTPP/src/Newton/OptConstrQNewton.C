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
#include <cstring>
#include <cmath>
#else
#include <string.h>
#include <math.h>
#endif

#include "OptConstrQNewton.h"
#include "cblas.h"
#include "ioformat.h"
#include <float.h>


using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;


namespace OPTPP {

//------------------------------------------------------------------------
//
//   Quasi-Newton Method member functions
//   checkDeriv()
//   updateH()
//------------------------------------------------------------------------

// static char* class_name = "OptConstrQNewton";

// Check the analytic gradient with FD gradient
int OptConstrQNewton::checkDeriv()
{ return checkAnalyticFDGrad(); }


//---------------------------------------------------------------------------- 
//
// Update Hessian using a Quasi-Newton update
//
//---------------------------------------------------------------------------- 
SerialSymDenseMatrix<int,double> OptConstrQNewton::updateH(SerialSymDenseMatrix<int,double>& Hk, int k) 
{

  double mcheps = DBL_EPSILON;
  double sqrteps = sqrt(mcheps);

  int i;

  NLP1* nlp = nlprob();
  int nr     = nlp->getDim();
  SerialDenseVector<int,double> grad(nr), xc(nlp->getXc().length());
  xc     = nlp->getXc();
  grad   = nlp->getGrad();

  SerialDenseVector<int,double> D(nr);
// BFGS formula
  
  if (k == 0) { // Initial Hessian is set equal to the Identity Matrix
    Hessian = 0.0;
//    D = sx.AsDiagonal()*sx.AsDiagonal();
//    D = sfx.AsDiagonal();
    double typx, xmax, gnorm;
//    double gamma;
    //gnorm = Norm2(grad);
    gnorm = sqrt(grad.dot(grad));

    // Initialize xmax, typx and D to default values
    xmax   = -1.e30; typx   =  1.0; D      =  1.0;

    for (i=0; i < nr; i++) xmax = max(xmax,fabs(xc(i)));
    if( xmax != 0.0) typx = xmax;
    if( gnorm!= 0.0) D    = gnorm/typx;
    if (debug_) {
      *optout << "UpdateH: gnorm0 = " << gnorm
	<< "typx = " << typx << "\n";
    }
    for (i=0; i < nr; i++) Hessian(i,i) = D(i);
    return Hessian;
  }
  
  SerialDenseVector<int,double> yk(nr), sk(nr), Bsk(nr);
  SerialSymDenseMatrix<int,double> Htmp(nr), Htmp2(nr);
  
  yk = grad;
  yk -= gprev;
  sk = xc;
  sk -= xprev;
  
  if (debug_) {
    Print(yk);
    Print(sk);
  }

  double gts =gprev.dot(sk);
  double yts = yk.dot(sk);
  
  double snorm = sqrt(sk.dot(sk));
  double ynorm = sqrt(yk.dot(yk));
  
  if (debug_) {
    *optout << "UpdateH: gts   = " << gts 
         << "  yts = " << yts << "\n";
    *optout << "UpdateH: snorm = " << snorm 
         << "  ynorm = " << ynorm << "\n";
  }

  if (yts <= sqrteps*snorm*ynorm) {
    if (debug_) {
      *optout << "UpdateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    Hessian = Hk; return Hk;
  }
  
  SerialDenseVector<int,double> res(nr);
  //res = yk - Hk*sk;
  SerialDenseVector<int,double> restmp(Hk.numRows());
  restmp.multiply(Teuchos::LEFT_SIDE,1.0,Hk,sk,0.0);
  res = yk;
  res -= restmp;
  if (res.normInf() <= sqrteps) {
    if (debug_) {
      *optout << "UpdateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    Hessian = Hk; return Hk;
  }
  
  // Bsk = Hk*sk;
  Bsk.multiply(Teuchos::LEFT_SIDE,1.0,Hk,sk,0.0);
  double sBs = sk.dot(Bsk);
  double etol = 1.e-8;

  if (sBs <= etol*snorm*snorm) {
    if (debug_) {
      *optout << "UpdateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    Hk = 0;
    for (i=0; i < nr; i++) Hk(i,i) = sx(i)*sx(i);
    Hessian = Hk; return Hk;
  }
  
// Otherwise update the Hessian approximation
  if (debug_) {
    //    *optout << "\nUpdateH: before update, k = " << k << "\n";
    //    FPrint(optout, Hk);
  }

//   //Htmp = - (Bsk * Bsk.t()) / sBs;
//    for(int i=0; i<nr; i++)
//     for(int j=0; j<=i; j++)
//       {Htmp(i,j) = -Bsk(i)*Bsk(j)/sBs;}

//   //Htmp = Htmp + (yk * yk.t()) / yts;
//   for(int i=0; i<nr; i++)
//     for(int j=0; j<nr; j++)
//       {Htmp2(i,j) = yk(i)*yk(j)/yts;} 
//    Htmp += Htmp2;

//    Htmp = Hk;
//    Htmp += Htmp;
//   Hk = Htmp;

 SerialDenseMatrix<int,double> Htmp3(Htmp.numRows(),Htmp.numCols());
  Htmp3.multiply(Teuchos::NO_TRANS,Teuchos::TRANS,-1/sBs,Bsk,Bsk,0.0);

  for(int i=0;i<Htmp3.numRows();i++)
    for(int j=0; j<=i;j++)
      {Htmp(i,j) = Htmp3(i,j);}
 
  // Htmp = Htmp + (yk * yk.t()) / yts;
  Htmp3.multiply(Teuchos::NO_TRANS,Teuchos::TRANS,1/yts,yk,yk,0.0);
   for(int i=0; i<nr; i++)
    for(int j=0; j<=i; j++)
      {Htmp2(i,j) = Htmp3(i,j);} 

   Htmp += Htmp2;


   Htmp += Hk;


  Hk = Htmp;
  // Htmp.Release(); 
  SerialDenseVector<int,double> Bgk(nr);
  //Bgk = Hk*grad;
  Bgk.multiply(Teuchos::LEFT_SIDE, 1.0, Hk, grad, 0.0);
  double gBg = grad.dot(Bgk);
  double gg  = grad.dot(grad);
  double ckp1= gBg/gg;
  if (debug_) {
    //    *optout << "\nUpdateH: after update, k = " << k << "\n";
    //    FPrint(optout, Hk);
    *optout << "UpdateH: sBs  = " << sBs << "\n";
    *optout << "UpdateH: ckp1 = " << ckp1 << "\n";
  }
  Hessian = Hk;
  return Hk;
}

} // namespace OPTPP
