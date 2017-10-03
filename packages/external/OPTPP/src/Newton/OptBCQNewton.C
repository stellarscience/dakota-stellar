//------------------------------------------------------------------------
// Copyright (C) 1996:
// Scientific computing department
// Sandia National Laboratories
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstring>
#else
#include <string.h>
#endif

#include "OptBCQNewton.h"
#include <float.h>
#include "cblas.h"
#include "ioformat.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

namespace OPTPP {

static const char* class_name = "OptBCQNewton";

//------------------------------------------------------------------------
// BCQNewton functions
// initHessian
// checkConvg
// checkDeriv
// initOpt
// printStatus
// stepTolNorm
// updateH
// computeSearch
// updateConstraints
// reset
//------------------------------------------------------------------------
void OptBCQNewton::initHessian()
{ 
  NLP1* nlp = nlprob();
  int   i,n = nlp->getDim();
  Hessian.reshape(n);
  Hessian = 0.0;
  for (i=0; i<n; i++) Hessian(i,i) = 1.0;
  return;
}

int OptBCQNewton::checkDeriv() // Check the analytic gradient with FD gradient
{ return checkAnalyticFDGrad(); }

  //double OptBCQNewton::stepTolNorm() const
  //{
  //  NLP1* nlp = nlprob();
  // SerialDenseVector<int,double> step(sx.AsDiagonal()*(nlp->getXc() - xprev));
  // return Norm2(step);
  //  SerialDenseVector<int,double> tmp(nlp->getXc().length());
  //  tmp = nlp->getXc();
  //  tmp -= xprev;
  //  SerialDenseVector<int,double> step(tmp.length());
  //  for(int i=0; i<tmp.length(); i++)
  //    {step(i) = tmp(i)*sx(i);}
  //  return sqrt(step.dot(step));
  //}

SerialSymDenseMatrix<int,double> OptBCQNewton::updateH(SerialSymDenseMatrix<int,double>& Hk, int k) 
{
  double mcheps = DBL_EPSILON;
  double sqrteps = sqrt(mcheps);

  NLP1* nlp = nlprob();
  int i, nr = nlp->getDim();
  SerialDenseVector<int,double> grad(nr), xc(nlp->getXc().length());
  xc     = nlp->getXc();
  grad   = nlp->getGrad();

  SerialDenseVector<int,double> D(nr);

  // BFGS formula
  if (k == 0) { 
    Hessian = 0.0;
    double typx, xmax, gnorm;

    // Initialize xmax, typx and D to default values
    xmax   = -1.e30; typx   =  1.0; D      =  1.0;

    gnorm = sqrt(grad.dot(grad));

    for (i=0; i < nr; i++) xmax = max(xmax,xc(i));
    if(xmax != 0.0) typx = xmax;
    if(gnorm!= 0.0) D    = gnorm/typx;
    if (debug_) {
      *optout << "updateH: gnorm0 = " << gnorm
	<< "typx = " << typx << "\n";
    }
    for (i=0; i < nr; i++) Hessian(i,i) = D(i);
    return Hessian;
  }
  
  // update the portion of H corresponding to the free variable list only

  SerialDenseVector<int,double> yk(nr), sk(nr), Bsk(nr);
  SerialDenseMatrix<int,double> Htmp(nr,nr);
  
  yk = grad;
  yk  -= gprev;
  sk = xc;
  sk -= xprev;
  for (i=0; i<nr; i++) {
    if (work_set(i) == true) {
      yk(i) = sk(i) = 0.0;
      //*optout << "fixed variable = " << i << "\n";
    }
  }
  
  double gts = gprev.dot(sk);
  double yts = yk.dot(sk);
  
  double snorm = sqrt(sk.dot(sk));
  double ynorm = sqrt(yk.dot(yk));
  
  if (debug_) {
    *optout << "updateH: gts   = " << gts 
         << "  yts = " << yts << "\n";
    *optout << "updateH: snorm = " << snorm 
         << "  ynorm = " << ynorm << "\n";
  }

  if (yts <= sqrteps*snorm*ynorm) {
    if (debug_) {
      *optout << "updateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "updateH: The BFGS update is skipped\n";
    }
    Hessian = Hk; return Hk;
  }
  
  SerialDenseVector<int,double> res(nr),tmp3(Hk.numRows());
  //res = yk - Hk*sk;
  tmp3.multiply(Teuchos::LEFT_SIDE, 1.0, Hk, sk, 0.0);
  res = yk;
  res -= tmp3;
  for (i=0; i<nr; i++) if (work_set(i) == true) res(i) = 0.0;
  if (res.normInf() <= sqrteps) {
    if (debug_) {
      *optout << "updateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "updateH: The BFGS update is skipped\n";
    }
    Hessian = Hk; return Hk;
  }
  
  //Bsk = Hk*sk;
  Bsk.multiply(Teuchos::LEFT_SIDE, 1.0, Hk, sk, 0.0);
  for (i=0; i<nr; i++) if (work_set(i) == true) Bsk(i) = 0.0;
  double sBs = sk.dot(Bsk);
  double etol = 1.e-8;

  if (sBs <= etol*snorm*snorm) {
    if (debug_) {
      *optout << "updateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "updateH: The BFGS update is skipped\n";
    }
    //D = sx.AsDiagonal()*sx.AsDiagonal();
    Hk = 0;
    for (i=0; i < nr; i++) Hk(i,i) = sx(i)*sx(i);
    Hessian = Hk; return Hk;
  }
  
  //Htmp = - (Bsk * Bsk.t()) / sBs;
 Htmp.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, -1/sBs, Bsk, Bsk, 0.0);

 //Htmp = Htmp + (yk * yk.t()) / yts;
 // Htmp = Hk + Htmp;
 // Hk << Htmp;
 // Htmp.Release();
 SerialDenseMatrix<int,double> tmp4(yk.length(),yk.length());
  for(int i=0; i<yk.length();i++)
    for(int j=0; j<yk.length();j++)
      {tmp4(i,j) = yk(i)*yk(j)/yts;}
 Htmp += tmp4;
 for(int i=0; i<Htmp.numRows(); i++)
   for(int j=0; j<=i; j++)
     {Hk(i,j) = Hk(i,j) + Htmp(i,j);}
 
  SerialDenseVector<int,double> Bgk(nr), ggrad(nr);
  //Bgk = Hk*grad;
  Bgk.multiply(Teuchos::LEFT_SIDE, 1.0, Hk, grad, 0.0);
  ggrad = grad;
  for (i=0; i<nr; i++) if (work_set(i) == true) ggrad(i) = 0.0;
  double gBg = ggrad.dot(Bgk);
  double gg  = ggrad.dot(ggrad);
  double ckp1= gBg/gg;
  if (debug_) {
    *optout << "\nupdateH: after update, k = " << k << "\n";
    *optout << "updateH: sBs  = " << sBs << "\n";
    *optout << "updateH: ckp1 = " << ckp1 << "\n";
  }
  Hessian = Hk;
  return Hk;
}

} // namespace OPTPP
