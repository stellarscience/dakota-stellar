//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last Modified December 2000
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <typeinfo>
#ifdef HAVE_STD
#include <cstring>
#include <ctime>
#else
#include <string.h>
#include <time.h>
#endif

#include "OptQNIPS.h"
#include <float.h>
#include "cblas.h"
#include "ioformat.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

namespace OPTPP {

int OptQNIPS::checkDeriv() // check the analytic gradient with FD gradient
{return GOOD;}

SerialSymDenseMatrix<int,double> OptQNIPS::updateH(SerialSymDenseMatrix<int,double>&Hk, int k)
{ 
  double mcheps  = DBL_EPSILON;
  double sqrteps = sqrt(mcheps);

  NLP1* nlp1 = nlprob();
  int ndim  = nlp1->getDim();
  double yts, snrm, ynrm, sBs, maxres;
  SerialDenseVector<int,double> xc, yk, sk, res, Bsk;
  SerialDenseVector<int,double> gradl_curr, gradl_prev, yzmultiplier;
  SerialDenseMatrix<int,double> Htmp(ndim,ndim);

  if(k == 0){
     initHessian();
     Hk = hessl;
     return Hk;
  }

  // Compute change in x and Lagrangian gradient 
  xc.resize(nlp1->getXc().length());
  xc = nlp1->getXc();
  sk.resize(xc.length());
  sk = xc;
  sk -=  xprev; 
  // yzmultiplier = y & z;

  yzmultiplier.resize(y.length()+z.length());  
  for(int i=0;i<y.length()+z.length();i++)
    {if(i<y.length())
	{yzmultiplier(i) = y(i);}
      else{yzmultiplier(i)=z(i-y.length());}
    }
  gradl_curr.resize(getGradL().length());
  gradl_curr   = getGradL();

  if( nlp->hasConstraints() )
    {// gradl_prev = gprev - constraintGradientPrev*yzmultiplier;
      gradl_prev.resize(gprev.length());
      Teuchos::SerialDenseVector<int,double> gtmp(gprev.length());
      gradl_prev = gprev;
      gtmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, constraintGradientPrev, yzmultiplier, 0.0);
      gradl_prev -= gtmp;

    }
  else
      gradl_prev = gprev;
  yk.resize(gradl_curr.length());
  yk = gradl_curr;
  yk -= gradl_prev; 

  // If yts too small, skip bfgs
  yts  = sk.dot(yk);
  //snrm = Norm2(sk);
  //ynrm = Norm2(yk);
  snrm = sqrt(sk.dot(sk));
  ynrm = sqrt(yk.dot(yk));

  if(yts <= sqrteps*snrm*ynrm){
    if (debug_) {
      *optout << "UpdateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    hessl = Hk;
    return Hk;
  }

  //res = yk - Hk*sk;
  res.resize(yk.length());
  res = yk;
  SerialDenseVector<int,double> tmp(Hk.numRows());
  tmp.multiply(Teuchos::LEFT_SIDE, 1.0, Hk, sk, 0.0);
  res -= tmp;
  maxres = res.normInf() ;
     
  if(maxres <= sqrteps){
    if (debug_) {
      *optout << "UpdateH: y-Hs = " << e(maxres,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    hessl = Hk;
    return Hk;
  }

  // If s'Hs too small, skip bfgs
  //Bsk = Hk*sk;
  Bsk.resize(Hk.numRows());
  Bsk.multiply(Teuchos::LEFT_SIDE, 1.0, Hk, sk, 0.0);
    sBs  = sk.dot(Bsk);

  if(sBs <= sqrteps*snrm*snrm){
    if (debug_) {
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    hessl = Hk;
    return Hk;
  }

  // Perfom BFGS update 
  //Htmp  = -(Bsk * Bsk.t()) / sBs;
  Htmp.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, -1/sBs, Bsk, Bsk, 0.0);
  //Htmp +=  (yk * yk.t()) / yts;
  SerialDenseMatrix<int,double> tmp2(yk.length(),yk.length());  
for(int i=0; i<yk.length();i++)
    for(int j=0; j<yk.length();j++)
      {tmp2(i,j) = yk(i)*yk(j)/yts;}
 Htmp += tmp2;
  
 //Htmp  = Hk + Htmp;
 for(int i=0; i<Htmp.numRows(); i++)
   for(int j=0; j<=i; j++)
     {Hk(i,j) = Hk(i,j) + Htmp(i,j);}

 //Hk << Htmp;
 // Htmp.Release();
  if (debug_) {
    *optout << "\nUpdateH: after update, k = " << k << "\n";
    *optout << "UpdateH: sBs  = " << sBs << "\n";
  }
  hessl = Hk;
  return Hk;
}

} // namespace OPTPP
