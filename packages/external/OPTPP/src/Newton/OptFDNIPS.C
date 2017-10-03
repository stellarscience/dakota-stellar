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

#include "OptFDNIPS.h"
#include "cblas.h"
#include "ioformat.h"
#include <float.h>

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;


namespace OPTPP {

int OptFDNIPS::checkDeriv() // check the analytic gradient with FD gradient

{return GOOD;}

SerialSymDenseMatrix<int,double> OptFDNIPS::updateH(SerialSymDenseMatrix<int,double>&Hk, int k)
{ 
  NLP1* nlp1 = nlprob();
  int ndim  = nlp1->getDim();
  SerialDenseVector<int,double> xc, xtmp, gradtmp, yzmultiplier,grad;
  SerialDenseMatrix<int,double> Htmp(ndim,ndim);
  double hi, hieps;

  double mcheps = DBL_EPSILON;
  SerialDenseVector<int,double> fcn_accrcy(nlp1->getFcnAccrcy().length());
  fcn_accrcy = nlp1->getFcnAccrcy();

  Htmp   = 0;
  xc.resize(nlp1->getXc().length());
  xc     = nlp1->getXc();
  // yzmultiplier = y & z;
  int alpha = y.length();
  int beta = z.length();
  yzmultiplier.resize(alpha+beta);

  for(int i=0;i<alpha+beta;i++)
    {if(i<alpha)
	{yzmultiplier(i) = y(i);}
      else{yzmultiplier(i) = z(i-alpha);}
    }

  // Get the gradient of the Lagrangian
  grad.resize(getGradL().length()); 
  grad  = getGradL();

  // Build approximation column by column 
  for(int j = 0; j < ndim ; j++){

     hieps = sqrt(max(mcheps,fcn_accrcy(j)));
     hi     = hieps*max(fabs(xc(j)), sx(j));
     hi     = copysign(hi, xc(j));
     xtmp.resize(xc.length());
     xtmp   = xc;
     xtmp(j)= xc(j) + hi;

     // Evaluate the gradient of the Lagrangian at the new point
     gradtmp.resize(grad.length());
     gradtmp  = nlp1->evalLagrangianGradient(xtmp,yzmultiplier,constrType); 
 
     /* Calculate jth column of Hessian using a first-order 
        finite-difference approximation  */
     //Htmp.Column(j) << (gradtmp - grad)/hi;
     for(int i=0;i<ndim;i++) {
       Htmp(i,j) =(gradtmp(i)- grad(i))/hi;
     }
  }

  // Symmetrize the Hessian Approximation
  for(int i=0;i<ndim;i++)
    for(int j=0;j<=i;j++)
      {		Hk(i,j) =(Htmp(i,j) + Htmp(j,i))/2.0;
      }
  // Hk   << (Htmp + Htmp.t())/2.0;
  // Htmp.Release();
  hessl = Hk;
  return Hk;
}

} // namespace OPTPP
