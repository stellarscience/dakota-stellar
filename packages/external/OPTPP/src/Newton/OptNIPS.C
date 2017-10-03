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

#include "OptNIPS.h"
#include "cblas.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialSymDenseMatrix;


namespace OPTPP {

void OptNIPS::initHessian()
{ 
  if (debug_) *optout << "OptNIPS::initHessian: \n";
  NLP2* nlp2 = nlprob2();
  Hessian   = nlp2->getHess();
  return;
}

  SerialSymDenseMatrix<int,double> OptNIPS::updateH(SerialSymDenseMatrix<int,double>&Hk, int k)
{ 
  NLP2* nlp2 = nlprob2();
  SerialDenseVector<int,double> xc(nlp2->getXc().length()), yzmultiplier;
  xc = nlp2->getXc();
  // Vertically concatenate vectors y and z
  //yzmultiplier = y & z;
  int alpha = y.length();
  int beta = z.length();
  yzmultiplier.resize(alpha+beta);
  for(int i=0;i<alpha+beta;i++)
    {if(i<alpha){yzmultiplier(i) = y(i);}
      else{yzmultiplier(i) = z(i-alpha);}
    }
  SerialDenseVector<int,double> yzmultiplier2(yzmultiplier.length());
  hessl   = nlp2->evalH(xc);

  if(nlp->hasConstraints()){
     CompoundConstraint* constraintah = nlp2->getConstraints();
     yzmultiplier2 = yzmultiplier;
     yzmultiplier2 *= -1;
     hessl += constraintah->evalHessian(xc, yzmultiplier2);
  }

  Hk    = hessl;
  return Hk;
}

void OptNIPS::printStatus(char *title) // set Message
{
  NLP2* nlp2 = nlprob2();

  *optout << "\n\n=========  " << title << "  ===========\n\n";
  *optout << "Optimization method       = " << method   << "\n";
  *optout << "Dimension of the problem  = " << nlp2->getDim()  << "\n";
  *optout << "No. equalities            = " << me      << "\n";
  *optout << "No. inequalities          = " << mi      << "\n";
  *optout << "Merit Function (0= NormFmu, 1 = Argaez, 2 = Vanderbei) = " 
         << mfcn      << "\n";
  *optout << "Return code               = " << ret_code << " (" 
         << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << nlp2->getFevals() << "\n";
  *optout << "No. gradient evaluations  = " << nlp2->getGevals() << "\n";

  if (debug_) {
    *optout << "\nHessian of the Lagrangian";
    FPrint(optout, hessl);
//  Compute eigenvalues of Hessian
 *optout << "Now computing eigenvalues of Hessian " << "\n";
Teuchos::LAPACK<int,double> lapack;
    int INFO;
    int LWORK = max(1,3*(hessl.numRows())-1);
    SerialDenseVector<int,double> D(hessl.numRows()), WORK(LWORK);;
    lapack.SYEV('N','U',hessl.numRows(),hessl.values(),hessl.numRows(),D.values(),WORK.values(),3*(hessl.numRows())-1,&INFO);
    // SVD(hessl, D);


    *optout << "\nEigenvalues of Hessian";
    FPrint(optout, D);
  }

  nlp2->fPrintState(optout, title);
  fPrintMultipliers(optout, title);

  tol.printTol(optout);
}

} // namespace OPTPP
