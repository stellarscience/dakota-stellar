//------------------------------------------------------------------------
// Copyright (C) 1996:
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstring>
#else
#include <string.h>
#endif

#include "OptBCNewton.h"
#include "cblas.h"
#include "ioformat.h"
#include <float.h>
#include "Teuchos_LAPACK.hpp"



using Teuchos::SerialDenseVector;
using Teuchos::SerialSymDenseMatrix;
using Teuchos::SerialDenseMatrix;

namespace OPTPP {

static const char* class_name = "OptBCNewton";

//------------------------------------------------------------------------
// BCNewton functions
// initHessian
// checkDeriv
// checkConvg
// initOpt
// printStatus
// stepTolNorm
// updateH
// computeSearch
// updateConstraints
// reset
//------------------------------------------------------------------------

int OptBCNewton::checkDeriv() // Check the analytic gradient with FD gradient
{
  NLP2* nlp = nlprob2();
  SerialSymDenseMatrix<int,double> Hess(dim), FDHess(dim), ErrH(dim);

  int retcode   = checkAnalyticFDGrad();
  double mcheps   = DBL_EPSILON;
  double third   = 0.3333333;
  double gnorm = nlp->getGrad().normInf();
  double eta   = pow(mcheps,third)*max(1.0,gnorm);
  FDHess       = nlp->FDHessian(sx); 
  Hess         = nlp->getHess();
  ErrH         = Hess;
  ErrH  -=  FDHess;
  double maxerr = ErrH.normInf();
  if(debug_){
     *optout <<"\nCheck_Deriv: Checking Hessian versus finite-differences\n";
     *optout << "maxerror = " << e(maxerr, 12,4) 
            << "tolerance =  " << e(eta, 12,4) << "\n";
  }
  if (maxerr > eta) retcode = BAD;
  return retcode;
}

void OptBCNewton::initHessian()
{ 
  if (debug_) *optout << class_name << "::initHessian: \n";
  NLP2* nlp = nlprob2();
  Hessian = nlp->getHess();
  return;
}


void OptBCNewton::printStatus(char *s) // set Message
{
  NLP2* nlp = nlprob2();
  *optout << "\n\n=========  " << s << "  ===========\n\n";
  *optout << "Optimization method       = " << method << "\n";
  *optout << "Dimension of the problem  = " << nlp->getDim()   << "\n";
  *optout << "No. of bound constraints  = " << nlp->getDim()   << "\n";
  *optout << "Return code               = " << ret_code << " ("
       << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << nlp->getFevals() << "\n";
  *optout << "No. gradient evaluations  = " << nlp->getGevals() << "\n";

  if (debug_) {
    *optout << "Hessian \n";
    Print(Hessian);
  }

  tol.printTol(optout);

  nlp->fPrintState(optout, s);
}

real OptBCNewton::stepTolNorm() const
{
  NLP2* nlp = nlprob2();
  // SerialDenseVector<int,double> step(sx.AsDiagonal()*(nlp->getXc() - xprev));
  //return Norm2(step);

  SerialDenseVector<int,double> tmp(nlp->getXc().length());
  tmp = nlp->getXc();
  tmp -= xprev;
  for(int i=0; i<tmp.length(); i++)
    {tmp(i) = tmp(i)*sx(i);}
  SerialDenseVector<int,double> step(tmp);
  return sqrt(step.dot(step));
}

SerialSymDenseMatrix<int,double> OptBCNewton::updateH(SerialSymDenseMatrix<int,double>&, int)
{
  return nlprob()->evalH();
}

} // namespace OPTPP
