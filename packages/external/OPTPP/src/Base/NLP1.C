//Norm2 function
//What is mem_grad?
//clarify line 105ish
//transpose

//JWG

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include "NLP1.h"
#include "TOLS.h"
#include "cblas.h"
#include "ioformat.h"
#include <float.h>

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

using namespace std;
using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;


namespace OPTPP {

//----------------------------------------------------------------------------
// Evaluate the Hessian using finite differences
// Assume that analytical gradients are available 
//----------------------------------------------------------------------------

SerialSymDenseMatrix<int,double> NLP1::FDHessian(SerialDenseVector<int,double>& sx) 
{
//  Tracer trace("NLP1::FDHessian");
  double mcheps =DBL_EPSILON;
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();

  int i;
  double hi, hieps;
  double xtmp;


  int nr = getDim();

  SerialDenseVector<int,double> gx(nr), grad(nr), xc(nr);
  SerialDenseMatrix<int,double> Htmp(nr,nr);
  SerialSymDenseMatrix<int,double> H(nr);
  SerialDenseVector<int,double> z(nr);		     
  CompoundConstraint* constraints = getConstraints();
  bool scaleStep = false;
  xc = getXc();
  gx = getGrad();

  for (i=0; i<nr; i++) {
    xc = NLP0::perturbX(i, xc, sx(i), *constraints, 
			fcn_accrcy(i), hi, scaleStep, ForwardDiff);
    grad = evalG(xc);
    z = grad;
    z -= gx;
    z*=(1/hi);
    for(int j=0;j<nr;j++)
      Htmp(j,i) = z(j);
    //Htmp[i] = z;
    //    xc(i) = xtmp;
 }
  for(i=0;i<nr;i++)
    for(int j=0;j<=i;j++)
      H(i,j) = (Htmp(j,i)+Htmp(i,j))/2.0;

  //H << (Htmp.t() + Htmp)/2.0;
 return H;
}

OptppArray<SerialSymDenseMatrix<int,double> > NLP1::CONFDHessian(SerialDenseVector<int,double>& sx) 
{
//  Tracer trace("NLP1::FDHessian");
  double mcheps =DBL_EPSILON;
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();

  int i, counter;
  double hi, hieps;
  double xtmp;

  int nr = getDim();

  SerialDenseVector<int,double> xc(nr);
  SerialDenseMatrix<int,double> grad(nr, ncnln), gx(nr, ncnln), Htmp(nr,nr);
  SerialSymDenseMatrix<int,double> H(nr);

  OptppArray<SerialSymDenseMatrix<int,double> > Hessian(ncnln);
	
    xc = getXc();
   gx = evalCG(xc);

    for (counter=0; counter< ncnln; counter++) {

      for (i=0; i<nr; i++) {

       hieps = sqrt(max(mcheps,fcn_accrcy(i) ));
        hi = hieps*max(fabs(xc(i)),sx(i));
       hi = copysign(hi,xc(i));
        xtmp = xc(i);
        xc(i) = xtmp + hi;
        grad = evalCG(xc);
        //Htmp.Column(i) << (grad - gx) / hi;
        for(int j=0;j<nr;j++)
  	{	
  	  Htmp(j,i) = grad(j,i);
  	  Htmp(j,i) -=gx(j,i);
  	  Htmp(j,i) *= (1/hi);
  	}
        xc(i) = xtmp;
     }

   for(i=0;i<nr;i++)
      for(int j=0;j<=i;j++)
        H(i,j) = (Htmp(j,i)+Htmp(i,j))/2.0;   
   // H << (Htmp.t() + Htmp)/2.0;

     Hessian[counter] = H;

    }
  return Hessian;
 }

//-------------------------------------------------------------------------
// Output Routines
//-------------------------------------------------------------------------

void NLP1::printState(const char * s) 
{ // Print out current state: x current, gradient and Function value
  cout << "\n\n=========  " << s << "  ===========\n\n";
  cout << "\n    i\t    xc \t\t grad  \t\t fcn_accrcy \n";
  for (int i=0; i<dim; i++) 
    cout << d(i,6) << e(mem_xc(i),12,4)<< "\t" << e(mem_grad(i),12,4) << "\t"
         << e(mem_fcn_accrcy(i),12,4) << "\n";
  cout <<"Function Value     = " << e(fvalue,12,4) << "\n";
  // double gnorm = Norm2(mem_grad);
  double gnorm = sqrt(mem_grad.dot(mem_grad));
  cout <<"Norm of gradient   = " << e(gnorm,12,4) << "\n";
  cout <<"\n\n==============================================\n\n";
}

void NLP1::fPrintState(ostream *nlpout, const char * s) 
{ // Print out current state: x current, gradient and Function value
  (*nlpout) << "\n\n=========  " << s << "  ===========\n\n";
  (*nlpout) << "\n    i\t    xc \t\t grad  \t\t fcn_accrcy \n";
  for (int i=0; i<dim; i++) 
    (*nlpout) << d(i,6) << e(mem_xc(i),12,4)<< "\t" 
              << e(mem_grad(i),12,4) << "\t"
              << e(mem_fcn_accrcy(i),12,4) << "\n";
  (*nlpout) <<"Function Value     = " << e(fvalue,12,4) << "\n";
  //double gnorm = Norm2(mem_grad);
  double gnorm = sqrt(mem_grad.dot(mem_grad));
  (*nlpout) <<"Norm of gradient   = " << e(gnorm,12,4) << "\n";
//  (*nlpout) <<"Function Accuracy  = " << e(mem_fcn_accrcy,12,4) << "\n";
  (*nlpout) <<"\n\n==============================================\n\n";
}

} // namespace OPTPP
