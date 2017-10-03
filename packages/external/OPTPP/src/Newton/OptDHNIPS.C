//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last Modified 01/29/2004 
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

#include "OptDHNIPS.h"
#include <float.h>
#include "OptppArray.h"
#include "cblas.h"
#include "ioformat.h"
#include "Teuchos_LAPACK.hpp"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;



namespace OPTPP {

void OptDHNIPS::reset()
{
   NLP2* nlp2 = nlprob2();
   int   n    = nlp2->getDim();
   nlp2->reset();
   OptimizeClass::defaultReset(n);
   indices = 0;
   HCk_ = 0;
}

void OptDHNIPS::nonLinearConstraintIndices(const SerialDenseVector<int,double>& types)
{ 
       for(int i = 0; i < types.numRows(); i++){
          if(types(i) == NLeqn || types(i) == NLineq)
	        indices.append(i);
       }
       return;
}

void OptDHNIPS::initHessian()
{ 
  NLP2* nlp2 = nlprob2();
  hessl      = nlp2->getHess();

  if(nlp2->hasConstraints()){
     CompoundConstraint* constraints = nlp2->getConstraints();
     int nlncons = constraints->getNumOfNLCons();

     if(nlncons){

       SerialSymDenseMatrix<int,double> Htmp(nlp2->getDim());
       Htmp  = 0.0; 

       SerialDenseVector<int,double> type(constraints->getConstraintType().length());
       type = constraints->getConstraintType();
       /* 
	* Map the Lagrange multipliers for the nonlinear constraints
        * to the corresponding constraints.
        */ 
       nonLinearConstraintIndices(type);


       /* 
        * Set initial constraint Hessian approximations to zero,
        * which guarantees a positive definite inital Hessian.
        * See M. Goldsmith (1999), "Sequential Quadratic Programming
        * Methods Based on Indefinite Hessian Approximations", p.66 
        */
       for(int i = 0; i < nlncons; i++)
	   HCk_.append(Htmp);
     }
  }
  return;
}

SerialSymDenseMatrix<int,double> OptDHNIPS::updateH(SerialSymDenseMatrix<int,double>&Hk, int k)
{ 
  if(k == 0){
      initHessian();
      Hk = hessl;
      return Hk;
  }

  double mcheps  = DBL_EPSILON;
  double sqrteps = sqrt(mcheps);

  NLP2* nlp2   = nlprob2();
  
  // Compute the Hessian of the objective function at the current point
  hessl       = nlp2->evalH();

  /*
   * Compute quasi-Newton approximation to the Hessian of each constraint 
   * NOTE: This is a very naive implementation which increases storage
   * costs by (nlncons + 1)(ndim*ndim).  Alternative methods include LBFGS, 
   * LSR1 or sparse methods.
   *
   * Compact representation of LBFGS, Nocedal and Wright p. 230
   * Bk = B - [BSk Yk]M^{-1}[BSk Yk]^T,
   * where B = B0, Sk = [ s0,....,sk-1], Yk = [y0,...,yk-1],
   * M = [ Sk^TBSk   Lk ; Lk^T   -Dk],
   * (Lk)ij = s_{i-1}^Ty_{j-1} for i > j,  Dk = diag{s^Ty} j = 0:k-1
   * NOTE: Requires nlncons factorizations of k x k matrix
   * and storage of 2 m x m matrices. Modifications needed for k > m.
   *
   * Compact representation of LSR1, Nocedal and Wright p. 232
   * Bk = B + (Yk - BSk)(Dk + Lk + Lk^T -Sk^TBSk)^{-1}(Yk - BSk)^T,
   * NOTE: Requires (nlncons) factorizations of k x k matrix
   * and storage of 2 ndim x k matrices.
   */
  if(nlp2->hasConstraints()){
      
     CompoundConstraint* constraints = nlp2->getConstraints();

     bool evaluate;
     int j, ndim, nlncons; 
     double rts, yts, rmnrm, rnrm, snrm, ynrm, maxres, restol, gamma;

     ndim    = nlp2->getDim();
     nlncons = constraints->getNumOfNLCons();

     SerialDenseVector<int,double> xc, yk, sk, res, Bsk, multipliers;
     SerialDenseMatrix<int,double> cg, cgprev;
     SerialSymDenseMatrix<int,double> Htmp(ndim);

     //multipliers = y & z;
     multipliers.resize(y.length()+z.length());
     for(int i=0;i<multipliers.length();i++)
       {if(i<y.length())
	   {multipliers(i) = y(i);}
	 else {multipliers(i) = z(i-y.length());}
       }

     gamma       = 1.0e8;

     // Compute the current step 
     xc.resize(nlp2->getXc().length());
     xc     = nlp2->getXc();
     sk.resize(xc.length());
     sk     = xc;
     sk  -= xprev; 
     cg.reshape(getConstraintGradient().numRows(),getConstraintGradient().numCols());
     cg     = getConstraintGradient();
     cgprev.reshape(constraintGradientPrev.numRows(),constraintGradientPrev.numCols());
     cgprev = constraintGradientPrev;

     for(j = 0; j < nlncons; j++){
        
       //yk   = cg.Column(indices[j-1]) - cgprev.Column(indices[j-1]);
       yk.resize(cg.numRows()); 
       for(int i=0;i<cg.numRows();i++)
	 {yk(i) = cg(i,indices[j]) - cgprev(i,indices[j]);}

       yts  = sk.dot(yk);
       // snrm = Norm2(sk);
       // ynrm = Norm2(yk);
       snrm = sqrt(sk.dot(sk));
       ynrm = sqrt(yk.dot(yk));
       //res    = yk - HCk_[j-1]*sk;
       res.resize(yk.length());
       res = yk;
       SerialDenseVector<int,double> tmp(HCk_[j].numRows());
       tmp.multiply(Teuchos::LEFT_SIDE,1.0, HCk_[j], sk, 0.0);
       res -= tmp;
       rts    = res.dot(sk);
       // rnrm   = Norm2(res);
       rnrm = sqrt(res.dot(res));
       //rmnrm  = (res*res.t()).Norm1();
       SerialDenseMatrix<int,double> tmp2(res.length(),res.length());
       tmp2.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, res, res, 0.0);
       rmnrm = tmp2.normOne();       

        maxres = res.normInf();
	restol = gamma*fabs(rts/ndim)*(1 + HCk_[j].normOne());

	evaluate = true;
     
        /* 
         * For more information about updating safeguards for
	 * the SR1 Hessian approximation; see A. R. Conn,
	 * N. I. M. Gould, and PH. L. Toint (1991), "Convergence of 
	 * quasi-Newton matrices generated by the symmetric rank one 
	 * update", Mathematical Programming 50:177-195.  and 
	 * M. Goldsmith (1999), "Sequential Quadratic Programming
         * Methods Based on Indefinite Hessian Approximations", pp. 60-62.
	 */
        if(fabs(rts) <= sqrteps*snrm*rnrm || restol < rnrm*rnrm){
	    evaluate = false;
            if (debug_) {
                *optout << "UpdateH: y-Hs = " << e(maxres,12,4) 
	  	        << " is too small\n";
                *optout << "UpdateH: The SR1 update is skipped\n";
            }
        }

	if(evaluate){
          // Perfom SR1 update 
	  // Htmp = HCk_[j-1] + (res * res.t()) / rts;
	  SerialSymDenseMatrix<int,double> tmp3(res.length(), res.length());
	  for (int i=0; i<res.length(); i++) {
	    for (int k=0; k<res.length(); k++)
	      tmp3(i,k) = tmp2(i,k)/rts;
	  }
	  Htmp = HCk_[j];
	  Htmp +=  tmp3;
          HCk_[j] = Htmp;
	  // Htmp.Release();
	}
	Htmp =  HCk_[j];
	Htmp *= multipliers(indices[j]);
	// hessl -= HCk_[j-1].scale(multipliers(indices[j-1]));
	hessl -= Htmp;
     }
  }
  Hk = hessl;
  return Hk; 
}


void OptDHNIPS::printStatus(char *title) // set Message
{
  NLP2* nlp2 = nlprob2();

  *optout << "\n\n=========  " << title << "  ===========\n\n";
  *optout << "Optimization method       = " << method   << "\n";
  *optout << "Dimension of the problem  = " << nlp2->getDim()  << "\n"; *optout << "No. equalities            = " << me      << "\n";
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
Teuchos::LAPACK<int,double> lapack;
 SerialDenseVector<int,double> D(hessl.numRows());
    //SVD(hessl, D);
    int alpha = hessl.numRows();
    int LWORK = max(1,3*alpha-1);
    SerialDenseVector<int,double> WORK(LWORK);
    int INFO;
 lapack.SYEV('N', 'L', alpha, hessl.values(), alpha, D.values(),WORK.values(),3*alpha-1, &INFO);
    *optout << "\nEigenvalues of Hessian";
    FPrint(optout, D);
  }

  nlp2->fPrintState(optout, title);
  fPrintMultipliers(optout, title);

  tol.printTol(optout);
}

} // namespace OPTPP
