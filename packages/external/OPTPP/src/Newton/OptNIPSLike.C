//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last Modified March 2003
//
// The OptNIPS algorithm is a C++ implementation of NIPSM, a nonlinear
// interior-point code developed under MATLAB by Amr El-Bakry at Rice
// University and NIPSF, a Fortran implementation of the same code written
// by Frederik Saaf.  Additional features include the merit functions 
// proposed by Miguel Argaez and Richard Tapia in "Global Convergence of a
// Primal-Dual Newton Interior-Point Method for Nonlinear Programming
// Using a Modified Augmented Lagrange Function" as well as 
// Robert Vanderbei and David Shanno in "An Interior-Point Algorithm For
// Nonconvex Nonlinear Programming".
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstring>
#include <ctime>
#else
#include <string.h>
#include <time.h>
#endif

#include <string>
#include <typeinfo>
#include <float.h>

using namespace std;

#include "OptNIPSLike.h"
#include "cblas.h"
#include "ioformat.h"
#include "Teuchos_LAPACK.hpp"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;


namespace OPTPP {

static const char* class_name = {"OptNIPSLike"};
//------------------------------------------------------------------------
//
// Now all the other functions that can be generalized
// to all NIPS cases
//
// setMeritFcn
// computeSearch2
// checkConvg
// checkDeriv
// computeStep
// initOpt
// initHessian
// initMultipliers
// optimize
// readOptInput
// printStatus
// reset
//------------------------------------------------------------------------

void OptNIPSLike::setMeritFcn (MeritFcn option)
{
    mfcn = option ;
    if(mfcn == ArgaezTapia){
      setCenteringParameter(0.2e0); setStepLengthToBdry(0.99995);
    } 
    else if(mfcn == NormFmu){
      setCenteringParameter(0.2e0); setStepLengthToBdry(0.8);
    }
    else if(mfcn == VanShanno){
      setCenteringParameter(0.1e0); setStepLengthToBdry(0.95);
    }
}


SerialDenseVector<int,double> OptNIPSLike::computeSearch2(SerialDenseMatrix<int,double>& J, SerialDenseVector<int,double>& rhs)
{  
   SerialDenseVector<int,double> result(rhs.length()); 
//   //result = J.i()*rhs;
//   Teuchos::SerialDenseSolver<int,double> My_Solver;
//   My_Solver.setMatrix(Teuchos::rcp(&J,false));
//   My_Solver.equilibrateMatrix();
//   My_Solver.setVectors(Teuchos::rcp(&result,false),Teuchos::rcp(&rhs,false));
//   My_Solver.equilibrateRHS();
//   My_Solver.solve();
//   std::cout<<"Result from solve = "<<result<<std::endl;
//   return result; 

   result = rhs;
   int INFO;
   int alpha = J.numRows();
   int *IPIV = new int[alpha];
   Teuchos::LAPACK<int,double> lapack;
   lapack.GESV(alpha,1,J.values(),alpha,IPIV,result.values(),result.length(),&INFO);

   if (IPIV != NULL)
     delete[] IPIV;

   return result;
}

int OptNIPSLike::checkConvg() // check convergence
{
  NLP1* nlp = nlprob();
  SerialDenseVector<int,double> xc(nlp->getXc());
  double ftol = tol.getFTol();
  double error, xnorm, rftol; 
  int convg_status;
  SerialDenseVector<int,double> Fzero(getGradL().length()+mi);
  
  Fzero = setupRHS(xc, 0.0);
  error = sqrt(.5*(Fzero.dot(Fzero)));
  //xnorm = Norm2(xc);
  xnorm = sqrt(xc.dot(xc));
  if (me > 0) xnorm  += sqrt(y.dot(y));
  if (mi > 0) xnorm  += (sqrt(z.dot(z)) + sqrt(s.dot(s)));
  rftol = ftol*(1.0+xnorm);
  if(error <= rftol){
     strcpy(mesg,"Algorithm converged - Norm of gradient less than relative gradient tolerance");
     *optout << "L2 norm = " << e(error,12,4) << "  " << "ftol = " 
             << e(ftol,12,4) << "\n";
     convg_status = 2;
  }
  else
     convg_status = 0;

  return convg_status;
}

int OptNIPSLike::checkDeriv() // check the analytic gradient with FD gradient
{return GOOD;}

int OptNIPSLike::computeStep(SerialDenseVector<int,double> sk)
//---------------------------------------------------------------------------- 
// 
// Compute a step along the direction sk a line search
//
//---------------------------------------------------------------------------- 
{
  NLP1* nlp    = nlprob();
  double lstol   = tol.getLSTol();
  double stpmax  = tol.getMaxStep();
  double stpmin  = tol.getMinStep();
  int  MAXITER = tol.getMaxBacktrackIter();

  if (trace) *optout << class_name << ": ComputeStep\n";


//------------------------------------------------------------------------
// Constrained version of backtrack.C part of the original OPT++ package
// All the parameters remain the same.  Changes include replacing f(x)
// with its constrained counterpart.
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
// Modified by P. J. Williams (pwillia@sandia.gov) September 2000
// to backtrack on constrained problems.
//------------------------------------------------------------------------
/****************************************************************************
 * subroutine backtrack
 * 
 *    FIND the next ITERATE BY LINE SEARCH. 

 * PARAMETERS 
 * ---------- 
 * N            --> DIMENSION OF PROBLEM 
 * X(N)         --> OLD ITERATE:   X[K-1] 
 * F            --> FUNCTION VALUE AT OLD ITERATE, F(X) 
 * G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE 
 * P(N)         --> NON-ZERO NEWTON STEP 
 * XPLUS(N)     <--  NEW ITERATE X[K] 
 * FPLUS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLUS) 
 * FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION 
 * MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED 

 * STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE 
 * STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES 
 *                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM 
 * SX(N)        --> DIAGONAL SCALING MATRIX FOR X 
 * IPR          --> DEVICE TO WHICH TO SEND OUTPUT 

 * INTERNAL VARIABLES 
 * ------------------ 
 * SLN              NEWTON LENGTH 
 * RELLENGTH        RELATIVE LENGTH OF NEWTON STEP 
 ***************************************************************************/
  /* Local variables */
  double a, b, disc, one;
  double t1, t2, t3, tlmbda, minlambda;
  double scl, rellength, sln;

  double lambda    =  1.0;
  double plmbda    = -1.0;
  double pcostplus = -1.0;

  int i, iter, n = nlp->getDim();
  SerialDenseVector<int,double> xc(n), xt(n), xyzs, d_lower(n), d_upper(n);
  SerialDenseVector<int,double> yt(me), zt(mi), st(mi), p(n+me+2*mi);
  double  costplus, tmp1, tmp2;
  double  lowerbd  = .1e0; 
  double  upperbd  = .5e0; 
  bool  interpolate = false;
  bool  constraintsExist = nlp->hasConstraints();
  bool  modeOverride = nlp->getModeOverride();

  xc      = nlp->getXc();
  d_lower = FLT_MAX;
  d_upper = FLT_MAX;

  p      = sk;
  //sln    = Norm2(p);
  sln = sqrt(p.dot(p));

  // Compute 2-norm of matrix with x,y,z,s on diagonal
  // and the relative length of the search direction
  rellength = 0.0;
  //  xyzs      = (xc & y) & (z & s);
  int alpha = y.length();
  int beta = z.length();
  int chi = s.length();
  xyzs.resize(alpha+beta+chi+n);
  for(i=0;i<n+alpha+beta+chi;i++)
    {
      if(i<n)
	xyzs(i)=xc(i);
      else if(i<n+alpha)
	xyzs(i) = y(i-n);
      else if(i<n+alpha+beta)
	xyzs(i) = z(i-n-alpha);
      else xyzs(i) = s(i-n-alpha-beta);
    }
  for (i = 0; i < n+me+2*mi; ++i) {
    tmp1      = fabs(p(i));
    tmp2      = max(fabs(xyzs(i)),1.0);
    rellength = max(rellength,tmp1)/tmp2;
  }

  minlambda = stpmin / rellength;

  CompoundConstraint* constraints = nlp->getConstraints();
  constraints->computeDistanceToBounds(xc, d_lower, d_upper);

  for (i=0; i<n; i++) {
    if ((p(i) > 0.0e0) && (d_upper(i)/p(i) < 1.0e0))
      stpmax = min(stpmax, d_upper(i));
    else if ((p(i) < 0.0e0) && (-d_lower(i)/p(i) < 1.0e0))
      stpmax = min(stpmax, d_lower(i));
  }
  
  if (sln >= stpmax) { // STEP LONGER THAN MAXIMUM ALLOWED
    scl = stpmax / sln;
    //p   = p   * scl;
    p.scale(scl);
    sln = stpmax;
 }
  
  if (debug_) {
    *optout <<"\n dirder   = " << dirder_   << "\n";
    *optout <<"sln         = " << sln       << "\n";
    *optout <<"search dir\n";
    for( i=0; i< n+me+2*mi; i++) *optout<< i << "\t" << p(i) << "\n";
  }
  
  /* CHECK IF NEW ITERATE SATISFACTORY.  GENERATE NEW LAMBDA IF NECESSARY. */
  iter = 0;

  while (iter < MAXITER) {
    iter++;
    // Compute candidate point
    //xt = xc + p.Rows(1,n)*lambda;
    for(i=0;i<n;i++)
      {xt(i) = xc(i) + p(i)*lambda;}
    if(me > 0)
      {for(i=0;i<me;i++){
	  yt(i) = y(i)  + p(n+i)*lambda;
	}
      }
    if(mi > 0){
      for(i=0;i<mi;i++)
	{zt(i) = z(i)  + p(n+me+i)*lambda;
	}
      for(i=0;i<mi;i++)
	{st(i) = s(i)  + p(n+me+mi+i)*lambda;
	}
    }

    costplus = merit(1,xt,yt,zt,st);

    if (debug_) {
      *optout
	<< "iter:  "        << iter 
	<< "  cost      = " << e(cost,12,4)
	<< "  costplus  = " << e(costplus,12,4)
        << "  lambda    = " << lambda << "\n";
    }

    
    // Is this an acceptable step ?
    if (costplus <= cost + dirder_ * lstol * lambda) { 
      // Yes !
       step_length = lambda;
       nlp->setX(xt);
       if(mfcn == NormFmu){
	 nlp->eval();
       }
       else {
	 if (!modeOverride)
	   nlp->evalG();
       }
       if(me > 0) setEqualityMultiplier(yt);
       if(mi > 0){
         setSlacks(st);
         setInequalityMultiplier(zt);
       }

       // Update the gradient of the Lagrangian
       SerialDenseVector<int,double> lgtmp(nlp->getGrad().length());
       lgtmp = nlp->getGrad();

       // Update the constraints 
       if( constraintsExist ){
	 if (modeOverride)
	   nlp->getConstraints()->evalCFGH(xt);

         SerialDenseVector<int,double> constraint_value(nlp->getConstraints()->evalResidual(xt));
	 //         SerialDenseVector<int,double> constraint_value(nlp->getConstraints()->evalResidual(xt).length());
	 //	 constraint_value = nlp->getConstraints()->evalResidual(xt);
         setConstraintResidual(constraint_value);
 	 SerialDenseVector<int,double> nl_values(nlp->getConstraints()->getNLConstraintValue().length());
	 nl_values = nlp->getConstraints()->getNLConstraintValue();
         nlp->setConstraintValue(nl_values);

	 int delta = yt.length();
	 int gamma = zt.length();
         SerialDenseVector<int,double> yzmultiplier(delta+gamma);
	 
	 for(i=0;i<delta+gamma;i++)
	   {if(i<delta)
	       {yzmultiplier(i) = yt(i);}
	     else{yzmultiplier(i)=zt(i-delta);}
	   }
	 // yzmultiplier = yt & zt;
         SerialDenseMatrix<int,double> cgradient(nlp->getConstraints()->evalGradient(xt));
	 //         SerialDenseMatrix<int,double> cgradient(nlp->getConstraints()->evalGradient(xt).numRows(),nlp->getConstraints()->evalGradient(xt).numCols());
	 //	 cgradient = nlp->getConstraints()->evalGradient(xt);
	 setConstraintGradient(cgradient);
	 // lgtmp -= cgradient*yzmultiplier;
	 SerialDenseVector<int,double> lgtemporary(cgradient.numRows());
	 lgtemporary.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,cgradient,yzmultiplier,0.0);
	 lgtmp -= lgtemporary;
       }
       setGradL(lgtmp);
       setCost(costplus);
       fcn_evals   = nlp->getFevals();
       grad_evals  = nlp->getGevals();
       backtracks += iter-1;
       if (iter == 1) return(Newton_Step); // Newton Step
       else return(Backtrack_Step); // Backtrack Step
    }

    // Insufficient decrease. Try to backtrack
    if(lambda < minlambda) { // Step size smaller than min allowed
       step_length = lambda;
       nlp->setX(xprev);
       nlp->setF(fprev);
       nlp->setGrad(gprev);
       fcn_evals   = nlp->getFevals();
       grad_evals  = nlp->getGevals();
       backtracks += iter-1;
       return (-1); // Error
    }  
    else if (interpolate) { // use quadratic/cubic interpolation 

       if (iter == 1) {/* FIRST BACKTRACK: QUADRATIC FIT */
	  if (debug_) *optout <<"First Backtrack\n";
	  tlmbda = -dirder_ / ((costplus - cost - dirder_) * 2.);
       }
       else {/* ALL SUBSEQUENT BACKTRACKS: CUBIC FIT */
   	  if (debug_) *optout << "More Backtrack\n";
	  t1 = costplus  - cost - lambda * dirder_;
	  t2 = pcostplus - cost - plmbda * dirder_;
	  t3 = 1. / (lambda - plmbda);
	  a  = t3 * (t1 / (lambda * lambda) - t2 / (plmbda * plmbda));
	  b  = t3 * (t2 * lambda / (plmbda * plmbda) - 
		   t1 * plmbda / (lambda * lambda));
	  disc = b * b - a * 3. * dirder_;
	  if (disc > b * b) {
	    // ONLY ONE POSITIVE CRITICAL POINT, MUST BE MINIMUM
	    one = copysign(1.0, a);
	    tlmbda = (-b + one * sqrt(disc)) / (a * 3.);
	  }
	  else {
	    // BOTH CRITICAL POINTS POSITIVE, FIRST IS MINIMUM 
	    one = copysign(1.0, a);
	    tlmbda = (-(double)b - one * sqrt(disc)) / (a * 3.);
	  }
       }

       if (tlmbda > lambda * upperbd) 
 	  tlmbda = lambda * upperbd;
       else if (tlmbda < lowerbd * lambda) 
 	  tlmbda = lambda * lowerbd;

       plmbda    = lambda;
       lambda    = tlmbda;
       pcostplus = costplus;
   }
   else{
      lambda = rho_ *lambda;
   }
  } 
// Too many iterations, reset to last acceptable iterate and corresponding f
  nlp->setX(xprev);
  nlp->setF(fprev);
  nlp->setGrad(gprev);
  backtracks+= iter;
  return (-1); 
}

void OptNIPSLike::initOpt()
{
  time_t t;
  char *c;

// get date and print out header

  t = time(NULL);
  c = asctime(localtime(&t));

  *optout << "**********************************************************\n";
  *optout << "OPT++ version " << OPT_GLOBALS::OPT_VERSION << "\n";
  *optout << "Job run at " << c << "\n";
  copyright();
  *optout << "**********************************************************\n";

  NLP1* nlp = nlprob();
  bool modeOverride = nlp->getModeOverride();

//
// Read in OPT++ input file if it exists
// Be aware that anything in the input file will
// override any variable set so far
//

  nlp->initFcn();
  readOptInput();

  if (debug_)
    nlp->setDebug();

  ret_code = 0;

  bool constraintsExist   = nlp->hasConstraints();

  if(constraintsExist){
    CompoundConstraint* constraints = nlp->getConstraints();
    SerialDenseVector<int,double> xstart(nlp->getXc().length());
    xstart = nlp->getXc();

    double feas_tol = tol.getCTol();
    if (modeOverride)
      constraints->evalCFGH(xstart);
    bool feasible = constraints->amIFeasible(xstart, feas_tol);
    if (!feasible) {
      *optout << "OptNIPSLike WARNING:  Initial guess not feasible.\n"
	      << "NIPS may be unable to make progress." << endl;

      if( getFeasibilityRecoveryFlag() )
        recoverFeasibility(xstart, constraints, feas_tol);

      SerialDenseVector<int,double> nl_values(nlp->getConstraints()->getNLConstraintValue());
      nl_values = nlp->getConstraints()->getNLConstraintValue();
      nlp->setConstraintValue(nl_values);
    }
  }

  if (ret_code == 0) {
    CompoundConstraint* constraints;
    double gradf_norm, half = 0.5;
    int i, ind, initIteration = 0, n = nlp->getDim();

    // Local Vectors
    SerialDenseVector<int,double> fscale, gradtmp, yzmultiplier;

    // Local Matrices
    SerialDenseMatrix<int,double> constraintGrad;
    SerialSymDenseMatrix<int,double> Hk(n);

    /* Reset number of constraints to zero
     * Otherwise, on subsequent calls to optimize
     * me = Cme , where C is the number of calls to optimize
     */
    me = 0;
    mi = 0;

    if(constraintsExist){
      constraints = nlp->getConstraints();
      int nCons   = constraints->getNumOfCons();
      constrType.resize(nCons);
      constraintResidual.resize(nCons);
      constraintGradient.reshape(n, nCons);
      constraintGradientPrev.reshape(n, nCons);
      constrType  = constraints->getConstraintType();

      for(i = 0; i < nCons; i++){
        if(constrType(i) == Leqn || constrType(i) == NLeqn)
	  me++;
        else
	  mi++;
      }
    }

    // Dimension check
    if( n < 1 ){
      if(debug_)
	*optout << "INITOPT:  MUST HAVE  n >= 1" << "\n" << flush;
      setMesg("MUST HAVE: # parameters >= 1");
      ret_code = -5;
      return;
    }
    if( me < 0 || me >= n){
      if(debug_)
	*optout << "INITOPT:  MUST HAVE 1 <= me <= n-1 " << "\n" << flush;
      setMesg("MUST HAVE: 0 <= # equality constraints < # parameters");
      ret_code = -5;
      return;
    }
    if( mi < 0 ){
      if(debug_)
	*optout << "INITOPT:  MUST HAVE  mi >= 0 " << "\n" << flush;
      setMesg("MUST HAVE: # inequality constraints >= 0");
      ret_code = -5;
      return;
    }

    xprev = nlp->getXc();

    if(constraintsExist){
      // Resize  vectors
      y.resize(me);
      s.resize(mi);
      z.resize(mi);

      // Evaluate constraints at the initial point
      if (modeOverride)
	constraints->evalCFGH(xprev);

      SerialDenseVector<int,double> constraintValue(constraints ->evalResidual(xprev));
      //      SerialDenseVector<int,double> constraintValue(constraints ->evalResidual(xprev).length());
      //      constraintValue = constraints ->evalResidual(xprev);
      setConstraintResidual(constraintValue);
      SerialDenseVector<int,double> nl_values(nlp->getConstraints()->getNLConstraintValue().length());
      nl_values = nlp->getConstraints()->getNLConstraintValue();
      nlp->setConstraintValue(nl_values);

      // Evaluate constraint gradients at the initial point
      //Helpful Comment
      constraintGrad.reshape(n, constraints->getNumOfCons());
      constraintGrad         = constraints->evalGradient(xprev);
      setConstraintGradient(constraintGrad);
      constraintGradientPrev = getConstraintGradient();
    }

    // Evaluate Function, gradient and compute initial Hessian
    nlp->eval();

    fprev  = nlp->getF();
    gprev  = nlp->getGrad();
    setFcnScale(fprev);
    fscale.resize(getFcnScale().length());
    fscale = getFcnScale();

    // Initialize Lagrange multipliers and slack variables
    if(constraintsExist){
      // gradf_norm   = Norm2(gprev);
      gradf_norm = sqrt(gprev.dot(gprev));
      if(me > 0) y = gradf_norm;
      for(i = 0; i < mi; i++){
        ind  = me + i;
        s(i) = max(half, constraintResidual(ind)); 
        z(i) = max(gradf_norm, .1); 
      }
      // 09/07/01 Least squares initialization of multipliers
      //if( mi > 0){
      //  z = initMultipliers(gprev, constraintGradientPrev);
      //}
    }

    // Initialize mu 
    if(mi > 0)   updateMu(initIteration);

    // Initialize beta_ to an "appropriately" large number if no inequalities 
    if(mi == 0)  beta_ = 100;

    // Concatenate vectors y and z to prevent passing multiple vectors
    //yzmultiplier =  y & z;
    yzmultiplier.resize(y.length()+z.length());
    for(i=0;i<y.length()+z.length();i++)
      {if(i<y.length())
	  {yzmultiplier(i) = y(i);}
	else{yzmultiplier(i) = z(i-y.length());}
      }

    // Evaluate the gradient of the Lagrangian at the initial point
    if (constraintsExist){
      gradtmp.resize(constraintGradient.numRows());
      gradtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-1.0,constraintGradient,yzmultiplier,0.0);
      gradtmp += gprev;
      
      // gradtmp = gprev - constraintGradient*yzmultiplier;
      setGradL(gradtmp);
    }
    else {
      setGradL(gprev);
    }

    // Print the initial state to the screen
    if (debug_)
      nlp->printState("Initial state");

    // Get Optimization Parameters
    *optout << "\n\t\t NIPS Iteration Summary\n";
    *optout << "\n\t" << method << " Method with Line Search \n ";
    *optout << "\n  Merit Function =  " << mfcn << " \n";
    *optout << "\n  Iter    F(x)       mu          alpha"
            << "      Merit    feval    btracks  Penalty\n\n";

    *optout << d(0,5) << " "  << e(fprev,12,4) << " " << e(mu_,12,4) 
            << " "    << "\n" << flush;

    hessl = updateH(Hk,initIteration);

    if (debug_) {
      nlp->fPrintState(optout, "OptNIPSLike: Initial Guess");
      *optout   << "index,  xc, grad \n";
      for(i=0; i<n; i++)
	*optout << i << e(xprev(i),12,4) << e(gprev(i),12,4) << "\n";
      FPrint(optout, hessl);
    }
  }
}

void OptNIPSLike::initHessian()
{ 
  NLP1* nlp = nlprob();
  int i, ndim  = nlp->getDim();

  if (WarmStart) {
    *optout << "OptNIPSLike::initHessian: Warm Start specified\n";
  }
  else {
    double typx, xmax, gnorm;
    SerialDenseVector<int,double> grad(ndim), xc(ndim);
    SerialDenseVector<int,double> D(ndim);
    xc     = nlp->getXc();
    grad   = nlp->getGrad();
    //gnorm  = Norm2(grad);
    gnorm = sqrt(grad.dot(grad));

    // Initialize xmax, typx and D to default values
    xmax   = -1.e30; typx   = 1.0; D      = 1.0;

    for (i = 0; i < ndim; i++) xmax = max(xmax,fabs(xc(i)));
    if(xmax != 0.0) typx   = xmax;
    if(gnorm!= 0.0) D      = gnorm/typx;
    if (debug_) {
      *optout << "OptNIPSLike::initHessian: gnorm0 = " << gnorm
	      << "  typx = " << typx << "\n";
    }
    hessl  = 0.0;
    for (i=0; i < ndim; i++) hessl(i,i) = D(i);
   }
}

SerialDenseVector<int,double> OptNIPSLike::initMultipliers(const SerialDenseVector<int,double>& df, 
                                          SerialDenseMatrix<int,double>& dcon)

{ 
    int i, gm, gn; 
    //double gnorm = Norm2(df);
    double gnorm = sqrt(df.dot(df));
    SerialDenseVector<int,double> b, sk, ztmp;
    SerialDenseMatrix<int,double> dg(dcon.numRows(), mi);
    double mcheps = DBL_EPSILON;
    
    // Drop columns of the constraint gradient corresponding to equality constraints
    if( me > 0){
        
      for(i=0;i<dcon.numRows();i++)
	for(int j=0;j<mi;j++)
	  {dg(i,j)=dcon(i,me+i);
	  }
    }
	    //dg.Columns(1,mi) = dcon.Columns(me+1,me+mi); 

    else        dg   = dcon;
    
    gm  = dg.numRows();
    gn  = dg.numCols();

    if (gm >= gn){
      SerialDenseMatrix<int,double> L(gn,gn),Htmp(gn,gn);
      SerialSymDenseMatrix<int,double> H(gn);

      // H << dg.t()*dg;
      Htmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, dg, dg, 0.0);
      for(i=0;i<gn;i++)
	for(int j=0;j<=i;j++)
	  {H(i,j) = Htmp(i,j);}
      // b =  dg.t()*df;
      b.resize(dg.numCols());
      b.multiply(Teuchos::TRANS,Teuchos::NO_TRANS, 1.0, dg, dg, 0.0);  

      L    = MCholesky(H);
      // ztmp = (L.t().i()*(L.i()*b));
      ztmp.resize(b.length());
      ztmp = b;
  
  //Using LAPACK to solve two triangular systems
Teuchos::LAPACK<int,double> lapack;
    int INFO;
   lapack.TRTRS('L','N','N',gn,1,L.values(),gn, ztmp.values(),gn, &INFO);
   lapack.TRTRS('L','T','N',gn,1,L.values(),gn, ztmp.values(),gn, &INFO);

   // Teuchos::SerialSpdDenseSolver<int,double> My_Solver;
   // int info = 0;
 
   // My_Solver.setMatrix(Teuchos::rcp(&H,false));
   // My_Solver.setVectors(Teuchos::rcp(&ztmp, false), Teuchos::rcp(&b,false));
   // My_Solver.equilibrateMatrix();
  //My_Solver.equilibrateRHS();
   // info = My_Solver.solve();
   // if(info != 0)
   // {return ztmp;}
    }
    else{
      SerialDenseMatrix<int,double> L(gm,gm),Htmp(gm,gm);
      SerialSymDenseMatrix<int,double> H(gm);

      // H << dg*dg.t();
      Htmp.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, dg, dg, 0.0);
      for(i=0;i<gn;i++)
	for(int j=0;j<=i;j++)
	  {H(i,j) = Htmp(i,j);}
       b =  df;  

      L  = MCholesky(H);
      // sk = (L.t().i()*(L.i()*b));
      // ztmp = dg.t()*sk;
      sk.resize(b.length());
      sk = b;
      return sk;
 
  //Using LAPACK to solve two triangular systems
Teuchos::LAPACK<int,double> lapack;
    int INFO;
   lapack.TRTRS('L','N','N',gm,1,L.values(),gm, sk.values(),gm, &INFO);
   lapack.TRTRS('L','T','N',gm,1,L.values(),gm,sk.values(),gm,&INFO);
   
  
   // Teuchos::SerialSpdDenseSolver<int,double> My_Solver;
   // int info = 0;
 
   // My_Solver.setMatrix(Teuchos::rcp(&H,false));
   // My_Solver.setVectors(Teuchos::rcp(&sk, false), Teuchos::rcp(&gprev,false));
   // My_Solver.equilibrateMatrix();
  //My_Solver.equilibrateRHS();
   // info = My_Solver.solve();
   // if(info != 0)
   // {return sk;}
 ztmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, dg, sk, 0.0);
    }

    // Safeguard against nonpositive values of z
    for(i = 0; i < mi; i++){
       if( ztmp(i) < mcheps ) 
         ztmp(i) = max(.1, gnorm);
    }
    return ztmp;
}
void OptNIPSLike::recoverFeasibility(SerialDenseVector<int,double> xinit, CompoundConstraint* constraints,
                   double ftol)
//---------------------------------------------------------------------------- 
//
// Given an infeasible starting point attempt to find a feasible starting point
// within a specified number of iterations
//
//---------------------------------------------------------------------------- 
{

  bool cc_flag = false;
  double ceps  = .1;
  int count, max_cc_iter;

  count = 1;
  max_cc_iter = getMaxFeasIter();
  
  NLP1* nlp = nlprob();

  while(count <= max_cc_iter && !cc_flag){
    constraints->computeFeasibleInequalities(xinit, ftol);
    constraints->computeFeasibleBounds(xinit, ceps);
    cc_flag = constraints->amIFeasible(xinit,ftol);
    count++;
  }
  nlp->setX(xinit);
  *optout << "\n Restoring feasibility with respect to inequalities ... \n";
  FPrint(optout, xinit);
}

void OptNIPSLike::optimize()
//---------------------------------------------------------------------------- 
//
// Given a nonlinear operator nlp find the minimizer using a
// Newton-like method
//
//---------------------------------------------------------------------------- 
{
  int k;
  int n = dim;
  int step_type;
  int convgd = 0;
  int maxiter, maxfev, fevals;
  double alpha_dmp = 1.0;

// Allocate local vectors 
  SerialDenseVector<int,double> sk(n + me + 2*mi), Fmu(n + me + 2*mi);
  SerialDenseVector<int,double> JtF, yzmultiplier; 

// Allocate local matrices
  SerialDenseMatrix<int,double> Jacobian(n+me+2*mi,n+me+2*mi);
  SerialSymDenseMatrix<int,double> Hk(n);

// Initialize iteration : evaluate Function, Gradient, and Hessian
  initOpt();

  if (ret_code == 0) {
    maxiter = tol.getMaxIter();
    maxfev  = tol.getMaxFeval();

    Hk = hessl;

    NLP1* nlp  = nlprob();

    // Check for convergence. Need to take into account that this is the
    // zeroth iteration
    convgd = checkConvg();
    if (convgd > 0) {
      iter_taken = 0;
      ret_code   = convgd;
      setReturnCode(ret_code);
      return;
    }
  
    for (k = 1; k <= maxiter; k++) {
      iter_taken   = k;
      
      //yzmultiplier = y & z;
      int zeta = y.length();
      int eta = z.length();
      yzmultiplier.resize(zeta+eta);
      for(int i=0;i<zeta+eta; i++)
	{if(i<zeta)
	    {yzmultiplier(i) = y(i);}
	  else {yzmultiplier(i) = z(i-zeta);}
	}
    


      // Select new perturbation mu
     updateMu(k);

      // Construct right-hand side (-PKKT) using new mu
    
      Fmu      = setupRHS(xprev, mu_);

      // Construct Jacobian matrix for the Newton system
      Jacobian = setupMatrix(xprev);

      // Compute the derivative of the cost fcn ||F|| 
      //JtF      = Jacobian.t()*Fmu;
      JtF.resize(Jacobian.numCols());
      JtF.multiply(Teuchos::TRANS,Teuchos::NO_TRANS, 1.0, Jacobian, Fmu, 0.0);
      SerialDenseVector<int,double> Fmu2(n+me+2*mi);
      Fmu2 = Fmu;

      // Solve for the Newton search direction
      try{
        Fmu2 *= -1;
	sk       = computeSearch2(Jacobian, Fmu2 );
      }
      catch(...){
        cout << "\n Algorithm terminated - Singular Jacobian \n";
        setMesg("Algorithm terminated - Singular Jacobian");
        //setReturnCode(-5);
        return;
      }

      // Dampen the step to ensure feasibility of the nonnegative iterates 
      if(mi > 0) alpha_dmp = dampenStep(sk);

      // Compute the directional derivative 
      computeDirDeriv(sk,xprev,JtF);

      if(dirder_ >= 0){
	if (debug_) {
	  *optout << "NIPS - Optimize: Warning Directional "
		  << "Derivative >= 0 "<< "\n";
	  *optout << "NIPS - Optimize: dirder = "  << e(dirder_,12,4) << endl;
	}
      }

      // 09/04/01 PJW
      // Placed the merit function call after the
      // computation of the directional derivative and penalty
      // parameter updates for Argaez-Tapia merit function.
      // This is consistent with the description of an nonlinear
      // interior-point algorithm provided in Argaez and Tapia, 1995.
      // In Vanderbei et al, "The current version of the code starts
      // the next iteration from the new point with the increased value
      // of beta.  In future work, we will test the efficiency gained
      // by simply using the new beta with the old step directions".

      // Evaluate the merit function
      setCost (  merit(0,xprev,y,z,s) ); 
      
      if ((step_type = computeStep(sk)) < 0) {
	*optout << "step_type = " << step_type << "\n";
	setMesg("Algorithm terminated - Number of backtracks exceeds the specified limit");
	ret_code = step_type;
        setReturnCode(ret_code);
	return;
      }

      // Update the Hessian of the Lagrangian
      hessl = updateH(Hk,k);
      Hk    = hessl;

      // Accept this step and update the nonlinear model
      xprev = nlp->getXc();
      fprev = nlp->getF();
      gprev = nlp->getGrad();
      constraintGradientPrev  = getConstraintGradient();
      updateModel(k, n, xprev);

      // Retrieve the number of function evaluations
      fevals     = nlp->getFevals();

      // Print iteration summary
      *optout 
        << d(k,5) << " " << e(fprev,12,4) << " " << e(mu_,12,4) 
	<< e(alpha_dmp*step_length,12,4)  << " " << e(cost,12,4)   
        << " " << d(fevals,4)  << " " << d(backtracks,3) 
        << " " << e(penalty_,10,2) <<  endl;

      // Test for algorithmic convergence
      convgd     = checkConvg();
      if (convgd > 0) {
	ret_code = convgd;
        setReturnCode(ret_code);
	return;
      }

      if (fevals > maxfev) break;
	}

    setMesg("OptNIPSLike:  Maximum number of iterations or fevals");
    ret_code = -4;
    setReturnCode(ret_code);
  }
  setReturnCode(ret_code);
}

//-----------------------------------------------------------------------
//  updateMu
//  merit
//  merit2
//  merit3
//  computeDirDeriv
//  dampenStep
//  setupRHS
//  setupMatrix
//-----------------------------------------------------------------------

void OptNIPSLike::updateMu(int k)
{  
  const double r     = 0.95e0; 

  const double zero  = 0.0e0;
  const double one   = 1.0e0;
  const double two   = 2.0e0;

  const double magic_constant     = 1.0e-2;

  int i;
  double t3, t4;
  double min_sz, squiggle;
  double hnormsq, sum, newmu;
  double dotprod, sigma, gamma = sigmin_; 

  SerialDenseVector<int,double> conresid(me+mi);
  NLP1* nlp = nlprob();
  SerialDenseVector<int,double> xc(nlp->getXc().length());
  xc        = nlp->getXc();


  if(mi > 0){
    dotprod = s.dot(z);
     sigma   = min(sigmin_,sw_*dotprod); 

     if(mfcn == NormFmu)
        setMu(sigma*dotprod/mi);
     else if(mfcn == ArgaezTapia){
        if(k == 0)
           setMu(sigma*dotprod/mi);
        else{
           sum       = zero;
	   hnormsq   = zero;
           conresid  = getConstraintResidual();
           for(i = 0; i < mi; i++){
              sum  += s(i)*z(i) + pow(mu_,2)/(s(i)*z(i));
	      conresid(me + i) -= s(i);
           }
           sum      -= 2*mu_*mi;
	   hnormsq   = conresid.dot(conresid);
           newmu     = hnormsq + sum;
	   if(newmu <= mu_/two)
	      setMu(mu_*magic_constant);
        }
     }
     else if(mfcn == VanShanno){
        min_sz    = s(1)*z(1);
        for(i = 0; i < mi; i++)
            min_sz  = min(min_sz, s(i)*z(i));
        t3        = dotprod/mi; 
        squiggle  = min_sz/t3;
	t4        = min((one-r)*(one-squiggle)/squiggle,two);
        setMu(gamma*(pow(t4,3))*t3);
    }
  }
}


double OptNIPSLike::merit(int flag, const SerialDenseVector<int,double>& xc, 
	            const SerialDenseVector<int,double>& yc, SerialDenseVector<int,double>& zc, 
		    SerialDenseVector<int,double>& sc){
//------------------------------------------------------------------------
//     Choose merit function to call
//
//     Parameters
//     flag 
//          1 NormFmu
//          2 Argaez, Tapia
//          3 Vanderbei, Shanno
//------------------------------------------------------------------------
      double  tcost;
      const double zero = 0.0e0;
      SerialDenseVector<int,double> Fmu;

      if (mfcn == NormFmu){ 
	 if(flag)
	    Fmu = setupRHS(xc,yc,zc,sc,zero);
	 else
	    Fmu = setupRHS(xc,zero);
         tcost  = .5*(Fmu.dot(Fmu));
      }
      else if(mfcn == ArgaezTapia) 
         tcost  = merit2(flag,xc,yc,zc,sc);
      else if(mfcn == VanShanno) 
         tcost  = merit3(flag,xc,zc,sc);
      else
         cout << "Merit: Error in merit flag choice, flag =  " << mfcn 
	      << endl;
      return tcost;
}

double OptNIPSLike::merit2(int flag, const SerialDenseVector<int,double>& xc, 
	                 const SerialDenseVector<int,double>& yc, 
                         SerialDenseVector<int,double>& zc, SerialDenseVector<int,double>& sc)
{
// --------------------------------------------------------------------
// Double-precision function to compute the Argaez-Tapia Merit function
//  
// In the case of simple bounds, the Argaez-Tapia merit fcn has the
// following form.
//    M(x,y,z)    = l(x,y,z) + rho*p(x,z), where
//    l(x,y,z)    = f(x) + h(x)'y - x'z
//    p(x,z)      = h(x)'h(x)/2  + x'z - mu(sum(log(xizi))
//
// Now, we consider the general nlp.
//    M(x,y,z,s,w)= l(x,y,z,s,w) + rho*p(x,z,s), where
//    l(x,y,z,s,w)= f(x) + h(x)'y -(g(x) - s)'w - s'z
//    p(x,z,s)    = h(x)'h(x)/2  + (g(x)-s)'(g(x) -s)/2 +s'z - mu(sum(log(sizi))
//
// From the KKT conditions, we see that w = z which reduces 
//    l(x,y,z)    = f(x) + h(x)'y - g(x)'z 
//
// OUTPUT:
//     merit2
// --------------------------------------------------------------------

// Local variables
  double t1, t2, t3, t4, sum;
  double phi, lagrangian, merit2;
  const double tiny  = 1.0e-30;

  NLP1* nlp             = nlprob();
  bool constraintsExist = nlp->hasConstraints();
  bool modeOverride     = nlp->getModeOverride();
  SerialDenseVector<int,double> conresid(me+mi), yzmultiplier;
  int theta = yc.length();
  int iota = zc.length();
  yzmultiplier.resize(theta+iota);
  for(int i=0;i<theta+iota;i++)
    {if(i<theta){yzmultiplier(i) = yc(i);}
      else{yzmultiplier(i) = zc(i-theta);}
    }
  // yzmultiplier = yc & zc;

// --------------------------------------------------------------------
//     Compute Penalty term
// --------------------------------------------------------------------

  t1 = 0.0e0;
  t2 = 0.0e0;
  t3 = 0.0e0;

  if(flag){
   // If within computeStep, computed values based on trial step 

    if( constraintsExist ){
      if (modeOverride) 
        nlp->getConstraints()->evalCFGH(xc);
      conresid  = nlp->getConstraints()->evalResidual(xc);
    }
    else
      conresid  = 0;

     // Compute Lagrangian
     if (modeOverride) {
       nlp->setX(xc);
       nlp->eval();
       lagrangian = nlp->getF();
     }
     else {
       lagrangian  = nlp->evalF(xc);
       nlp->setF(lagrangian);
     }
     if( constraintsExist ) 
       lagrangian -=  conresid.dot(yzmultiplier);
  }
  else{
     // Use previously computed values 
     lagrangian = nlp->getF();
     if( constraintsExist ) { 
	 conresid    = getConstraintResidual();
	 lagrangian -= conresid.dot(yzmultiplier);
     }
  }

  if(me > 0){
    sum = 0.0;
    for(int i = 0; i < me; i++)
      sum += pow(conresid(i),2);
    t1  = sum/2.0; 
  }
  if(mi > 0){
    sum = 0.0;
    t2  = sc.dot(zc);
    for(int i = 0; i < mi; i++){
      conresid(me+i)-=sc(i);
      t1  += pow(conresid(me+i),2)/2.0;
      t4   = max(tiny,sc(i)*zc(i));
      sum += log(t4); 
    }
    t3  = mu_*sum;
  }

  phi     = t1 + t2 - t3;
  merit2  = lagrangian + penalty_ *phi;
      
  if(debug_){      
     *optout <<  "merit2:" << "\n";
     *optout << " lagrangian    t1    t2    t3     mu   rho phi merit2 "
       << "\n";
     *optout << e(lagrangian,12,4) << " "  << e(t1,12,4)  << " " << e(t2,12,4)
       << " " << e(t3,12,4)   << " "    << e(mu_,12,4)
       << " " << e(penalty_,12,4) << " "    << e(phi,12,4)
       << " " << e(merit2,12,4)        << "\n";
         
     *optout <<  "merit2: s, z" << "\n";
     for(int i = 0; i < mi; i++){
        *optout << i << " " << e(s(i),12,4)  << " " << e(z(i),12,4)
   	      << "\n ";
     *optout << flush;
     }
  }

  return merit2;
}


double OptNIPSLike::merit3(int flag, const SerialDenseVector<int,double>& xc, 
                         SerialDenseVector<int,double>& zc, SerialDenseVector<int,double>& sc)
{
// --------------------------------------------------------------------
// Double-precision function to compute the Vanderbei-Shanno Merit function
//
// Authors: Juan Meza and Patty Hough
//
// INPUTS:
//     f
//     h(me)
//     s(mi)
//     y(me)
//     z(mi)
//     mu
//     rho
//     
// OUTPUTS:
//     merit3
// --------------------------------------------------------------------

// Local variables
  int i;
  double sumlog;
  double t2, t3, tcost;
  const double tiny  = 1.0e-30;
  const double zero  = 0.0e0;
  const double two   = 2.0e0;

  SerialDenseVector<int,double> conresid;
  NLP1* nlp = nlprob();
  bool modeOverride = nlp->getModeOverride();

  if(flag){ 
    if (modeOverride){
      nlp->setX(xc);
      nlp->eval();
      tcost = nlp->getF();
    }
    else {
      tcost = nlp->evalF(xc);
      nlp->setF(tcost);
    }
  }
  else
     tcost  = nlp->getF();
           
  sumlog    = zero;

  if(nlp->hasConstraints()){

    if(flag){
      if (modeOverride)
        nlp->getConstraints()->evalCFGH(xc);
      conresid.resize(nlp->getConstraints()->getNumOfCons());
      conresid  = nlp->getConstraints()->evalResidual(xc);
    }
    else
      conresid.resize(getConstraintResidual().length());
      conresid = getConstraintResidual();

     for(i = 1; i <= mi; i++){
        t2     = max(tiny,sc(i));
        sumlog =  sumlog + log(t2);
        conresid(me+i) -= sc(i);
     }

     t3        = conresid.dot(conresid);
     tcost    += ( (beta_/two)*(t3) );
     tcost    -= (mu_ * sumlog);
  }

  return tcost;
}

void OptNIPSLike::computeDirDeriv(SerialDenseVector<int,double>& sk, const SerialDenseVector<int,double>& xc,
                                  SerialDenseVector<int,double>& derivative)
{
   double sum1, sum2, term1, term2, temp, penaltybd = 2.0;
   double ftol   = tol.getFTol();
   NLP1* nlp   = nlprob();
   int i, ind, n = nlp->getDim();
   bool constraintsExist = nlp->hasConstraints();
   SerialDenseVector<int,double> conresid(me+mi), gradf(n), tgradl(n);

/*  
 *  NormFmu directional derivative
 *  gradM'(dv) =  -2(J'F)'(J.inv()F) 
 */
   if (mfcn == NormFmu)
     dirder_ = derivative.dot(sk);
/*  
 *  Argaez-Tapia Merit Function and directional derivative
 *  gradM'(dv) =  gradl'(dx) + h(x)'(dy) - s'(dz)
 *           -rho*(h(x)'h(x) + (g(x)-s)'(g(x)-s)/2 )
 *           -rho*(s'z - 2mumi +mu*mu*e'(SZ^(-1))e))
 */
   else if (mfcn == ArgaezTapia){ 
      sum1     = 0.0e0;
      sum2     = 0.0e0;
      term1    = 0.0e0;
      term2    = 0.0e0;
      tgradl   = getGradL();

      /* gradl'(dx) */
      for(i = 0; i < n ; i++)
           term1  +=  tgradl(i)*sk(i); 

      if( constraintsExist ){  

         conresid = getConstraintResidual();

         // g(x)'(dz) 
         for(i = 0; i < mi ; i++){
            ind    =  n + me + i;
            term1 -= conresid(me+i)*sk(ind);
         }

         /* The penalty terms
          *  mu^2*e'(SZ^(-1))e 
          */
         for(i = 0; i < mi ; i++){
            conresid(me+i) -= s(i);
            sum2  += pow(mu_,2)/(s(i)*z(i));
         }

         // s'z  - 2mumi 
         sum2     +=  s.dot(z); 
         sum2     -=  2*mu_*mi;

         // h(x)'h(x) + (g(x)- s)'(g(x) - s) 
         sum1      =  conresid.dot(conresid);

         // s'z - 2mumi + mu^2*e'(SZ^(-1))e + h(x)'h(x) + (g(x) - s)'(g(x) - s) 
         term2     = sum1 + sum2;
      }

      dirder_ = term1 - penalty_*(term2);

      if (debug_)
	*optout << "\n Directional directive:  " << dirder_ << flush;
       
      if( fabs(term2) > ftol)
        temp  = term1/fabs(term2);
      else
        temp  = term1; 

      if(temp + penaltybd > penalty_) 
        penalty_ = temp + penaltybd;

      if (debug_)
       *optout << "\n Update Penalty : " << penalty_ << flush;
   }   
/*  
    Vanderbei-Shanno directional derivative
    gradM'(dv) =  gradf'(dx) + -mu*(sum(ds/s)) - beta( (g(x) - s)'(g(x) - s)
*/
   else if(mfcn == VanShanno){ 
      sum1     = 0.0e0;
      sum2     = 0.0e0;

      gradf    = nlp->getGrad();
      for(i = 0; i < n ; i++)
         sum1 += gradf(i)*sk(i);

      if( constraintsExist ){  

         conresid = getConstraintResidual();

         for(i = 0; i < mi; i++){
            ind   =  n + me + mi + i;
            sum1 -= mu_*sk(ind)/s(i);
	    conresid(me+i)-=s(i);
         }

         sum2     = conresid.dot(conresid);

      }

      dirder_ = sum1 - beta_*sum2;
      if(dirder_ > 0.0){
	 if(sum2 > ftol){ 
            beta_ = 10*(fabs(sum1)/sum2);
            if(debug_){
               *optout << "||constraints||:  " << sum2 << "\n";
               *optout << "Directional directive:  " << dirder_ << "\n";
               *optout << "Beta :         " << beta_ << "\n";
            }
         }
      }
   }
}

double OptNIPSLike::dampenStep(SerialDenseVector<int,double>& step)
{  
 // Local variables
  double alphaHat, dotprod, q, minq, tau;
  const double zero = 0.0e0;
  const double one  = 1.0e0;
  int i, ind, n      = dim;

  // Compute "approach to the boundary" parameter tau
  dotprod  = s.dot(z);

  if(mfcn == NormFmu)
    tau  = max(taumin_, (1 - sw_*dotprod));
  else
     tau = taumin_;

  minq = one;  alphaHat = one;

  // Constraint Violations
  for(i = 0; i <  mi ; i++){
     ind = n + me + mi + i;
     if(s(i) <= zero){
        if(step(ind) <= zero) 
	   step(ind)  = zero;
     }
     else{
        q    = step(ind)/s(i);
	minq = min(minq, q);
     }
  }

  // Inequality Multipliers 
  for(i = 0; i <  mi ; i ++){
     ind = n + me +  i;
     if(z(i) <= zero){
        if(step(ind) <= zero) 
	   step(ind)  = zero;
     }
     else{
        q    = step(ind)/z(i);
	minq = min(minq, q);
     }
  }

  if(minq < 0.0)
     alphaHat = min(-tau/minq,one);

// Scale the step
  step*= alphaHat;
  if (debug_)
    *optout << "\n dampenStep: alphaHat = " << e(alphaHat,14,6) << flush;

  return alphaHat;
}

SerialDenseVector<int,double> OptNIPSLike::setupRHS(const SerialDenseVector<int,double>& x, double mu)
{  
  //Local variables
  int ind;
  NLP1* nlp = nlprob();

  //Local vectors
  SerialDenseVector<int,double> conresid(me+mi), szmu(mi), rhs(getGradL().length());

  rhs         = getGradL();

  if(nlp->hasConstraints()){

     conresid = getConstraintResidual();

     for(int i = 0; i < mi; i++){
        ind      = me + i;
        conresid(ind)-= s(i); // g(x) - s
        szmu(i)  = s(i)*z(i) - mu;
     }

     //rhs &= conresid;
     int omega = rhs.length();
     rhs.resize(omega+conresid.length());
     for(int i=omega;i<omega+conresid.length();i++)
       {rhs(i) = conresid(i-omega);}
	
     if(mi > 0){// rhs &= szmu;
       int psi = rhs.length();
       rhs.resize(psi+szmu.length());
       for(int i=psi;i<psi+szmu.length();i++)
	 {rhs(i) = szmu(i-psi);}
     }

  }

  return rhs;
}

SerialDenseVector<int,double> OptNIPSLike::setupRHS(const SerialDenseVector<int,double>& xt, 
                               const SerialDenseVector<int,double>& yt,
                               const SerialDenseVector<int,double>& zt,
                               const SerialDenseVector<int,double>& st,
			       double  mu)
{  
  //Local variables
  int ind;
  NLP1* nlp = nlprob();
  bool constraintsExist = nlp->hasConstraints();
  bool modeOverride = nlp->getModeOverride();

 //Local vectors
  SerialDenseVector<int,double> conresid(me+mi), szmu(mi), rhs, trhs, yzmultiplier;

 // Vector compatibility check
  if(yt.numRows() != me || zt.numRows() != mi || st.numRows()!= mi){
    if(debug_){
       *optout
        << "The equality multiplier   contains " << yt.numRows() << " elts. " 
	<< "\n";
       *optout
        << "The inequality multiplier contains " << zt.numRows() << " elts. "
	<< "\n";
       *optout
        << "The slack vector has               " << st.numRows() << " elts. "
	<< "\n";
       *optout << flush;
    }
  }

  if(constraintsExist){

    if (modeOverride)
      nlp->getConstraints()->evalCFGH(xt);

    conresid  = nlp->getConstraints()->evalResidual(xt);

    for(int i = 0; i < mi; i++){
      ind      = me + i;
      conresid(ind)-= st(i); // g(xt) - st
      szmu(i)  = st(i)*zt(i) - mu;
    }
    trhs.resize(conresid.length());
    trhs = conresid;
    if(mi > 0){// trhs &= szmu;
      int phi = trhs.length();
      trhs.resize(phi+szmu.length());
      for(int i=phi;i<phi+szmu.length();i++)
	{trhs(i) = szmu(i-phi);}
    }

  }

  // 08/16/01 PJW - Explanation of code fragment location
  //  Computing the gradient of the constraint prior to the residual
  //  leads to duplicate function evaluations for DAKOTA
  //  Concatenate the Lagrange multipliers for ease of gradient evaluation
  //
  yzmultiplier.resize(yt.length()+zt.length());
  for(int i=0;i<yt.length()+zt.length();i++)
    {if(i<yt.length())
	{yzmultiplier(i) = yt(i);}
      else{yzmultiplier(i) = zt(i-yt.length());}
    }
  //yzmultiplier = yt & zt;

  SpecOption SpecPass      = nlp->getSpecOption(); 
  // Turn the speculative gradient option off
  nlp->setSpecOption(NoSpec);
  if (modeOverride){
    nlp->setX(xt);
    nlp->eval();
  }
  SerialDenseVector<int,double> grad(xt.length());   
  grad    = nlp->evalG(xt);
  SerialDenseVector<int,double> rhstemp;
  rhs.resize(grad.length());
  rhs = grad; 
  // Reset the speculative gradient option to its original value
  nlp->setSpecOption(SpecPass);

  if(constraintsExist){
    SerialDenseMatrix<int,double> constraintGrad(xt.length(),nlp->getConstraints()->getNumOfCons());
    constraintGrad    = nlp->getConstraints()->evalGradient(xt);
     rhstemp.resize(constraintGrad.numRows());
     rhstemp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,constraintGrad,yzmultiplier,0.0);    
     // rhs -= constraintGradient*yzmultiplier;
     rhs -= rhstemp;    
     //rhs &= trhs;
     int upsilon = rhs.length();
     rhs.resize(upsilon+trhs.length());
     for(int i=upsilon;i<upsilon+trhs.length();i++)
	   {rhs(i) = trhs(i-upsilon);}
      
  }

  return rhs;
}

SerialDenseMatrix<int,double> OptNIPSLike::setupMatrix(const SerialDenseVector<int,double>& x)
{  
/* mat  = [ hessl(n,n)       -gradh(n,me)     -gradg(n,mi)    zeros(n,mi)    ]
          [ gradh'(me,n)     zeros(me,me)    zeros(me,mi)    zeros(me,mi)   ]
          [ gradg'(mi,n)     zeros(mi,me)    zeros(mi,mi)    -diag(ones(mi))]
          [ zeros(mi,n)      zeros(mi,me)    diag(s)         diag(z)        ] */


  NLP1* nlp = nlprob();
  int n     = nlp->getDim();

 // Let's build the Jacobian 
  int mdim = n + me + 2*mi;
  SerialDenseMatrix<int,double> jacobian(mdim,mdim),testit(mdim,mdim);
  SerialDenseMatrix<int,double> D(mi,mi), S(mi,mi), Z(mi,mi);
  D = 0; S = 0; Z = 0;
  for(int i = 0; i < mi ; i++){
     D(i,i) = -1;
     S(i,i) = s(i);
     Z(i,i) = z(i);
  }

  jacobian = 0;
  testit = 0;

  if((me + mi) == 0){
    // jacobian = hessl;
    for(int i=0;i<hessl.numRows();i++)
      {for(int j=0; j<=i;j++)
	  {jacobian(i,j) = hessl(i,j);}
	for(int j=i+1;j<mdim;j++)
	  {jacobian(i,j) = hessl(j,i);}
      }
  }
//   else{
//      SerialDenseMatrix<int,double>  temp = getConstraintGradient();

//      if( mi > 0){
//         // The First Row 
//        jacobian.SubMatrix(1, n, 1, n)             = hessl;
//        jacobian.SubMatrix(1, n, n+1, n+me+mi)     = -temp;

//        jacobian.SubMatrix(n+1,  n+me+mi, 1, n) = temp.t(); 
//        jacobian.SubMatrix(n+me+1, n+me+mi, n+me+mi+1, mdim ) = D;

//        // The Last Row 
//        jacobian.SubMatrix(n+me+mi+1, mdim,  n+me+1,  n+me+mi)  = S;
//        jacobian.SubMatrix(n+me+mi+1, mdim,  n+me+mi+1, mdim )  = Z;
//     }
//     else if(me > 0){
//        // The First Row 
//        jacobian.SubMatrix(1, n, 1, n)          = hessl;
//        jacobian.SubMatrix(1, n, n+1, n+me)     = -temp; 
//        jacobian.SubMatrix(n+1, n+me, 1, n)     = temp.t(); 
//     }
//   }
  else{
    SerialDenseMatrix<int,double> temp(getConstraintGradient().numRows(),getConstraintGradient().numCols());
    temp = getConstraintGradient();
    if(mi>0)
      {
	for(int i=0;i<mdim;i++)
	for(int j=0;j<mdim;j++)
	  {if(0<=j && j<n)
	      {if(0<=i && i<n)
		  {jacobian(i,j) = hessl(i,j);
		    testit(i,j) += 1;}
 		else if(n<=i && i<n+me+mi)
 		  {jacobian(i,j) = temp(j,i-n);
 		    testit(i,j) += 1;}
	      }
	    if(n<=j && j<n+me+mi)
	      {if(0<=i && i<n)
		  {jacobian(i,j) = -temp(i,j-n);
		    testit(i,j) += 1;}
		if(n+me+mi<=i && i<mdim && n+me<=j && j<n+me+mi)
		  {jacobian(i,j) = S(i-n-me-mi,j-n-me);
		    testit(i,j) += 1;}
	      }
	    if(n+me+mi<=j && j<mdim)
	      {if(n+me<=i && i<n+me+mi)
		  {jacobian(i,j) = D(i-n-me,j-n-me-mi);
		    testit(i,j) += 1;}
		else if(n+me+mi<=i && i<mdim)
		  {jacobian(i,j) = Z(i-n-me-mi,j-n-me-mi);
		    testit(i,j) += 1;}
	      }
	  }
      } else if(me>0)
      {
	for(int i=0;i<n+me;i++)
	  {for(int j=0;j<n+me;j++)
	      {if(0<=j && j<n)
		  {if(0<=i && i<n)
		      {jacobian(i,j) = hessl(i,j);}
		  
		else if(n<=i && i<n+me)
      		  {jacobian(i,j) = temp(j,i-n);}
		  }
	       else if(n<=j && j<n+me)
		 {if(0<=i && i<n)
		     {jacobian(i,j) = -temp(i,j-n);}
		 }
	      }
	  }
      }

    
  }	  
  if (debug_) {
    Print(hessl);
    Print(jacobian);
  }

  return jacobian;
}

void OptNIPSLike::printStatus(char *title) // set Message
{
  NLP1* nlp = nlprob();

  *optout << "\n\n=========  " << title << "  ===========\n\n";
  *optout << "Optimization method       = " << method   << "\n";
  *optout << "Dimension of the problem  = " << nlp->getDim()  << "\n";
  *optout << "No. equalities            = " << me      << "\n";
  *optout << "No. inequalities          = " << mi      << "\n";
  *optout << "Merit Function (0= NormFmu, 1 = Argaez, 2 = Vanderbei) = " 
         << mfcn      << "\n";
  *optout << "Return code               = " << ret_code << " (" 
         << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << nlp->getFevals() << "\n";
  *optout << "No. gradient evaluations  = " << nlp->getGevals() << "\n";
  *optout << "No. backtracks in lnsrch  = " << backtracks << "\n";
  *optout << flush;

  if (debug_) {
    *optout << "\nHessian of the Lagrangian";
    FPrint(optout, hessl);
    
//  Compute eigenvalues of Hessian
    *optout << "Now computing eigenvalues of Hessian " << "\n";
Teuchos::LAPACK<int,double> lapack;
    int INFO;
    int LWORK = max(1,3*(hessl.numRows())-1);
    SerialDenseVector<int,double> D(hessl.numRows()), WORK(LWORK);;
    // SVD(hessl, D);
    lapack.SYEV('N','L',hessl.numRows(),hessl.values(),hessl.numRows(),D.values(),WORK.values(),3*(hessl.numRows())-1,&INFO);
    *optout << "\nEigenvalues of Hessian";
    FPrint(optout, D);
  }

  nlp->fPrintState(optout, title);
  fPrintMultipliers(optout, title);

  tol.printTol(optout);

}

void OptNIPSLike::reset()
{
  NLP1* nlp = nlprob();
  int n     = nlp->getDim();
  nlp->reset();
  OptimizeClass::defaultReset(n);
  me = mi = 0;
}

void OptNIPSLike::readOptInput() // Read opt.input file if it exists
{
  NLP1* nlp = nlprob();

/* A VERY simple routine for reading the optimization parameters
 * We should really make this more general, but as a first pass this
 * will have to do.
 * 
 * The input file should be of the form keyword = value
 * where keyword is one of the following
 * 
 * search      = linesearch
 * diff_option = forward
 * max_iter    = 100
 * maxfeval    = 1000
 * grad_tol    = 1.e-6
 * fcn_tol     = 1.e-9
 * max_step    = 100.0
 * fcn_accrcy  = 1.e-9
 * merit_fcn   = elbakry
 * hess_option = finitediff
 *
 */
  
  int  index, max_iter, max_feval;
  double grad_tol, con_tol, fcn_tol, max_step, fcn_accrcy, sigma_val, tau_val;

  char token[80], ignore[80], equals[1];
//
// Keywords allowed
//
  string keyword;
  string cdebug("debug");
  string cdiff_option("diff_option");
  string cmerit_fcn("merit_fcn");
  string cfcn_accrcy("fcn_accrcy");
  string cfcn_tol("fcn_tol");
  string ccon_tol("con_tol");
  string cgrad_tol("grad_tol");
  string cmaxfeval("maxfeval");
  string cmaxiter("max_iter");
  string cmax_step("max_step");
  string csearch("search");
  string csteplength_to_bdry("tau");
  string ccentering_parameter("sigma");

  string diff_option, debug_flag;
  string merit_fcn,   search;
  SearchStrategy search_strategy = LineSearch;
  MeritFcn    merit = ArgaezTapia;

  int keyword_count = 0;

// 
// Default name of input file
//
  const char *opt_input  = {"opt.input"};

//
// Open opt.input file and check to see if we succeeded
//

  ifstream optin(opt_input);
  if (!optin.rdbuf()->is_open()) {
    if (debug_) {
      *optout << "OptNIPSLike::ReadOptInput: No opt.input file found\n";
      *optout << "OptNIPSLike::ReadOptInput: Default values will be used\n";
    }
    return;
  }

  if (debug_) *optout << "OptNIPSLike::ReadOptInput: Reading opt.input file\n";

  optin >> token;

  while (!optin.eof()) {

    keyword = token;
    keyword_count++;

    if (keyword == cdiff_option) {

      optin >> equals >> token;
      diff_option = token;

      if ( diff_option == "forward")
	nlp->setDerivOption(ForwardDiff);
      else if ( diff_option == "backward")
	nlp->setDerivOption(BackwardDiff);
      else if ( diff_option == "central")
	nlp->setDerivOption(CentralDiff);
    }    
    else if (keyword == cdebug) {
      optin >> equals >> token;
      debug_flag = token;
      if ( debug_flag == "true") {
	setDebug();
	nlp->setDebug();
      }
    }    
    else if (keyword == cfcn_accrcy) {
      //optin >> equals >> fcn_accrcy;
      optin >> equals >> index >> fcn_accrcy;
      nlp->setFcnAccrcy(index, fcn_accrcy);
    }    
    else if (keyword == cfcn_tol) {
      optin >> equals >> fcn_tol;
      setFcnTol(fcn_tol);
    }    
    else if (keyword == ccon_tol) {
      optin >> equals >> con_tol;
      setConTol(con_tol);
    }    
    else if (keyword == cgrad_tol) {
      optin >> equals >> grad_tol;
      setGradTol(grad_tol);
    }    
    else if (keyword == cmaxfeval) {
      optin >> equals >> max_feval;
      setMaxFeval(max_feval);
    }    
    else if (keyword == cmaxiter) {
      optin >> equals >> max_iter;
      setMaxIter(max_iter);
    }
    else if (keyword == cmax_step) {
      optin >> equals >> max_step;
      setMaxStep(max_step);
    }
    else if (keyword == csteplength_to_bdry) {
      optin >> equals >> tau_val;
      setStepLengthToBdry(tau_val);
    }
    else if (keyword == ccentering_parameter) {
      optin >> equals >> sigma_val;
      setCenteringParameter(sigma_val);
    }
    else if (keyword == cmerit_fcn) {
      optin >> equals >> token;
      merit_fcn = token;
      if ( merit_fcn == "normfmu")
	merit = NormFmu;
      else if ( merit_fcn == "argaeztapia")
	merit = ArgaezTapia;
      else if ( merit_fcn == "vanderbei")
	merit = VanShanno;
      setMeritFcn(merit);
    }
    else if (keyword == csearch) {
      optin >> equals >> token;
      search = token;
      if ( search == "trustregion")
	search_strategy = TrustRegion;
      else if ( search == "linesearch")
	search_strategy = LineSearch;
      else if ( search == "trustpds")
	search_strategy = TrustPDS;
      setSearchStrategy(search_strategy);
    }
    else {
      *optout << "Unrecognized keyword '" << keyword << "'. "
	<< "Skipping the rest of this line\n";
      optin.getline(ignore, sizeof(ignore));
    }
  optin >> token;
  }

  *optout << "\n\n======  Summary of input file  ======\n\n";

  *optout << csearch      << " = " << search << "\n";
  *optout << cdiff_option << " = " << diff_option << "\n";
  *optout << cmerit_fcn   << " = " << merit_fcn << "\n";
  *optout << cmaxiter     << " = " << max_iter << "\n";
  *optout << cmaxfeval    << " = " << max_feval << "\n";
  *optout << cgrad_tol    << " = " << grad_tol << "\n";
  *optout << cfcn_tol     << " = " << fcn_tol << "\n";
  *optout << ccon_tol     << " = " << con_tol << "\n";
  *optout << cmax_step    << " = " << max_step << "\n";
  SerialDenseVector<int,double> fcnacc(nlp->getFcnAccrcy().length());
  fcnacc = nlp->getFcnAccrcy();
  for(int i = 0; i< fcnacc.numRows(); i++)
     *optout << cfcn_accrcy  << " = " << fcnacc(i) << "\n";
  *optout << ccentering_parameter  << " = " << sigma_val << "\n";
  *optout << csteplength_to_bdry  << " = " <<  tau_val << "\n";

  tol.printTol(optout);

}

} // namespace OPTPP
