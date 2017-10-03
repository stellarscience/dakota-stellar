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
#include <ctime>
#else
#include <string.h>
#include <time.h>
#endif

#include "OptBaQNewton.h"
#include <float.h>
#include "cblas.h"
#include "ioformat.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;


namespace OPTPP {

//------------------------------------------------------------------------
// Notes: Some of these functions are first declared in opt.h as virtual
//        functions for the abstract base class Optimize.
//------------------------------------------------------------------------
// initOpt
// optimize
// checkInnerConvg
// checkConvg
// updateBarrierMultiplier
// compute_Barrier_Gradient
// compute_Barrier_Hessian
// computeSearch2
// computeStep
// computeMaxStep
// 2calarNewton
//------------------------------------------------------------------------

//--------------------------------------------------------------------------- 
// Initialize the other barrier parameters
//--------------------------------------------------------------------------- 
void OptBaQNewton::initOpt()
{
  // initialize mu - the multiplier for the barrier term

  //mu = 1.0e-2;
  mu = 1.0e-1; // This value was used in opt++.1.6.6

  // perform whatever that was done in OptBCNewtonLike::initOpt

  OptBCNewtonLike::initOpt();

  if (ret_code == 0) {
    // get the nonlinear problem 

    NLP1* nlp = nlprob();
    int   n = nlp->getDim();

    if (debug_)
      nlp->setDebug();

    // Compute the barrier function value

    double   fvalue = nlp->getF();
    SerialDenseVector<int,double> xc(nlp->getXc().length());
    xc = nlp->getXc();
    fvalue_barrier  = compute_Barrier_Fvalue(fvalue,xc);

    // Compute the barrier gradient 

    SerialDenseVector<int,double> local_grad(nlp->getGrad().length());
    local_grad = nlp->getGrad();
    grad_barrier.resize(n);
    grad_barrier = compute_Barrier_Gradient(local_grad,xc);
  }
}

//--------------------------------------------------------------------------- 
// Given a nonlinear operator nlp find the minimizer using the Newton's method
//--------------------------------------------------------------------------- 
void OptBaQNewton::optimize()
{
  int     inner_convgd, outer_convgd, step_type;
  int     inner_iter_taken, outer_iter_taken;

  // get the nonlinear problem and its dimension

  NLP1* nlp = nlprob();
  int     n = nlp->getDim();

  // declare a local search vector of length n 

  SerialDenseVector<int,double>    search_vector(n); 

  // Initialize Function , Gradient, and Hessian

  initOpt();

  if (ret_code == 0) {
    // Initialize other iteration parameters 

    outer_convgd     = false;
    outer_iter_taken = 0;
    iter_taken       = 0;

    // The main outer loop

    while (!outer_convgd) {
    
      // initialize for inner iterations

      outer_iter_taken++; 
      inner_convgd = false;
      inner_iter_taken = 0;
      fprev_outer = nlp->getF();

      while (!inner_convgd) {
	Hessian = updateH(Hessian,inner_iter_taken);
	inner_iter_taken++; 
	if (debug_) 
	  *optout << "OptBaQNewton::Optimize: iteration count = " 
		  << iter_taken << "\n";
	iter_taken++;

	// temporarily put variables aside to accommodate new data

	setAsideCurrentVariables();

	// calculate the search direction

	search_vector = computeSearch2(Hessian,grad_barrier);

	// compute the step length using quadratic-logarithmic interpolation

	step_type = computeStep(search_vector);
	if (debug_) *optout << "step_type = " << step_type << "\n";

	// if successful, accept the step; otherwise terminate inner iterations

	if (step_type < 0) inner_convgd = true;
	else {
	  acceptStep(iter_taken, step_type);
	  inner_convgd = checkInnerConvg(outer_iter_taken);
	} 
      } // while - inner loop 

      // Compute the next mu and check for convergence

      updateBarrierMultiplier();
      outer_convgd = checkConvg();

    } // while - outer loop
  }
}

//--------------------------------------------------------------------------- 
// Check for convergence in the inner iterations
//--------------------------------------------------------------------------- 
int OptBaQNewton::checkInnerConvg(int iter) 
{
  NLP1*     nlp = nlprob();
  SerialDenseVector<int,double> xc(nlp->getXc().length());
  xc = nlp->getXc();
  double      dtmp, epik, xnorm, gnorm;

  epik  = pow(10.0e0,-(iter+1.0e0));
  epik  = max(1.0e-5,epik);
  //xnorm = Norm2(xc);
  xnorm = sqrt(xc.dot(xc));
  dtmp  = max(1.0e0,xnorm);
  //gnorm = Norm2(grad_barrier);
  gnorm = sqrt(grad_barrier.dot(grad_barrier));
  dtmp  = gnorm / dtmp;
  if (debug_) *optout << "CheckInnerConvg : " << dtmp << " < " << epik << " ? \n";
  if (dtmp < epik) return true; else return false;
}

//--------------------------------------------------------------------------- 
// Check for convergence in the outer iterations
//--------------------------------------------------------------------------- 
int OptBaQNewton::checkConvg() 
{
  NLP1*      nlp = nlprob();
  SerialDenseVector<int,double> xc(nlp->getXc().length());
  xc = nlp->getXc();
  SerialDenseVector<int,double> grad(nlp->getGrad()); 
  grad = nlp->getGrad();
  SerialDenseVector<int,double> upper(nlp->getConstraints()->getUpper().length());
  upper = nlp->getConstraints()->getUpper();
  SerialDenseVector<int,double> lower(nlp->getConstraints()->getLower().length());
  lower = nlp->getConstraints()->getLower();
  double       fvalue, deltaf, rftol;
  double       gnorm, xnorm, q1, q2, qtmp;
  int          i, n = nlp->getDim();

  // Test 1. function tolerance

  if (mu < 1.0e-12){
    strcpy(mesg,"Algorithm terminated - barrier term is less than tolerance");
    return 3;
  }
  fvalue = nlp->getF();
  deltaf = fprev_outer - fvalue;
  if (deltaf == 0.0) return 0;

  rftol = 1.0e-6 * (1.0e0 + fabs(fprev));
  if (deltaf <= rftol) {
    *optout << "CheckConvg: deltaf = " << e(deltaf,12,4) 
         << " rftol = " << e(rftol,12,4) << "\n";
    return 1;
  }
  
  // Test 2. gradient tolerance 

  // xnorm = Norm2(xc);
  xnorm = sqrt(xc.dot(xc));
  for (i=0; i<n; i++) {
    if (fabs(xc(i)-lower(i))<1.0e-4 || fabs(upper(i)-xc(i))<1.0e-4) 
      grad(i) = 0.0;
  }
  //gnorm = Norm2(grad_barrier);
  gnorm = sqrt(grad_barrier.dot(grad_barrier));
  q1    = gnorm / (1.0e0 + xnorm);
  if (debug_) *optout << "CheckConvg: gnorm/(1+xnorm) = " << e(q1,12,4) << "\n"; 
  q2    = FLT_MAX;
  for (i=0; i<n; i++) {
    qtmp = xc(i) - lower(i); q2 = (qtmp < q2) ? qtmp : q2;
    qtmp = upper(i) - xc(i); q2 = (qtmp < q2) ? qtmp : q2;
  }
  q2 = - q2;
  qtmp = max(q1, q2);
  if (qtmp < 1.0e-4) {
    strcpy(mesg,"Algorithm converged - Norm or gradient is less than gradient tolerance ");
    return 2;
  } 
  return 0;
}

//--------------------------------------------------------------------------- 
// Update the Lagrange multipliers and mu
//--------------------------------------------------------------------------- 
void OptBaQNewton::updateBarrierMultiplier()
{
  NLP1*        nlp = nlprob();
  SerialDenseVector<int,double> xc(nlp->getXc().length());  
  xc = nlp->getXc();
  SerialDenseVector<int,double> lower(nlp->getConstraints()->getLower().length());
  lower = nlp->getConstraints()->getLower();
  SerialDenseVector<int,double> upper(nlp->getConstraints()->getUpper().length());
  upper = nlp->getConstraints()->getUpper();
  int          i, n  = nlp->getDim();
  double       maxmu, dtmp;

  maxmu = 10.0;
  for (i=0; i<n; i++) {
    if (lower(i) != -FLT_MAX) {
      dtmp = (xc(i) - lower(i)) / mu;
      if (dtmp < 0.0) maxmu = min(maxmu, 1.0e0/dtmp);
    }
  }
  for (i=0; i<n; i++) {
    if (upper(i) != FLT_MAX) {
      dtmp = (upper(i) - xc(i)) / mu;
      if (dtmp < 0.0) maxmu = min(maxmu, 1.0e0/dtmp);
    }
  }
  maxmu = mu / min(maxmu, 1.0e1);
  mu = maxmu;
  *optout << "UpdateBarrierMultiplier: new mu = " << mu << "\n";
}

//--------------------------------------------------------------------------- 
// Compute the barrier part of the function value
//--------------------------------------------------------------------------- 
double OptBaQNewton::compute_Barrier_Fvalue(double fcurrent, SerialDenseVector<int,double> &xc)
{
  NLP1*         nlp = nlprob();
  int            i, n = nlp->getDim();
  SerialDenseVector<int,double> upper(nlp->getConstraints()->getUpper().length());
  upper = nlp->getConstraints()->getUpper();
  SerialDenseVector<int,double> lower(nlp->getConstraints()->getLower().length());
  lower = nlp->getConstraints()->getLower();
  double        dtmp1, dtmp2, fval;

  fval = fcurrent;
  for (i=0; i<n; i++) {
    if (lower(i) != -FLT_MAX) {
      dtmp1 = xc(i) - lower(i);
      dtmp1 = log(dtmp1);
    } else dtmp1 = 0.0;
    if (upper(i) !=  FLT_MAX) {
      dtmp2 = upper(i) - xc(i);
      dtmp2 = log(dtmp2);
    } else dtmp2 = 0.0;
    fval = fval - mu * (dtmp2 + dtmp1);
  }
  return fval;
}

//--------------------------------------------------------------------------- 
// Compute the barrier part of the gradient 
//--------------------------------------------------------------------------- 
SerialDenseVector<int,double> OptBaQNewton::compute_Barrier_Gradient(SerialDenseVector<int,double> &ingrad,
						    SerialDenseVector<int,double> &xc)
{
  NLP1*         nlp = nlprob();
  int            i, n = nlp->getDim();
  SerialDenseVector<int,double> upper(nlp->getConstraints()->getUpper().length());
  upper = nlp->getConstraints()->getUpper();
  SerialDenseVector<int,double> lower(nlp->getConstraints()->getLower().length());
  lower = nlp->getConstraints()->getLower();
  SerialDenseVector<int,double>  gk(n);
  double        dtmp1, dtmp2;
  
  gk = ingrad;
  for (i=0; i<n; i++) {
    if (lower(i) != -FLT_MAX) dtmp1 = 1.0 / (xc(i) - lower(i)); else dtmp1=0.0;
    if (upper(i) !=  FLT_MAX) dtmp2 = 1.0 / (upper(i) - xc(i)); else dtmp2=0.0;
    gk(i) = gk(i) + mu * (dtmp2 - dtmp1);
  }
  return gk;
}

//--------------------------------------------------------------------------- 
// Compute the barrier part of the Hessian 
//--------------------------------------------------------------------------- 
SerialSymDenseMatrix<int,double> OptBaQNewton::compute_Barrier_Hessian(SerialSymDenseMatrix<int,double> &H,
						      SerialDenseVector<int,double> &xc)
{
  NLP1*           nlp = nlprob();
  int              i, n = nlp->getDim();
  SerialDenseVector<int,double> upper(nlp->getConstraints()->getUpper().length());  
  upper = nlp->getConstraints()->getUpper();
  SerialDenseVector<int,double> lower(nlp->getConstraints()->getLower().length());   
  lower = nlp->getConstraints()->getLower();
  double          dtmp1, dtmp2;
  SerialSymDenseMatrix<int,double> H2(n);

  H2 = H;
  for (i=0; i<n; i++) {
    if (lower(i) != -FLT_MAX) {
      dtmp1  = xc(i) - lower(i);
      dtmp1  = 1.0 / (dtmp1 * dtmp1);
    } else dtmp1 = 0.0;
    if (upper(i) != FLT_MAX) {
      dtmp2  = upper(i) - xc(i);
      dtmp2  = 1.0 / (dtmp2 * dtmp2);
    } else dtmp2 = 0.0;
    H2(i,i) = H2(i,i) + mu * (dtmp1 - dtmp2);
  }
  return H2;
}

//--------------------------------------------------------------------------- 
// Compute the Search direction 
//--------------------------------------------------------------------------- 
SerialDenseVector<int,double> OptBaQNewton::computeSearch2(SerialSymDenseMatrix<int,double> &H, SerialDenseVector<int,double> &g)
{
  NLP1*               nlp = nlprob();
  int                 n   = nlp->getDim();
  SerialDenseVector<int,double>        sk(n);
  SerialDenseMatrix<int,double> L(n,n);

  L = MCholesky(H);
  //  sk = -(L.t().i()*(L.i()*g));
  //g *= -1;
  sk = g;
  sk *= -1;
//Using LAPACK to solve two triangular systems
Teuchos::LAPACK<int,double> lapack;
    int INFO;
   lapack.TRTRS('L','N','N',n,1,L.values(),n, g.values(),n, &INFO);
   lapack.TRTRS('L','T','N',n,1,L.values(),n,g.values(),n,&INFO);
   return sk;
   // Teuchos::SerialSpdDenseSolver<int,double> My_Solver;
   // int info = 0;
 
   // My_Solver.setMatrix(Teuchos::rcp(&H,false));
   // My_Solver.setVectors(Teuchos::rcp(&sk, false), Teuchos::rcp(&gprev,false));
   // My_Solver.equilibrateMatrix();
  //My_Solver.equilibrateRHS();
   // info = My_Solver.solve();
   // if(info != 0)
   // {return sk;}

   //return sk;
}

//--------------------------------------------------------------------------- 
// Compute the step length along pk
// The algorithm goes as follow :
//  1. initialize the interval of uncertainty
//  2. Calculate alpha(0) based on feasibility and the barrier function
//        alpha(0) = max(ComputeMaxStep(pk)+mu/(g'*pk), 0.5*ComputeMaxStep(pk))
//
//  Reference : "Line search procedures for the logarithmic barrier function"
//              by Murray and Wright (SIAM J. Optimization, May 1994)
//--------------------------------------------------------------------------- 
int OptBaQNewton::computeStep(SerialDenseVector<int,double> pk)
{
  NLP1*       nlp = nlprob();
  int          n = nlp->getDim();
  double       alpha_bar, inner_gp, alpha_b, alpha_bar_plus, alpha;
  double       alpha_u, fplus, initslope, ftol, fnext, b, d, y;
  double       inner_gpnew, a, c, dtmp1, dtmp2, dtmp3;
  SerialDenseVector<int,double> gplus(n), gnext(n);
  SerialDenseVector<int,double> xc(nlp->getXc().length()),xplus(n);
  xc = nlp->getXc();

  // Initialization

  ftol     = tol.getFTol();
  alpha_u  = 1.0;

  // compute alpha_bar (max step that can be taken w/o violating constraints)

  alpha_bar = computeMaxStep(pk);
  if (debug_)
    *optout << "ComputeStep : max alpha that can be taken = " 
         << alpha_bar << "\n";

  // choose a reasonable step on the barrier function based on alpha_bar

  inner_gp = grad_barrier.dot(pk);
  alpha_bar_plus = alpha_bar + mu / inner_gp;
  if (alpha_bar < FLT_MAX && alpha_bar_plus < 0.0) 
     alpha_b = max(alpha_bar_plus,0.5*alpha_bar);
  else if (alpha_bar < FLT_MAX && alpha_bar_plus >= 0.0) 
     alpha_b = 0.95 * alpha_bar;
  else alpha_b = FLT_MAX;
  if (debug_)
    *optout << "ComputeStep : best alpha that can be taken = " 
	 << alpha_b << "\n";

  // set initial upper bound for step length 

  alpha = (alpha_b < alpha_u) ? alpha_b : alpha_u;
  if (debug_) *optout << "ComputeStep : initial alpha = " << alpha << "\n";

  // Solve quadratic-logarithmic function to find optimal alpha
  // Q(x) = a + bx + cx^2 - mu * log(d-x)

  SerialDenseVector<int,double> AnotherTemp(pk.length());
  AnotherTemp = pk.scale(alpha);
  xplus = xc;
  xplus +=  AnotherTemp;
  fnext = nlp->evalF(xplus);
  fplus = compute_Barrier_Fvalue(fnext,xplus);

  initslope = -grad_barrier.dot(grad_barrier);
  if (fplus < fvalue_barrier + initslope * ftol) {
    nlp->setX(xplus);
    nlp->setF(fnext);
    nlp->evalG();
    fcn_evals   = nlp->getFevals();
    grad_evals  = nlp->getGevals();
    step_length = alpha;
    return 0 ;
  } 

  gnext = nlp->evalG(xplus);
  gplus = compute_Barrier_Gradient(gnext,xplus);
  inner_gpnew = gplus.dot(pk);

  if (debug_) {
    *optout << "ComputeStep : fval (old, new) = " << fvalue_barrier << " " 
         << fplus << "\n";
    *optout << "ComputeStep : g'p  (old, new) = " << inner_gp << " " 
         << inner_gpnew << "\n";
  }
  y = scalarNewton(fvalue_barrier, inner_gp, fplus, inner_gpnew, alpha);
  if (debug_) *optout << "ComputeStep : y = " << y << "\n";
  if (y == 1) return -1;

  d = alpha / (1.0e0 - y);
  c = (inner_gpnew - inner_gp + mu / d - mu / (d - alpha)) / (2.0e0 * alpha);
  b = inner_gp - mu / d;
  a = fvalue_barrier + mu * log(d);
  if (debug_) 
    *optout << "ComputeStep : a,b,c,d = " << a << " " << b 
	 << " " << c << " " << d << "\n";
  dtmp1 = 2.0 * c * d - b;
  dtmp2 = dtmp1 * dtmp1 + 8.0 * c * (mu + b * d);
  dtmp2 = sqrt(dtmp2);
  dtmp3 = 4 * c;
  if (c == 0.0) {
    *optout << "ComputeStep: error - divide by 0. \n";
    return -1;
  }
  alpha = (dtmp1 - dtmp2) / dtmp3;
  if (debug_) { 
    *optout << "ComputeStep : alpha chosen    = " << alpha << "\n";
    *optout << "ComputeStep : the other alpha = " << (dtmp1+dtmp2)/dtmp3 << "\n";
  }

  // Check to see if step OK
  SerialDenseVector<int,double> YetAnotherTemp(pk.length());
  YetAnotherTemp = pk.scale(alpha);
  xplus = xc;
  xplus += YetAnotherTemp;
  //xplus = xc + pk.scale(alpha);
  fnext = nlp->evalF(xplus);
  fplus = compute_Barrier_Fvalue(fnext,xplus);
  if (fplus < fvalue_barrier + initslope * ftol) {
    nlp->setX(xplus);
    nlp->setF(fnext);
    nlp->evalG();
    fcn_evals   = nlp->getFevals();
    grad_evals  = nlp->getGevals();
    step_length = alpha;
    return 0 ;
  } else {
    setMesg("Algorithm terminated - No longer able to compute step with sufficient decrease.");
    return -1;
  }
}

//--------------------------------------------------------------------------- 
// Compute the maximum step allowed along the search direction sk
// before we hit a constraint
//--------------------------------------------------------------------------- 
double OptBaQNewton::computeMaxStep(SerialDenseVector<int,double> &sk)
{
  NLP1* nlp = nlprob();
  int i, n = nlp->getDim();
  double gamma=FLT_MAX, delta;
  SerialDenseVector<int,double> lower(nlp->getConstraints()->getLower().length());
  lower = nlp->getConstraints()->getLower();
  SerialDenseVector<int,double> upper(nlp->getConstraints()->getUpper().length());
  upper = nlp->getConstraints()->getUpper();
  SerialDenseVector<int,double> xc(nlp->getXc().length());
  xc    = nlp->getXc();

  double feas_tol = 1.e-3;

  for (i=0; i<n; i++) {
    if      (sk(i) > 0.0e0) {
      delta = (upper(i)-xc(i)) / sk(i);
      if (delta <= feas_tol) {
	if (debug_)
	  *optout << "OptBaQNewton::ComputeMaxStep: variable " << i 
		  << " hits upper constraint.\n";
      }
    } else if (sk(i) < 0.0e0) {
      delta = (lower(i)-xc(i)) / sk(i);
      if (delta <= feas_tol) {
	if (debug_)
	  *optout << "OptBaQNewton::ComputeMaxStep: variable " << i 
		  << " hits lower constraint.\n";
      }
    }
    if (delta < 0.0e0) delta = 0.0e0;
    gamma = min(gamma,delta);
  }
  if (debug_) 
    *optout << "OptBaQNewton::ComputeMaxStep: maximum step allowed = " 
         << gamma << "\n";
  return gamma;
}

//---------------------------------------------------------------------------- 
// Use Newton's method to find the root of the equation
// f(z) = ln(z) + 0.5 * (1/z - z) - kappa = 0 where
//    kappa = (0.5*alpha*(phi1'+phi2')-(phi2-phi1)) / mu
//--------------------------------------------------------------------------- 
double OptBaQNewton::scalarNewton(double phi1, double phi1_prime, double phi2, 
			          double phi2_prime, double alpha)
{
  int    convgd = 0;
  double y=1.0e-6, fval, f_prime, kappa;

  if (debug_) {
    *optout << "ScalarNewton: phi1       = " << phi1 << "\n";
    *optout << "ScalarNewton: phi1_prime = " << phi1_prime << "\n";
    *optout << "ScalarNewton: phi2       = " << phi2 << "\n";
    *optout << "ScalarNewton: phi2_prime = " << phi2_prime << "\n";
    *optout << "ScalarNewton: alpha      = " << alpha << "\n";
  }
  kappa = (0.5 * alpha * (phi1_prime + phi2_prime)- phi2 + phi1) / mu;
  if (debug_) *optout << "ScalarNewton: kappa = " << kappa << "\n";
  if (kappa <= 0.0) {
    *optout << "ScalarNewton: Error - interpolant inadequate. \n";
    return 1.0e0;
  }
  while (convgd == 0) {
    fval  = log(y) + 0.5 * (1/y - y) - kappa;
    if (fabs(fval) < 1.0e-4) convgd = 1;
    else {
      f_prime = 1/y - 1/(2*y*y) - 0.5;
      y       = y - fval / f_prime;
    } 
  } 
  if (debug_) *optout << "ScalarNewton: y, f       = " << y << " " << fval << "\n";
  return y;
} 

//--------------------------------------------------------------------------- 
// Accept the step and update the parameters
//--------------------------------------------------------------------------- 
void OptBaQNewton::acceptStep(int iter, int steptype)
{
  OptBCNewtonLike::acceptStep(iter, steptype);

  NLP1*      nlp = nlprob();
  SerialDenseVector<int,double> xc(nlp->getXc().length());
  xc  = nlp->getXc();
  SerialDenseVector<int,double> gg(nlp->getGrad().length());
  gg  = nlp->getGrad();
  double       fv  = nlp->getF();

  grad_barrier   = compute_Barrier_Gradient(gg,xc);
  fvalue_barrier = compute_Barrier_Fvalue(fv,xc);
}

//--------------------------------------------------------------------------- 
// update the previous variables
//--------------------------------------------------------------------------- 
void OptBaQNewton::setAsideCurrentVariables()
{
  NLP1* nlp = nlprob();
  xprev       = nlp->getXc();
  fprev       = nlp->getF();
  gprev       = nlp->getGrad();
  fprev_barrier = compute_Barrier_Fvalue(fprev,xprev);
  gprev_barrier = compute_Barrier_Gradient(gprev,xprev);
  fvalue_barrier = fprev_barrier;;
  grad_barrier   = gprev_barrier;
}

//--------------------------------------------------------------------------- 
SerialSymDenseMatrix<int,double> OptBaQNewton::updateH(SerialSymDenseMatrix<int,double>& Hk, int k) 
{

  double mcheps = DBL_EPSILON;
  double sqrteps = sqrt(mcheps);

  int i;

  NLP1*  nlp = nlprob();
  int nr     = nlp->getDim();
  SerialDenseVector<int,double> grad(nr), gradtmp(nr), xc(nlp->getXc().length());
  xc     = nlp->getXc();
  gradtmp   = nlp->getGrad();
  grad = compute_Barrier_Gradient(gradtmp,xc);

  SerialDenseVector<int,double> D(nr);

  // BFGS formula
  
  if (k == 0) { 
    Hessian = 0.0;
    double typx, xmax, gnorm;

    // Initialize xmax, typx and D to default values
    xmax   = -1.e30; typx   =  1.0; D      =  1.0;

    //gnorm = Norm2(grad);
    gnorm = sqrt(grad.dot(grad));

    for (i=0; i < nr; i++) xmax = max(xmax,xc(i));
    if(xmax != 0.0) typx = xmax;
    if(gnorm!= 0.0)    D = gnorm/typx;
    if (debug_) {
      *optout << "UpdateH: gnorm0 = " << gnorm
	   << " typx = " << typx << "\n";
    }
    for (i=0; i < nr; i++) Hessian(i,i) = D(i);
    return Hessian;
  }
  
  SerialDenseVector<int,double> yk(nr), sk(nr), Bsk(nr);
  SerialDenseMatrix<int,double> Htmp(nr,nr);
  
  yk = grad;
  yk -= gprev_barrier;
  sk = xc; 
  sk -=  xprev;
  
  double gts = gprev_barrier.dot(sk);
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
  
  SerialDenseVector<int,double> res(nr),tmp(Hk.numRows());
  //res = yk - Hk*sk;
  tmp.multiply( Teuchos::LEFT_SIDE, 1.0, Hk, sk, 0.0);
  res = yk;
  res -= tmp;
  if (res.normInf() <= sqrteps) {
    if (debug_) {
      *optout << "UpdateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    Hessian = Hk; return Hk;
  }
  
  //Bsk = Hk*sk;
  Bsk.multiply(Teuchos::LEFT_SIDE, 1.0, Hk, sk, 0.0);

  double sBs = sk.dot(Bsk);
  double etol = 1.e-8;

  if (sBs <= etol*snorm*snorm) {
    if (debug_) {
      *optout << "UpdateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    // D = sx.AsDiagonal()*sx.AsDiagonal();
    Hk = 0;
    for (i=0; i < nr; i++) Hk(i,i) = sx(i)*sx(i);
    Hessian = Hk; return Hk;
  }
  
  // Otherwise update the Hessian approximation

  // Htmp = - (Bsk * Bsk.t()) / sBs;
  Htmp.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, -1/sBs, Bsk, Bsk, 0.0);
  //Htmp = Htmp + (yk * yk.t()) / yts;
  SerialDenseMatrix<int,double> tmp2(yk.length(),yk.length());
  for(int i=0; i<yk.length();i++)
    for(int j=0; j<yk.length();j++)
      {tmp2(i,j) = yk(i)*yk(j)/yts;}
 Htmp += tmp2;
 for(int i=0; i<Htmp.numRows(); i++)
   for(int j=0; j<=i; j++)
     {Hk(i,j) = Hk(i,j) + Htmp(i,j);}
 //Htmp = Hk + Htmp;
 //Hk << Htmp;
 // Htmp.Release(); 
  SerialDenseVector<int,double> Bgk(nr);
  //Bgk = Hk*grad;
  Bgk.multiply(Teuchos::LEFT_SIDE, 1.0, Hk, grad, 0.0);
  double gBg = grad.dot(Bgk);
  double gg  = grad.dot(grad);
  double ckp1= gBg/gg;
  if (debug_) {
    *optout << "\nUpdateH: after update, k = " << k << "\n";
    *optout << "UpdateH: sBs  = " << sBs << "\n";
    *optout << "UpdateH: ckp1 = " << ckp1 << "\n";
  }
  Hessian = Hk;
  return Hk;
}

} // namespace OPTPP
