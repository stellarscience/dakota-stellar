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
#include <cstdio>
#include <cstring>
#include <ctime>
#else
#include <stdio.h>
#include <string.h>
#include <time.h>
#endif

#include <string>

using namespace std;

#include "OptBCNewtonLike.h"
#include "cblas.h"
#include "ioformat.h"

#include <float.h>
#include "Teuchos_LAPACK.hpp"



using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

#ifdef OPTPP_HAVE_MPI
#include "mpi.h"
#endif

namespace OPTPP {

static const char* class_name = "OptBCNewtonLike";

//------------------------------------------------------------------------
//
//   Constrained Newton Base Class functions
//   Notes:
//   These functions are first declared in OptBCNewtonLike.h as
//   virtual functions for the abstract base class OptBCNewtonLike
//   Therefore we need to define them so that the derived classes
//   can be instantiated by the end user.  Of course the user
//   can provide his/her own version of these.
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//
// First the default functions
// defaultacceptStep
// defaultcomputeSearch
//
//------------------------------------------------------------------------
void OptBCNewtonLike::defaultAcceptStep(int iter, int step_type)
{
// Successful step
// Print out iteration summary and anything else
// you might want to do before the next iteration

  if (trace) 
    *optout << "\n***** OptBCNewtonLike:defaultacceptStep\n";

  NLP1* nlp = nlprob();
  int n     = nlp->getDim();

  static const char *steps[] = {"C", "D", "N", "B"};
  SerialDenseVector<int,double> xc(n), grad(n);
  double fvalue, gnorm;

  xc     = nlp->getXc();
  mem_step   = xc;
  mem_step -= xprev;;
  step_length = sqrt(mem_step.dot(mem_step));

  fvalue = nlp->getF();

  grad   = nlp->getGrad();
  gnorm  = sqrt(grad.dot(grad));
  
  if (debug_) {
    *optout << "\n\t xc \t\t\t   grad \t\t   step\n";
    for(int i=0; i<n; i++)
      *optout << i <<  e(xc(i),24,16) << e(grad(i),24,16) 
	   << e(mem_step(i),24,16) << "\n";
    *optout << "\nHessian";
    Print(Hessian);

//  Compute eigenvalues of Hessian
Teuchos::LAPACK<int,double> lapack;
    SerialDenseVector<int,double> D(n);
    //EigenValues(Hessian, D);
 int LWORK = max(1,3*n-1);
    SerialDenseVector<int,double> WORK(LWORK);
    int INFO;
    lapack.SYEV('N', 'L', n, Hessian.values(),n, D.values(),WORK.values(),3*n-1, &INFO);

    *optout << "\nEigenvalues of Hessian";
    Print(D);

    *optout << "\n***************************************";
    *optout << "***************************************\n";


  }
//  Iteration summary
// 
  if(step_type >= 0){
  *optout 
    << d(iter,5)  << " " << e(fvalue,12,4) << " " << e(gnorm,12,4) << " "
    << e(step_length,12,4) << "  " << steps[step_type] << " " 
    << d(fcn_evals,5) << " " << d(grad_evals,5) << endl;
  }
  else{
  *optout 
    << d(iter,5)  << " " << e(fvalue,12,4) << " " << e(gnorm,12,4) << " "
    << e(step_length,12,4) << "  " << " "  << " " 
    << d(fcn_evals,5) << " " << d(grad_evals,5) << endl;
  }
}

SerialDenseVector<int,double> OptBCNewtonLike::defaultComputeSearch(SerialSymDenseMatrix<int,double>& H)
{  
  NLP1*	nlp = nlprob();
  int   i, j, ncnt=0, *index_array, n = nlp->getDim();

  SerialDenseVector<int,double>          gg(n), sk2(n), sk(n);
  SerialSymDenseMatrix<int,double>       H1;
  SerialDenseMatrix<int,double> L;

  // set up index_array to count the number of free variables

  index_array = new int[n+1];
  for (i=1; i<=n; i++) index_array[i] = 0;
  for (i=1; i<=n; i++) 
    if (work_set(i-1) == false) index_array[i] = ++ncnt;
  if (ncnt != (n-nactive)) {
    *optout << "Number of fixed and free variables do not correspond. \n";
    exit(-1);
  }

  // Form the projected Hessian

  H1.reshape(ncnt);
  for (i=1; i<=n; i++) {
     for (j=1; j<=n; j++) 
       if (index_array[i] != 0 && index_array[j] != 0)  
	  H1(index_array[i]-1,index_array[j]-1) = H(i-1,j-1);
  }

  // Form the projected gradient

  gg.resize(ncnt);
  for (i=1; i<=n; i++) 
    if (index_array[i] != 0) gg(index_array[i]-1) = gprev(i-1); 

  // Solve (H1 * sk2 = - gg) for projected search direction sk2 

  L.reshape(ncnt,ncnt);
  sk2.resize(ncnt);
  if (ncnt == 1) sk2(0) = - gg(0) / H1(0,0);
  else if (ncnt > 1) {
    L   = MCholesky(H1);
    //sk2 = -(L.t().i()*(L.i()*gg));
    Teuchos::LAPACK<int,double> lapack;
    sk2 = gg;
    sk2 *= -1;
    int INFO;
    lapack.TRTRS('L','N','N',ncnt,1,L.values(),ncnt, sk2.values(),ncnt, &INFO);
    lapack.TRTRS('L','T','N',ncnt,1,L.values(),ncnt,sk2.values(),ncnt,&INFO);
    //return sk2;

    // Teuchos::SerialSpdDenseSolver<int,double> My_Solver;
    // int info = 0;
 
    // My_Solver.setMatrix(Teuchos::rcp(&H1,false));
    // My_Solver.setVectors(Teuchos::rcp(&sk2, false), Teuchos::rcp(&gg,false));
    // My_Solver.equilibrateMatrix();
    //My_Solver.equilibrateRHS();
    // info = My_Solver.solve();
    // if(info != 0)
    // {return sk2;}

  }


  // Form search direction sk from from projected search direction sk2 

  for (i=0; i<n; i++) sk(i) = 0.0;
  for (i=1; i<=n; i++) 
    if (index_array[i] != 0) sk(i-1) = sk2(index_array[i]-1);

  // Sanitation and return

  delete [] index_array; 
  return sk;
}

//------------------------------------------------------------------------
//
// Now all the other functions that can be generalized
// to all constrained Newton cases
//
// checkConvg
// CheckDeriv
// computeStep
// InitOpt
// initTrustRegionSize
// Optimize
//------------------------------------------------------------------------

//---------------------------------------------------------------------------- 

int OptBCNewtonLike::checkAnalyticFDGrad() 
{
  int i, n = dim, retcode = GOOD; 
  double eta, gnorm, maxerr, third;
  SerialDenseVector<int,double> error(n), fd_grad(n), grad(n);

  double mcheps = DBL_EPSILON;

  NLP1* nlp = nlprob();
  SerialDenseVector<int,double> xc(nlp->getXc().length());
  xc = nlp->getXc();
  double fx = nlp->getF();
  SpecOption tmpSpec = nlp->getSpecOption();

  nlp->setSpecOption(NoSpec);
  fd_grad   = nlp->FDGrad(sx, xc, fx, fd_grad); // Evaluate gradient using finite differences
  nlp->setSpecOption(tmpSpec);
  grad      = nlp->getGrad();  	// Now use the analytical functions
  third     = 0.33333;
  gnorm     = grad.normInf();
  eta       = pow(mcheps,third)*max(1.0,gnorm);

  if(debug_){
     *optout << "Check_Deriv: Checking gradients versus finite-differences\n";
     *optout << "    i    gradient     fd grad       error\n";
     for (i=0; i<n; i++) {
        error(i) = fabs(grad(i)-fd_grad(i));
        *optout << d(i,5) << e(grad(i),12,4) 
               << e(fd_grad(i),12,4) << e(error(i),12,4) << "\n";
      }
  }
  maxerr = error.normInf();           
  if(debug_){
     *optout << "maxerror = " << e(maxerr, 12,4) 
            << "tolerance =  " << e(eta, 12,4) << "\n";
  }
  if (maxerr > eta) retcode = BAD;
  return retcode;
}

int OptBCNewtonLike::checkConvg() // Check convergence
{
  NLP1* nlp = nlprob();
  SerialDenseVector<int,double> xc(nlp->getXc());
  int   i, n = nlp->getDim();

  // Test 1. step tolerance 
  double step_tol = tol.getStepTol();
  double snorm = stepTolNorm();
  //double xnorm =  Norm2(xc);
  double xnorm = sqrt(xc.dot(xc));
  double stol  = step_tol*max(1.0,xnorm);
  if (snorm  <= stol) {
    strcpy(mesg,"Algorithm converged - Norm of last step is less than step tolerance");
    *optout << "checkConvg: snorm = " << e(snorm,12,4) 
      << "  stol = " << e(stol,12,4) << "\n";
    return 1;
  }
  
  // Test 2. function tolerance
  double ftol = tol.getFTol();
  double fvalue = nlp->getF();
  double rftol = ftol*max(1.0,fabs(fvalue));
  double deltaf = fprev - fvalue;
  if (deltaf <= rftol) {
    strcpy(mesg,"Algorithm converged - Difference in successive fcn values less than tolerance");
    *optout << "checkConvg: deltaf = " << e(deltaf,12,4) 
         << "  ftol = " << e(ftol,12,4) << "\n";
    return 2;
  }

  // Test 3. gradient tolerance 
  SerialDenseVector<int,double> grad(nlp->getGrad());
  double gtol = tol.getGTol();
  double rgtol = gtol*max(1.0,fabs(fvalue));
  for (i=0; i<n; i++) if (work_set(i) == true) grad(i) = 0.0;
  //double gnorm = Norm2(grad);
  double gnorm = sqrt(grad.dot(grad));
  if (gnorm <= rgtol) {
    strcpy(mesg,"Algorithm converged - Norm of gradient is less than gradient tolerance");
    *optout << "checkConvg: gnorm = " << e(gnorm,12,4) 
      << "  gtol = " << e(rgtol, 12,4) << "\n";
    return 3;
  }
  
  // Test 4. absolute gradient tolerance 
  if (gnorm <= gtol) {
    strcpy(mesg,"Algorithm converged - Norm of gradient is less than gradient tolerance");
    *optout << "checkConvg: gnorm = " << e(gnorm,12,4) 
      << "  gtol = " << e(gtol, 12,4) << "\n";
    return 4;
  }
  
  // Nothing to report 
  return 0;
}


int OptBCNewtonLike::checkDeriv() // Check the analytic gradient with FD gradient
{return GOOD;}

//---------------------------------------------------------------------------- 
// Compute the maximum step allowed along the search direction sk
// before we hit a constraint
//--------------------------------------------------------------------------- 
double OptBCNewtonLike::computeMaxStep(SerialDenseVector<int,double> &sk)
{
  NLP1* nlp = nlprob();
  int i, n = nlp->getDim();
  double gamma=FLT_MAX, delta=FLT_MAX;
  SerialDenseVector<int,double> lower(nlp->getConstraints()->getLower().length());
  lower = nlp->getConstraints()->getLower();
  SerialDenseVector<int,double> upper(nlp->getConstraints()->getUpper());
  upper = nlp->getConstraints()->getUpper();
  SerialDenseVector<int,double> xc(nlp->getXc().length());
  xc    = nlp->getXc();

  double snorm = sqrt(sk.dot(sk));
  double feas_tol = 1.e-3;

  for (i=0; i<n; i++) {
    if (work_set(i) == false) {
      if      (sk(i) > 0.0e0) {
        delta = (upper(i)-xc(i)) / sk(i);
        if (delta <= feas_tol) {
	  if (debug_)
	    *optout << "Hit an upper constraint for variable " << i << "\n";
	}
      }
      else if (sk(i) < 0.0e0) {
        delta = (lower(i)-xc(i)) / sk(i);
        if (delta <= feas_tol) {
	  if (debug_)
	    *optout << "Hit a  lower constraint for variable " << i << "\n";
	}
      }
      gamma = min(gamma,delta);
    }
  }
  if (debug_)
    *optout << "computeMaxStep: maximum step allowed = " << gamma*snorm << "\n";
  return gamma*snorm;
}

//---------------------------------------------------------------------------- 
// 
// Compute a step along the direction sk using either a line search
// or model trust-region approach
//
//---------------------------------------------------------------------------- 
int OptBCNewtonLike::computeStep(SerialDenseVector<int,double> sk)
{
  NLP1* nlp = nlprob();
  real stp_length = 1.0;
  real stptmp;
  real lstol  = tol.getLSTol();
//  real xtol   = tol.getStepTol();
//  real gtol   = tol.getGTol();
  real stpmax = tol.getMaxStep();
  real stpmin = tol.getMinStep();
  int  step_type;
  int  itnmax = tol.getMaxBacktrackIter();

  if (trace) *optout << class_name << ": computeStep\n";

//
// Compute the maximum step allowed
//
  stptmp = computeMaxStep(sk);
  stpmax = min(stpmax, stptmp);

  if (strategy == TrustRegion) {
    SerialSymDenseMatrix<int,double> H(Hessian.numRows());
    H = Hessian;
    step_type = trustregion(nlp, optout, H, sk, sx, TR_size, stp_length, 
			    stpmax, stpmin);
  }
  else if (strategy == LineSearch) {
    step_type = linesearch(nlp, optout, sk, sx, &stp_length, stpmax, stpmin,
			   itnmax, lstol);
  }
  else if (strategy == TrustPDS) {
    SerialSymDenseMatrix<int,double> H(Hessian.numRows());
    H = Hessian;
    step_type = trustpds(nlp, optout, H, sk, sx, TR_size, stp_length, 
			    stpmax, stpmin, searchSize);
  }
  else {
    return(-1);
  }
  
  if (step_type < 0) {
    setMesg("Algorithm terminated - No longer able to compute step with sufficient decrease");
    ret_code = -1;
    setReturnCode(ret_code);
    return(ret_code);
  }
  fcn_evals   = nlp->getFevals();
  grad_evals  = nlp->getGevals();
  step_length = stp_length;
  return(step_type);
}

void OptBCNewtonLike::initHessian()
{ 
  int i;
  NLP1* nlp = nlprob();
  int ndim = nlp->getDim();

  if (WarmStart) {
    *optout << "OptBCNewtonLike::initHessian: Warm Start specified\n";
  }
  else {
    double typx, xmax, gnorm;
    SerialDenseVector<int,double> grad(ndim), xc(ndim);
    xc     = nlp->getXc();
    grad   = nlp->getGrad();
    gnorm = sqrt(grad.dot(grad));
    SerialDenseVector<int,double> D(ndim);

    // Initialize xmax, typx and D to default values
    xmax   = -1.e30; typx   =  1.0; D      =  1.0;

    for (i=0; i < ndim; i++) xmax = max(xmax,xc(i));
    if(xmax != 0.0) typx = xmax;
    if(gnorm!= 0.0) D    = gnorm/typx;
    if (debug_) {
      *optout << "OptBCNewtonLike::initHessian: gnorm0 = " << gnorm
	<< "  typx = " << typx << "\n";
    }
    Hessian = 0.0;
    for (i=0; i < ndim; i++) Hessian(i,i) = D(i);
   }
}

void OptBCNewtonLike::initOpt()
{
  double gnorm;
  NLP1* nlp = nlprob();
  int n = nlp->getDim();

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

  if(nlp->hasConstraints()){
    CompoundConstraint* constraints = nlp->getConstraints();
    SerialDenseVector<int,double> xstart(nlp->getXc().length());
    xstart = nlp->getXc();
    double feas_tol = tol.getCTol();
    bool feasible = constraints->amIFeasible(xstart, feas_tol);
    if (!feasible) {
      *optout << "OptBCNewtonLike WARNING:  Initial guess not feasible.\n"
	      << "BCNewton may be unable to make progress." << endl;
    }
  }

  if (ret_code == 0) {
    // evaluate Function, gradient and compute initial Hessian

    nlp->eval();

    xprev = nlp->getXc();
    fprev = nlp->getF();
    gprev = nlp->getGrad();
    gnorm = sqrt(gprev.dot(gprev));
  
    //  SerialSymDenseMatrix<int,double> Hk(n);
    //  Hessian = updateH(Hk,0);

    initHessian();
    setFcnScale(fprev);

    nlp->fPrintState(optout, "Initial state");

    if(strategy == TrustRegion) {
      *optout << "\n\t\t" << method << " Method with Trust Regions\n";
      TR_size = getTRSize();
      if (TR_size == 0.0) TR_size = getGradMult()*gnorm;
      *optout << "\t\t Initial Trust Region = " << e(TR_size,12,4) << "\n";
    }
    else if(strategy == TrustPDS) {
      *optout << "\n\t\t" << method << " Method with Trust Region / PDS\n";
      TR_size = getTRSize();
      if (TR_size == 0.0) TR_size = getGradMult()*gnorm;
      *optout << "\t\t Initial Trust Region = " << e(TR_size,12,4) << "\n";
    }
    else  
      *optout << "\n\t\t" << method << " Method with Line Search\n";

    *optout << "\n  Iter      F(x)       ||grad||     "
	    << "||step||      f/g\n\n"
	    << d(0,5) << " " << e(fprev,12,4) << " " << e(gnorm,12,4) << endl;
    if (debug_) {
      nlp->fPrintState(optout, "BCNewtonLike: Initial Guess");
      *optout << "xc, grad, step\n";
      for(int i=0; i<n; i++)
	*optout << i << e(xprev(i),24,16) << e(gprev(i),24,16) << "\n";
      Print(Hessian);
    }
    updateConstraints(0);
  }
}

double OptBCNewtonLike::initTrustRegionSize() const
{ 
  double init_tr;
//
// return minimum of 100||x||, tolerance default, or Maxstep
//
  init_tr = 100.0*sqrt(xprev.dot(xprev));
  init_tr = min(init_tr, tol.getTRSize());    
  init_tr = min(init_tr, tol.getMaxStep());    

  return init_tr;
}

void OptBCNewtonLike::optimize()
//---------------------------------------------------------------------------- 
//
// Given a nonlinear operator nlp find the minimizer using a
// Newton-like method
//
//---------------------------------------------------------------------------- 
{
  int k;
  int convgd = 0;
  int maxiter, maxfev, myfevals, fevals, step_type;

// Allocate local vectors 

  int n = dim;
  SerialDenseVector<int,double> sk(n);
  SerialSymDenseMatrix<int,double> Hk(n);
  NLP1* nlp = nlprob();

// Initialize iteration
// Evaluate Function, Gradient, and Hessian

  initOpt();

  if (ret_code == 0) {
    Hk = Hessian;
    maxiter = tol.getMaxIter();
    maxfev  = tol.getMaxFeval();

    // Check for convergence. Need to take into account that this is the
    // zeroth iteration
    //  convgd = objfcn.Check_Convg(tol,stat);
    //  if (convgd > 0) {
    //    stat.ret_code = convgd;
    //    return;
    //  }
  
    for (k=1; k <= maxiter; k++) {

      iter_taken = k;
      if (debug_)
	*optout << " **** OptBCNewtonLike : iteration count = " << k << "\n";

      //  Compute search direction

      try{
        sk = computeSearch(Hk);
      }
      catch(...){
	std::cout<<"OptBCNewtonLike.C"<<std::endl;
        cout << "\n Algorithm terminated - Singular Jacobian \n";
        setMesg("Algorithm terminated - Singular Jacobian");
        return;
      }

      //  attempt to take a step in the direction sk from the current point. 
      //  The default method is to use a trust region

      if ((step_type = computeStep(sk)) >= 0) {
	acceptStep(k, step_type);
	convgd    = checkConvg();
        m_nconvgd = convgd;
      }

      //  Update Constraints

      ret_code = updateConstraints(step_type);

      //  Error checking 

      if (ret_code <= 0) { // constraints have not been modified
	if (step_type<0 && convgd==0) {//not converged and cannot take a step 
	  ret_code = step_type;
          setReturnCode(ret_code);
	  *optout << "OptBCNewtonLike : cannot take a step \n";
	  return;
	} else if (convgd > 0) { // converged
          setReturnCode(convgd);
	  *optout << "OptBCNewtonLike : convergence achieved. \n";
	  return;
	}
      }

      myfevals = nlp->getFevals();

#ifdef OPTPP_HAVE_MPI

      char buffer[MPI_MAX_ERROR_STRING];
      int error, resultlen, flag;

      error = MPI_Initialized(&flag);
      if (error != MPI_SUCCESS)
	{
	  MPI_Error_string(error, buffer, &resultlen);
	  printf("\nOptNewtonLike: MPI Error - %s\n", buffer);
	  strcpy(mesg, "Algorithm aborted - MPI generated error\n");
	  ret_code = -14;
          setReturnCode(ret_code);
	}

      if (flag == 0) {
	fevals = myfevals;
      }
      else{
	error = MPI_Allreduce(&myfevals, &fevals, 1, MPI_INT, MPI_MAX,
			      MPI_COMM_WORLD);
	if (error != MPI_SUCCESS)
	  {
	    MPI_Error_string(error, buffer, &resultlen);
	    printf("\nOptNewtonLike: MPI Error - %s\n", buffer);
	    strcpy(mesg, "Algorithm aborted - MPI generated error\n");
	    ret_code = -15;
            setReturnCode(ret_code);
	  }
      }

#else

      fevals = myfevals;

#endif

      if (fevals > maxfev) break;

      //  if not converged, update the Hessian 

      if (convgd <= 0 || ret_code > 0) {
	Hessian = updateH(Hk,k);
	Hk = Hessian;
	xprev = nlp->getXc();
	fprev = nlp->getF();
	gprev = nlp->getGrad();
      }
    }

    setMesg("Algorithm terminated - Number of iterations or fevals exceeds the specified limit");
    ret_code = -4;
    setReturnCode(ret_code);
  }
}

void OptBCNewtonLike::printStatus(char *s) // set Message
{
  NLP1* nlp = nlprob();

  *optout << "\n\n=========  " << s << "  ===========\n\n";
  *optout << "Optimization method       = " << method << "\n";
  *optout << "Dimension of the problem  = " << nlp->getDim()  << "\n";
  *optout << "No. of bound constraints  = " << nlp->getDim()  << "\n";
  *optout << "Return code               = " << ret_code << " ("
       << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << nlp->getFevals() << "\n";
  *optout << "No. gradient evaluations  = " << nlp->getGevals() << "\n";

  if (debug_) {
    Print(Hessian);
//  Compute eigenvalues of Hessian
*optout << "Now computing eigenvalues of Hessian " << "\n";
Teuchos::LAPACK<int,double> lapack;
    int alpha = Hessian.numRows();
    SerialDenseVector<int,double> D(alpha);
    //SVD(Hessian, D);
int LWORK = max(1,3*alpha-1);
     SerialDenseVector<int,double> WORK(LWORK);
     int INFO;
    lapack.SYEV('N', 'L', alpha,Hessian.values(), alpha, D.values(),WORK.values(),3*alpha-1, &INFO);


    *optout << "\nEigenvalues of Hessian";
    Print(D);
  }

  nlp->fPrintState(optout, s);
  tol.printTol(optout);

}

void OptBCNewtonLike::readOptInput() // Read opt.input file if it exists
{
  NLP1* nlp = nlprob();

/* A VERY simple routine for reading the optimization parameters
 * We should really make this more general, but as a first pass this
 * will have to do.
 * 
 * The input file should be of the form keyword = value
 * where keyword is one of the following
 * 
 * search      = trustregion
 * diff_option = forward
 * max_iter    = 100
 * max_feval   = 1000
 * grad_tol    = 1.e-6
 * fcn_tol     = 1.e-9
 * max_step    = 100.0
 * fcn_accrcy  = 1.e-9
 *
 */
  
  int  index, max_iter, max_feval, backtrack_iter;
  real grad_tol,  fcn_tol, max_step, fcn_accrcy, backtrack_tol;

  char token[80], ignore[80], equals[1];
//
// Keywords allowed
//
  string keyword;
  string cdiff_option("diff_option");
  string cfcn_accrcy("fcn_accrcy");
  string cfcn_tol("fcn_tol");
  string cgrad_tol("grad_tol");
  string cmaxfeval("maxfeval");
  string cmaxiter("max_iter");
  string cmax_step("max_step");
  string csearch("search");
  string cbacktrack_iter("backtrack_iter");
  string cbacktrack_tol("backtrack_tol");

  string diff_option;
  string search;
  SearchStrategy s = TrustRegion;

  int keyword_count = 0;

// 
// default name of input file
//
  const char *opt_input  = {"opt.input"};

//
// Open opt.input file and check to see if we succeeded
//

  ifstream optin(opt_input);
  if (!optin.rdbuf()->is_open()) {
    *optout << "readOptInput: No opt.input file found\n";
    *optout << "readOptInput: default values will be used\n";
    return;
  }

  *optout << "readOptInput: Reading opt.input file\n";

  optin >> token;

  while (!optin.eof()) {

    keyword = token;
    keyword_count++;

//debug    *optout << "keyword = " << keyword << "\n";

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
    else if (keyword == cfcn_accrcy) {
      //optin >> equals >> fcn_accrcy;
      //nlp->setFcnAccrcy(fcn_accrcy);
      optin >> equals >> index >> fcn_accrcy;
      nlp->setFcnAccrcy(index, fcn_accrcy);
    }    
    else if (keyword == cfcn_tol) {
      optin >> equals >> fcn_tol;
      setFcnTol(fcn_tol);
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
    else if (keyword == csearch) {
      optin >> equals >> token;
      search = token;
      if ( search == "trustregion")
	s = TrustRegion;
      else if ( search == "linesearch")
	s = LineSearch;
      setSearchStrategy(s);
    }
    else if (keyword == cbacktrack_iter) {
      optin >> equals >> backtrack_iter;
      tol.setMaxBacktrackIter(backtrack_iter);
      *optout << cbacktrack_iter      << " = " << backtrack_iter << "\n";
    }
    else if (keyword == cbacktrack_tol) {
      optin >> equals >> backtrack_tol;
      tol.setLSTol(backtrack_tol);
      *optout << cbacktrack_tol      << " = " << backtrack_tol << "\n";
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
  *optout << cmaxiter     << " = " << max_iter << "\n";
  *optout << cmaxfeval    << " = " << max_feval << "\n";
  *optout << cgrad_tol    << " = " << grad_tol << "\n";
  *optout << cfcn_tol     << " = " << fcn_tol << "\n";
  *optout << cmax_step    << " = " << max_step << "\n";
  SerialDenseVector<int,double> fcnacc(nlp->getFcnAccrcy().length());
  fcnacc = nlp->getFcnAccrcy();
  for(int i = 0; i< fcnacc.numRows(); i++)
     *optout << cfcn_accrcy  << " = " << fcnacc(i) << "\n";

  tol.printTol(optout);

}

void OptBCNewtonLike::reset()
{
   NLP1* nlp = nlprob();
   int   n   = nlp->getDim();
   nlp->reset();
   OptimizeClass::defaultReset(n);
   nactive  = 0;
   work_set = false;
}

int OptBCNewtonLike::updateConstraints(int step_type)
{
  NLP1*       	nlp = nlprob();
  int          	n = nlp->getDim(), ret_flag=0;
  int          	i, j, j2, k, *new_active, actcnt=0, notnew;
  double       	reduced_grad_norm, feas_tol=1.0e-12;
  SerialDenseVector<int,double> 	lower(n), upper(n), xc(n), gg(n);

  // initialization

  lower      = nlp->getConstraints()->getLower();
  upper      = nlp->getConstraints()->getUpper();
  xc         = nlp->getXc();
  new_active = new int[n];

// cpjw - Use current gradient info or rather gradient at xc
  //gg = nlp->evalG(xc);

  // Add variables to the working set

  for (i=0; i<n; i++) {
    if ((fabs(upper(i)-xc(i))<feas_tol) || (fabs(lower(i)-xc(i))<feas_tol)) { 
      if (work_set(i) == false) {
        new_active[actcnt++] = i; work_set(i) = true; nactive++;
        *optout << "OptBCNewtonLike : variable added to working set : " << i << "\n";
      }
    } 
  }

  // Delete variables from the active set 
  // First compute the norm of the reduced gradient

  int    jdel=-1;
  double ggdel=0;

  gg = nlp->getGrad();
  for (i=0; i<n; i++) if(work_set(i) == true) gg(i) = 0.0;
  //reduced_grad_norm = Norm2(gg);
  reduced_grad_norm = sqrt(gg.dot(gg));
  if (m_nconvgd > 0 || step_type < 0) {
    gg = nlp->getGrad();
    ret_flag = -1;
    *optout << "OptBCNewtonLike : reduced_grad_norm = " << reduced_grad_norm << "\n";
    for (i=0; i<n; i++) {
      notnew = true;
      for (j=0; j<actcnt; j++) if (new_active[j] == i) notnew = false;
      if (work_set(i) == true && notnew) 
        if (((fabs(upper(i)-xc(i))<feas_tol) && gg(i)>feas_tol) ||
            ((fabs(lower(i)-xc(i))<feas_tol) && gg(i)<(-feas_tol))) {
	 if (fabs(gg(i)) > ggdel) {jdel = i; ggdel = fabs(gg(i)); }
	}
    }
    if (jdel != -1) {
      work_set(jdel) = false; nactive--;
      *optout << "OptBCNewtonLike : variable deleted from working set : " << jdel << "\n";
      ret_flag = 1;
    }
  }
  if (nactive > 0) *optout << "OptBCNewtonLike: Current working set  \n";
  k = 1;
  for (i=1; i<=nactive; i+=10) {
    *optout << " ----- variable index: ";
    j2 = min(i*10,nactive);
    for (j=(i-1)*10+1; j<=j2; j++) {
      while (work_set(k-1) == false) k++;
      *optout << d(k-1,6)  << "\t" << xc(k-1); k++;
    }
    *optout << "\n ";
  }
  return ret_flag;
}

} // namespace OPTPP
