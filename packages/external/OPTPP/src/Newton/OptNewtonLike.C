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

#include "OptNewtonLike.h"
#include "cblas.h"
#include "ioformat.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"
#include "Teuchos_RCP.hpp"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;


#ifdef OPTPP_HAVE_MPI
#include "mpi.h"
#endif

namespace OPTPP {

static const char* class_name = {"OptNewtonLike"};

//------------------------------------------------------------------------
//
//   Newton-Like Base Class functions
//   Notes:
//   These functions are first declared in opt.h as
//   virtual functions for the abstract base class optimize
//   Therefore we need to define them so that the derived classes
//   can be instantiated by the end user.  Of course the user
//   can provide his/her own version of these.
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//
// First the default functions
// defaultAcceptStep
// defaultComputeSearch
//
//------------------------------------------------------------------------

void OptNewtonLike::defaultAcceptStep(int iter, int step_type)
{
// Successful step
// Print out iteration summary and anything else
// you might want to do before the next iteration

  double condh;

  if (trace) 
    *optout << "\n***** OptNewtonLike:defaultAcceptStep\n";

  NLP1* nlp = nlprob();
  int n     = nlp->getDim();

  static const char *steps[] = {"C", "D", "N", "B"};
  SerialDenseVector<int,double> xc(n), grad(n);
  double fvalue, gnorm;

  xc     = nlp->getXc();
  mem_step   = xc;
  mem_step -= xprev;
  // step_length = Norm2(mem_step);
  step_length = sqrt(mem_step.dot(mem_step));

  fvalue = nlp->getF();

  grad   = nlp->getGrad();
  // gnorm  = Norm2(grad);
  gnorm = sqrt(grad.dot(grad));
  
  if (debug_) {
    *optout << "\n\t xc \t\t\t   grad \t\t   step\n";
    for(int i=0; i<n; i++)
      *optout << i <<  e(xc(i),24,16) << e(grad(i),24,16) 
	   << e(mem_step(i),24,16) << "\n";
    *optout << "\nHessian";
    FPrint(optout, Hessian);

//  Compute eigenvalues of Hessian

   Teuchos::LAPACK<int,double> lapack;
    int alpha = Hessian.numRows();
    SerialDenseVector<int,double> D(alpha);
    char JOBZ = 'N';
    char UPLO = 'U';
    int LDA = max(1,alpha);
    int LWORK = max(1,3*alpha-1);
    SerialDenseVector<int,double> WORK(LWORK);
    int INFO;
    lapack.SYEV('N', 'L', alpha,Hessian.values(), alpha, D.values(),WORK.values(),3*alpha-1, &INFO);


    // EigenValues(Hessian, D);
    *optout << "\nEigenvalues of Hessian";
    FPrint(optout, D);
    condh = abs(D(n-1)/ D(0));

    *optout << "Reciprocal Condition Number of H = " << condh << "\n";
    *optout << "\n***************************************";
    *optout << "***************************************\n";


  }
//  Iteration summary
// 
  if(step_type >= 0){
  *optout 
    << d(iter,5)  << " " << e(fvalue,12,4) << " " << e(gnorm,12,4) << " "
    << e(step_length,12,4) << "  " << steps[step_type] << " " 
    << d(fcn_evals,5) << " " << d(grad_evals,5) << "\n" << flush;
  }
  else{
  *optout 
    << d(iter,5)  << " " << e(fvalue,12,4) << " " << e(gnorm,12,4) << " "
    << e(step_length,12,4) << "  " << "  " << " " 
    << d(fcn_evals,5) << " " << d(grad_evals,5) << "\n" << flush;
  }

}
SerialDenseVector<int,double> OptNewtonLike::defaultComputeSearch(SerialSymDenseMatrix<int,double>& H)
{  
  NLP1* nlp = nlprob();
  int n     = nlp->getDim();

  SerialDenseVector<int,double> sk(n),skhat(n);
  SerialDenseMatrix<int,double> L(n,n);
  
  L = MCholesky(H);
  
  sk = gprev;
  sk *= -1;
  //sk = -(L.t().i()*(L.i()*gprev));
   //Using LAPACK to solve two triangular systems
    int INFO;
    // std::cout<<"L.values = "<<L.values()<<std::endl;
    //std::cout<<"gprev.values = "<<gprev.values()<<std::endl;

Teuchos::LAPACK<int,double> lapack;
  
   lapack.TRTRS('L','N','N',n,1,L.values(),n, sk.values(),n, &INFO);
   lapack.TRTRS('L','T','N',n,1,L.values(),n,sk.values(),n,&INFO);


  return sk;
  // Teuchos::SerialSpdDenseSolver<int,double> My_Solver;
  // int info = 0;
 
  // My_Solver.setMatrix(Teuchos::rcp(&H,false));
  // My_Solver.setVectors(Teuchos::rcp(&sk, false), Teuchos::rcp(&gprev,false));
  // My_Solver.equilibrateMatrix();
  //My_Solver.equilibrateRHS();
  //info = My_Solver.solve();
  //if(info != 0)
  //{return sk;}

}

//------------------------------------------------------------------------
//
// Now all the other functions that can be generalized
// to all Newton cases
//
// checkAnalyticFDGrad
// checkConvg
// checkDeriv
// computeStep
// initOpt
// initHessian
// initTrustRegionSize
// optimize
// readOptInput
// reset
// printStatus
//------------------------------------------------------------------------

int OptNewtonLike::checkAnalyticFDGrad()
{
  double mcheps = DBL_EPSILON;

  int i;
  int retcode = GOOD;
  int n = dim;

  double third = 0.33333;
  SerialDenseVector<int,double> error(n);

  NLP1* nlp = nlprob();
  SerialDenseVector<int,double> xc(nlp->getXc().length());
  xc = nlp->getXc();
  double fx = nlp->getF();
  SpecOption tmpSpec = nlp->getSpecOption();
  SerialDenseVector<int,double> fd_grad(n);
  nlp->setSpecOption(NoSpec);
  fd_grad = nlp->FDGrad(sx, xc, fx, fd_grad);    // evaluate gradient using finite differences
  nlp->setSpecOption(tmpSpec);
  SerialDenseVector<int,double> grad(nlp->getGrad());        // Now use the analytical functions

  double gnorm = grad.normInf();
  double eta   = pow(mcheps,third)*max(1.0,gnorm);

  *optout << "checkDeriv: checking gradients versus finite-differences\n";
  *optout << "    i    gradient     fd grad       error\n";
  for (i=0; i<n; i++) {
    error(i) = fabs(grad(i)-fd_grad(i));
    *optout << d(i,5) << e(grad(i),12,4)
           << e(fd_grad(i),12,4) << e(error(i),12,4);
  }
  double maxerr = error.normInf();
  *optout << "maxerror = " << e(maxerr, 12,4)
    << "tolerance =  " << e(eta, 12,4) << "\n";
  if (maxerr > eta) retcode = BAD;

  return retcode;
}

int OptNewtonLike::checkConvg() // check convergence
{
  NLP1* nlp = nlprob();
  SerialDenseVector<int,double> xc(nlp->getXc());

// Test 1. step tolerance 

  double step_tol = tol.getStepTol();
  double snorm = stepTolNorm();
  // double xnorm =  Norm2(xc);
  double xnorm = sqrt(xc.dot(xc));
  double stol  = step_tol*max(1.0,xnorm);
  if (snorm  <= stol) {
    strcpy(mesg,"Algorithm converged - Norm of step is less than step tolerance");
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
    strcpy(mesg,"Algorithm converged - Difference in successive fcn values is less than tolerance");
    *optout << "checkConvg: deltaf = " << e(deltaf,12,4) 
         << "  ftol = " << e(ftol,12,4) << "\n";
    return 2;
  }
  

// Test 3. gradient tolerance 

  SerialDenseVector<int,double> grad(nlp->getGrad());
  double gtol = tol.getGTol();
  double rgtol = gtol*max(1.0,fabs(fvalue));
  // double gnorm = Norm2(grad);
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

int OptNewtonLike::checkDeriv() // check the analytic gradient with FD gradient
{return GOOD;}

int OptNewtonLike::computeStep(SerialDenseVector<int,double> sk)
//---------------------------------------------------------------------------- 
// 
// Compute a step along the direction sk using either a line search
// or model trust-region approach
//
//---------------------------------------------------------------------------- 
{
  NLP1* nlp = nlprob();
  double stp_length = 1.0;
  double lstol  = tol.getLSTol();
  double stpmax = tol.getMaxStep();
  double stpmin = tol.getMinStep();
  int  step_type;
  int  itnmax = tol.getMaxBacktrackIter();

  if (trace) *optout << class_name << ": ComputeStep\n";

  if (strategy == TrustRegion) {
    SerialSymDenseMatrix<int,double> H(Hessian.numRows());
    H = Hessian;
    step_type = trustregion(nlp, optout, H, sk, sx, TR_size, stp_length, 
			    stpmax, stpmin);
    if (step_type < 0)
      Hessian = H;
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

void OptNewtonLike::initOpt()
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
  ret_code = 0;

  if (debug_)
    nlp->setDebug();

  if(nlp->hasConstraints()){
    cerr << "Error: Newton's method does not support bound, linear, or "
         << "nonlinear constraints.\n       Please select a different method "  
         << "(OptFDNIPS, OptQNIPS, OptNIPS)\n       for constrained problems." 
         << endl; 
    abort_handler(-1); 
  }

  if (ret_code == 0) {
    // evaluate Function, gradient and compute initial Hessian

    nlp->eval();

    xprev = nlp->getXc();
    fprev = nlp->getF();
    gprev = nlp->getGrad();
    // gnorm = Norm2(gprev);
    gnorm = sqrt(gprev.dot(gprev));
    //  SerialSymDenseMatrix<int,double> Hk(n);
    //  Hessian = updateH(Hk,0);

    initHessian();
    setFcnScale(fprev);

    // get optimization parameters

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
	    << d(0,5) << " " << e(fprev,12,4) << " " << e(gnorm,12,4) << "\n";

    if (debug_) {
      nlp->fPrintState(optout, "OptNewtonLike: Initial Guess");
      *optout << "xc, grad, step\n";
      for(int i=0; i<n; i++)
	*optout << i << e(xprev(i),24,16) << e(gprev(i),24,16) << "\n";
      FPrint(optout, Hessian);
    }
  }
}
void OptNewtonLike::initHessian()
{ 
  int i;
  NLP1* nlp = nlprob();
  int ndim = nlp->getDim();

  if (WarmStart) {
    *optout << "OptNewtonlike::initHessian: Warm Start specified\n";
  }
  else {
    double typx, xmax, gnorm;
    SerialDenseVector<int,double> grad(ndim), xc(ndim);
    xc     = nlp->getXc();
    grad   = nlp->getGrad();
    // gnorm  = Norm2(grad);
    gnorm = sqrt(grad.dot(grad));
    SerialDenseVector<int,double> D(ndim);

    // Initialize xmax, typx and D to default values
    xmax   = -1.e30; typx   =  1.0; D      =  1.0;
    
    for (i=0; i < ndim; i++) xmax = max(xmax,fabs(xc(i)));
    if( xmax != 0.0) typx = xmax;
    if( gnorm!= 0.0) D = gnorm/typx;
    if (debug_) {
      *optout << "OptNewtonlike::initHessian: gnorm0 = " << gnorm
	<< "  typx = " << typx << "\n";
    }
    //std::cout<<"D = "<<D<<std::endl;
    Hessian = 0.0;
    for (i=0; i < ndim; i++) Hessian(i,i) = D(i);
    //std::cout<<"Hessian After D "<<Hessian<<std::endl;
   }
}
double OptNewtonLike::initTrustRegionSize() const
{ 
  double init_tr;
//
// return minimum of 100||x||, tolerance default, or Maxstep
//
  // init_tr = 100.0*Norm2(xprev);
  init_tr = 100.0*sqrt(xprev.dot(xprev)); 
 init_tr = min(init_tr, tol.getTRSize());    
  init_tr = min(init_tr, tol.getMaxStep());    

  return init_tr;
}

void OptNewtonLike::optimize()
//---------------------------------------------------------------------------- 
//
// Given a nonlinear operator nlp find the minimizer using a
// Newton-like method
//
//---------------------------------------------------------------------------- 
{
  int k;
  int maxiter, maxfev, myfevals, fevals;
  int convgd = 0;
  int step_type;

// Allocate local vectors 

  int n = dim;
  SerialDenseVector<int,double> sk(n);
  SerialSymDenseMatrix<int,double> Hk(n);

// Initialize iteration
// evaluate Function, Gradient, and Hessian

  initOpt();

  if (ret_code == 0) {
    maxiter = tol.getMaxIter();
    maxfev  = tol.getMaxFeval();

    Hk = Hessian;

    // check for convergence. Need to take into account that this is the
    // zeroth iteration
    //  convgd = objfcn.check_Convg(tol,stat);
    //  if (convgd > 0) {
    //    stat.ret_code = convgd;
    //    return;
    //  }
  
    for (k=1; k <= maxiter; k++) {

      iter_taken = k;

      //  Solve for the Newton direction
      //  H * step = -grad;

      sk = computeSearch(Hk);

      //  ComputeStep will attempt to take a step in the direction sk 
      //  from the current point. 
      //  The default method is to use a trust region

      if ((step_type = computeStep(sk)) < 0) {
	*optout << "step_type = " << step_type << "\n";
	setMesg("Algorithm terminated - No longer able to compute step with sufficient decrease");
	ret_code = step_type;
        setReturnCode(ret_code);
	return;
      }

      //  Accept this step and update the nonlinear model

      acceptStep(k, step_type);

      //  Test for Convergence

      convgd = checkConvg();

     if (convgd > 0) {
	ret_code = convgd;
        setReturnCode(ret_code);
	return;
      }

      NLP1* nlp = nlprob();
      myfevals = nlp->getFevals();

#ifdef OPTPP_HAVE_MPI

      char buffer[MPI_MAX_ERROR_STRING];
      int error, resultlen, flag;

      // Check to see if MPI has been initialized.

      error = MPI_Initialized(&flag);
      if (error != MPI_SUCCESS)
	{
	  MPI_Error_string(error, buffer, &resultlen);
	  printf("\nOptNewtonLike: MPI Error - %s\n", buffer);
	  strcpy(mesg, "Algorithm aborted - MPI generated error\n");
	  ret_code = -14;
	  setReturnCode(ret_code);
	}

      // If it has, obtain the MAX # of fevals per processor via a
      // REDUCE operation in order to check stopping criteria.

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

      // Update state
      Hessian = updateH(Hk,k);
      Hk = Hessian;

      xprev = nlp->getXc();
      fprev = nlp->getF();
      gprev = nlp->getGrad();

      updateModel(k, n, xprev);
   }

    setMesg("Algorithm terminated - Number of iterations or fevals exceeds the specified limit");
    ret_code = -4;
    setReturnCode(ret_code);
  }
}

void OptNewtonLike::printStatus(char *s) // set Message
{
  NLP1* nlp = nlprob();

  *optout << "\n\n=========  " << s << "  ===========\n\n";
  *optout << "Optimization method       = " << method << "\n";
  *optout << "Dimension of the problem  = " << nlp->getDim()  << "\n";
  *optout << "Return code               = " << ret_code << " ("
       << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << nlp->getFevals() << "\n";
  *optout << "No. gradient evaluations  = " << nlp->getGevals() << "\n";

  if (debug_) {
    *optout << "\nHessian";
    FPrint(optout, Hessian);
//  Compute eigenvalues of Hessian
    *optout << "Now computing eigenvalues of Hessian " << "\n";
     int alpha = Hessian.numRows();
     SerialDenseVector<int,double> D(alpha);
     int LWORK = max(1,3*alpha-1);
     SerialDenseVector<int,double> WORK(LWORK);
     int INFO;
//     lapack.SYEV(JOBZ, UPLO, alpha, Hessian, LDA, D, WORK, LWORK, INFO);
Teuchos::LAPACK<int,double> lapack;
    lapack.SYEV('N', 'L', alpha,Hessian.values(), alpha, D.values(),WORK.values(),3*alpha-1, &INFO);

    // SVD(Hessian, D);
    *optout << "\nEigenvalues of Hessian";
    FPrint(optout, D);
  }

  tol.printTol(optout);

  nlp->fPrintState(optout, s);
}

void OptNewtonLike::reset() // Reset parameters 
{
   NLP1* nlp = nlprob();
   int   n   = nlp->getDim();
   nlp->reset();
   OptimizeClass::defaultReset(n);
   grad_evals = 0;
   TR_size    = 0.0;
}

void OptNewtonLike::readOptInput() // Read opt.input file if it exists
{
  NLP1* nlp = nlprob();

/* A VERY simple routine for reading the optimization parameters
 * We should doublely make this more general, but as a first pass this
 * will have to do.
 * 
 * The input file should be of the form keyword = value
 * where keyword is one of the following
 * 
 * search      = trustregion
 * diff_option = forward
 * max_iter    = 100
 * maxfeval    = 1000
 * grad_tol    = 1.e-6
 * fcn_tol     = 1.e-9
 * max_step    = 100.0
 * fcn_accrcy  = 1.e-9
 *
 */
  
  int  index, max_iter, max_feval, backtrack_iter;
  double grad_tol,  fcn_tol, max_step, fcn_accrcy, backtrack_tol;

  char token[80], ignore[80], equals[1];
//
// Keywords allowed
//
  string keyword;
  string cdebug("debug");
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

  string diff_option, debug_flag;
  string search;
  SearchStrategy s = TrustRegion;

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
      *optout << "OptNewtonLike::ReadOptInput: No opt.input file found\n";
      *optout << "OptNewtonLike::ReadOptInput: Default values will be used\n";
    }
    return;
  }

  if (debug_) *optout << "OptNewtonLike::ReadOptInput: Reading opt.input file\n";

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
      else if ( search == "trustpds")
	s = TrustPDS;
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
  fcnacc  = nlp->getFcnAccrcy();
  for(int i=0; i < fcnacc.length(); i++)
     *optout << cfcn_accrcy  << " = " << fcnacc(i) << "\n";

  tol.printTol(optout);

}

} // namespace OPTPP
