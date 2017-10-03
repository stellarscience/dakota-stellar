#ifndef Opt_h
#define Opt_h

/*---------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000, there is a non-exclusive license for use of this
 work by or on behalf of the U.S. Government.

 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#include <fstream>
#ifdef HAVE_STD
#include <cstdlib>
#include <cstring>
// WJB - ToDo: eradicate using directives from header files
using std::strcpy;
using std::exit;
#else
#include <string.h>
#endif

#include "globals.h"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

#include "NLP.h"
#include "NLF.h"
#include "TOLS.h"

namespace OPTPP {

//-------------------------------------------------------------------------
// Various Optimization methods and related support routines
//-------------------------------------------------------------------------

inline void abort_handler(int code)
{
  if (code > 1) // code = 2 (Cntl-C signal), 0 (normal), & -1/1 (abnormal)
    std::cout << "Signal Caught!" << std::endl;
 
  // Clean up
  std::cout << std::flush; // flush cout or ofstream redirection
  std::cerr << std::flush; // flush cerr or ofstream redirection
  std::exit(code);
}

inline void opt_default_update_model(int /*k*/, int /*dim*/, Teuchos::SerialDenseVector<int,double> /*x*/) {
  /*  clog << "opt_default_update_model: " 
       << "Iter =    "   << k 
       << ", dim  =    " << dim
       << ", x(1) =    " << x(1) << "\n"; */
}

int trustregion(NLP1*, std::ostream*, Teuchos::SerialSymDenseMatrix<int,double>&, 
		Teuchos::SerialDenseVector<int,double>&, 
		Teuchos::SerialDenseVector<int,double>&, double&, double&, double stpmax = 1.e3,
		double stpmin = 1.e-9);

int trustpds(NLP1*, std::ostream*, Teuchos::SerialSymDenseMatrix<int,double>&, 
	     Teuchos::SerialDenseVector<int,double>&, 
	     Teuchos::SerialDenseVector<int,double>&, double&, double&, double stpmax = 1.e3,
	     double stpmin = 1.e-9, int searchSize = 64);

int dogleg(NLP1*, std::ostream*, Teuchos::SerialSymDenseMatrix<int,double>&, 
	   Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseVector<int,double>&,
           Teuchos::SerialDenseVector<int,double>&, double&, double&, double);

int pdsstep(NLP1*, std::ostream*, Teuchos::SerialSymDenseMatrix<int,double>&, 
	    Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseVector<int,double>&, 
	    Teuchos::SerialDenseVector<int,double>&, double&, double&, double, double&, bool, int);

int linesearch(NLP1*, std::ostream*, Teuchos::SerialDenseVector<int,double>&, 
	       Teuchos::SerialDenseVector<int,double>&, double *,
	       double stpmax = 1.e3, double stpmin = 1.e-9,
 	       int itnmax = 5, double ftol = 1.e-4, double xtol = 2.2e-16, 
	       double gtol = 0.9);

int backtrack(NLP1*, std::ostream*, Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseVector<int,double>&, double *,
	      int itnmax = 5, double ftol = 1.e-4, 
              double stpmax = 1.e3, double stpmin = 1.e-9);

int mcsrch(NLP1*, Teuchos::SerialDenseVector<int,double>&, std::ostream*, double *,
	   int itnmax = 5, double ftol = 1.e-4, double xtol = 2.2e-16, 
	   double gtol = 0.9, double stpmax = 1.e3, 
	   double stpmin = 1.e-9);

int mcstep(double *, double *, double *, double *, double *, 
	   double *, double *, double  , double  , bool *, 
	   double, double, int *);

 Teuchos::SerialDenseMatrix<int,double> PertChol(Teuchos::SerialSymDenseMatrix<int,double>&, double, 
                              double&);
 Teuchos::SerialDenseMatrix<int,double> MCholesky(Teuchos::SerialSymDenseMatrix<int,double>&);

/**
 *
 * Opt is the Base Optimization Class
 * All other Optimization classes are derived from this one
 *
 */

class OptimizeClass {

private:
  int x_optout_fd;

protected:
  /// Dimension of the problem
  int  dim;			
  /// Various tolerances assoc. with the problem
  TOLS tol;			
  /// Diagonal Scaling Matrix for x
  Teuchos::SerialDenseVector<int,double> sx;       	
  /// Diagonal Scaling Matrix for f
  Teuchos::SerialDenseVector<int,double> sfx;      	
  /// Previous iterate
  Teuchos::SerialDenseVector<int,double> xprev;	 	
  /// Objective function value at xprev
  real         fprev;	 	
  /// Current step direction
  Teuchos::SerialDenseVector<int,double> mem_step;	
  /// Length of step direction
  real         step_length;	
  virtual real stepTolNorm() const { return step_length; }

  /// What method is being used
  char method[80];  
  /// Optional message
  char mesg[80];    
  /// Return code from Optimization class
  int  ret_code;    
  /// Number of iterations taken 
  int  iter_taken;  
  /// Number of function evaluations taken
  int  fcn_evals;   
  /// Number of bactracks in linesearch  
  int  backtracks;  
  /// Print debug statements
  bool  debug_;	    
  int  trace;
  /// Compute time per iteration
  double iter_time;
  /// Total compute time 
  double total_time;

  /// User defined function to call after each nonlinear iteration
  UPDATEFCN  update_fcn;  

  std::filebuf file_buffer;
  /// Output file 
  std::ostream *optout;
  /// Output file success
  int     optout_fd;


/**
 * Provide default implementation of reset 
 */
  virtual void defaultReset(int n) 
  { 
     sfx.resize(n);
     sx.resize(n);
     xprev.resize(n);
     sx    = 1.0; 
     sfx   = 1.0; 
     xprev = 0.0;
     fcn_evals = backtracks = 0; 
  }

/**
 * Provide default implementation of AcceptStep 
 */
  virtual void defaultAcceptStep(int, int) 
    { if (debug_) *optout << "Optimize: AcceptStep not implemented yet.\n";};

/**
 * Provide default implementation of UpdateModel
 */
  virtual void defaultUpdateModel(int k, int ndim, Teuchos::SerialDenseVector<int,double> x) 
    {update_fcn(k, ndim, x);};

/**
 * Write copyright information to the screen
 */
  void copyright()
  {
    //pdh  did a quick fix here...path should not be relative
    //     for now, took out cerr statement and put else block
    //     around the rest.  
     char str[255];
     std::ifstream in("../../include/abbrev_copyright.h");
     if (!in){
       //        std::cerr << "Cannot open input file.\n";
     }
     else {
       while (in) {
	 in.getline(str,255);
	 if(in) *optout << str<< std::endl;
       }
       in.close();
     }
  }

public:
/**
 * Default Constructor
 * @see OptimizeClass(int n)
 * @see OptimizeClass(TOLS t)
 * @see OptimizeClass(int n, TOLS t)
 */
  OptimizeClass(): x_optout_fd(-1), dim(0), debug_(0), trace(0) {
    optout = new std::ostream(&file_buffer);
    file_buffer.open("OPT_DEFAULT.out", std::ios::out);
    if (!file_buffer.is_open() || !optout->good()) {
      std::cout << "OptimizeClass:: Can't open default output file\n";
      optout_fd = 0;
    }
    update_fcn = &opt_default_update_model;
    tol.setDefaultTol();
  }

/**
 * @param n an integer argument
 */
  OptimizeClass(int n): x_optout_fd(-1), dim(n), sx(n), sfx(n), xprev(n),
    fcn_evals(0), backtracks(0), debug_(0), trace(0)      {
    optout = new std::ostream(&file_buffer);
    file_buffer.open("OPT_DEFAULT.out", std::ios::out);
    if (!file_buffer.is_open() || !optout->good()) {
      std::cout << "OptimizeClass:: Can't open default output file\n";
      optout_fd = 0;
    }
    update_fcn = &opt_default_update_model;
    sx  = 1.0; sfx = 1.0; xprev = 0.0; 
    tol.setDefaultTol(); 
  }
  
/**
 * @param t a TOLS object
 */
  OptimizeClass(TOLS t): x_optout_fd(-1), dim(0), tol(t), debug_(0), trace(0){
    optout = new std::ostream(&file_buffer);
    file_buffer.open("OPT_DEFAULT.out", std::ios::out);
    if (!file_buffer.is_open() || !optout->good()) {
      std::cout << "OptimizeClass:: Can't open default output file\n";
      optout_fd = 0;
    }
    update_fcn = &opt_default_update_model;
    sx  = 1.0; sfx = 1.0; xprev = 0.0; 
  }
  
/**
 * @param n an integer argument
 * @param t a TOLS object
 */
  OptimizeClass(int n, TOLS t): x_optout_fd(-1), dim(n), tol(t), sx(n),sfx(n),
      xprev(n), fcn_evals(0), backtracks(0), debug_(0), trace(0){
    optout = new std::ostream(&file_buffer);
    file_buffer.open("OPT_DEFAULT.out", std::ios::out);
    if (!file_buffer.is_open() || !optout->good()) {
      std::cout << "OptimizeClass:: Can't open default output file\n";
      optout_fd = 0;
    }
      update_fcn = &opt_default_update_model;
      sx  = 1.0; sfx = 1.0; xprev = 0.0;
    }

  virtual ~OptimizeClass() { cleanup(); if (optout != NULL) delete optout;}
  void  cleanup() {optout->flush();};

// set various properties

  /// Set message 
  void setMesg(const char *s)   {std::strcpy(mesg,s);}     
  /// Set method of choice
  void setMethod(const char *s) {std::strcpy(method,s);}   
  /// Set maximum steplength
  void setMaxStep(real x) {tol.setMaxStep(x);}  
  /// Set minimum steplength
  void setMinStep(real x) {tol.setMinStep(x);}  
  /// Set step tolerance
  void setStepTol(real x) {tol.setStepTol(x);}	
  /// Set function tolerance
  void setFcnTol(real x)  {tol.setFTol(x);}	
  /// Set constraint tolerance
  void setConTol(real x)  {tol.setCTol(x);}	
  /// Set gradient tolerance
  void setGradTol(real x) {tol.setGTol(x);}	
  /// Set linesearch tolerance
  void setLineSearchTol(real x) {tol.setLSTol(x);} 
  /// Set maximum outer iterations
  void setMaxIter(int k)  {tol.setMaxIter(k);}	
  /// Set maximum backtrack iterations
  void setMaxBacktrackIter(int k)  {tol.setMaxBacktrackIter(k);} 
  /// Set maximum fcn evaluations
  void setMaxFeval(int k) {tol.setMaxFeval(k);} 
  /// Set update model 
  void setUpdateModel(UPDATEFCN u) {update_fcn = u;}

  /// Set step scale
  void setXScale(Teuchos::SerialDenseVector<int,double> x)  {sx = x;}  	
  /// Set function scale
  void setFcnScale(Teuchos::SerialDenseVector<int,double> x) {sfx = x;}  	
  /// Set function scale
  void setFcnScale(double x) {sfx = x;}		

  /// Set return code
  void setReturnCode(int val) {ret_code = val;}		

  /**
   * @return Problem dimension
   */
  int          getDim()      const {return dim;}	
  /**
   * @return Iterations 
   */
  int          getIter()     const {return iter_taken;}	
  /**
   * @return Return Code
   */
  int          getReturnCode()  const {return ret_code;}	
  /**
   * @return Reason algorithm terminated 
   */
  char*        getMesg()        {return mesg;}	
  /**
   * @return Previous iterate
   */
  Teuchos::SerialDenseVector<int,double> getXPrev()    const {return xprev;}
  /**
   * @return Scaling vector used for x
   */
  Teuchos::SerialDenseVector<int,double> getXScale()   const {return sx;}
  /**
   * @return Scaling vector used for f(x)
   */
  Teuchos::SerialDenseVector<int,double> getFcnScale() const {return sfx;}
  /**
   * @return Output Filename 
   */
  std::ostream*    getOutputFile() {return optout;};

  int setOutputFile(const char *filename, int append) { 

    if (x_optout_fd == -1) {  // Change the default output file
      file_buffer.close();
      if (append)
         file_buffer.open(filename, std::ios::out|std::ios::app);
      else
         file_buffer.open(filename, std::ios::out);
      if (!file_buffer.is_open() || !optout->good()) {
	std::cout << "OptimizeClass::setOutputFile: Can't open " << filename 
		  << std::endl;
	optout_fd = 0;
      }
      else
	optout_fd = 1;
    }
    else {
      std::cout << "OptimizeClass::setOutputFile: File already attached\n";
      optout_fd = 1;
    }
    return optout_fd;
  }

  int setOutputFile(int FileDescriptor) { 

    optout_fd = FileDescriptor;
    std::cerr << "setOutputFile(int FileDescriptor) no longer supported.\n"
	      << "Please use setOutputFile(const char *filename, int append)"
	      << "or setOutputFile(std::ostream& fout)."
	      << std::endl;
    optout_fd = 0;

    return optout_fd;
  }

  int setOutputFile(std::ostream& fout) { 

    optout->rdbuf(fout.rdbuf());
    if (!optout->good()) {
      std::cout << "OptimizeClass::setOutputFile: Can't open file." 
		<< std::endl;
      optout_fd = 0;
    }
    else
      optout_fd = 1;

    return optout_fd;
  }

 /// Set debug flag to true
  void setDebug()         {debug_ = true;}   
/**
 * @return Debug parameter
 */
  bool Debug()            {return debug_;}

/**
 * @note Pure virtual functions 
 * @note Each derived class must define these functions for themselves
 */
  
  virtual void  acceptStep(int, int )  = 0;
  virtual int   checkConvg()           = 0;
  virtual Teuchos::SerialDenseVector<int,double> computeSearch(Teuchos::SerialSymDenseMatrix<int,double>& ) = 0;
  virtual void  optimize()             = 0;
  virtual void  reset()                = 0;
  virtual void  readOptInput()         = 0;
  virtual void  printStatus(char *)    = 0;
  virtual void  updateModel(int, int, Teuchos::SerialDenseVector<int,double>)       = 0;
};

} // namespace OPTPP

#endif
