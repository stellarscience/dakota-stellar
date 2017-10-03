#ifndef NLP_h
#define NLP_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.

 P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 ----------------------------------------------------------------------*/


#include "NLPBase.h"
#include "OptppSmartPtr.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

namespace OPTPP {

/**
 * NLP is a handle class for NLPBase.
 * This class is an interface to NLP0-2 for evaluating nonlinear functions.
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date Last modified 03/2007
 */


class NLP{

private:
  SmartPtr<NLPBase> ptr_;  ///< Pointer to an NLPBase object

public:
  /**
   * Default Constructor
   */
  NLP();
  /**
   * @param base pointer to an NLPBase object 
   */
  NLP(NLPBase* base); 

  /// Set the ith component of the vector x 
  void setX(const int i, const real& x);

  /// Set the current point
  void setX(const Teuchos::SerialDenseVector<int,double>& x);

  /// Set the function value
  void setF(const real& fx);

  /**
   * e = 1, simple backtracking linesearch is used in opt. algorithm
   *
   * e = 0, More-Thuente linesearch is used in opt. algorithm
   */
  void setIsExpensive(const int e);

  /// Set the ith component of the function accuracy 
  void setFcnAccrcy(const int i, const real& accrcy);

  /// Set the function accuracy  
  void setFcnAccrcy(const Teuchos::SerialDenseVector<int,double>& accrcy);

  /**
   * @return Problem Dimension
   */
  int  getDim()         const;

  /**
   * @return Number of function evaluations taken 
   */
  int  getFevals()      const;

  /**
   * @return   1,  Function evaluation is expensive
   *
   * @return  = 0,  Function evaluation is inexpensive
   */
  int  getIsExpensive() const;

  /**
   * @return The current value of function 
   */
  real getF()           const;

  /**
   * @return User-specified function accuracy
   */
  Teuchos::SerialDenseVector<int,double> getFcnAccrcy()   const;

  /**
   * @return The current value of x 
   */
  Teuchos::SerialDenseVector<int,double> getXc()  const;

  /**
   * @return CPU time used 
   */
  real getFcnTime()     const;

  /**
   * @return Total number of constraints
   */
  int getNumOfCons() const; 

  /**
   * @return Total number of nonlinear constraints
   */
  int getNumOfNLCons() const; 

  /**
   * @return  = true,  Problem is constrained
   *
   * @return  = false, Problem is unconstrained
   */
  bool hasConstraints() ; 

  /// Print value of constraints to the screen 
  void printConstraints(); 

  /// Set debug parameter = true
  void setDebug();

  /**
   * @return  = true,  debug output statements printed 
   *
   * @return  = false, debug output statements are not printed 
   */
  bool getDebug() const;

  /// Reset parameter values 
  void reset();  			

  /// Initialize selected function 
  void initFcn();  			

  /// Evaluate the function 
  real evalF();  			

  /// Evaluate the function at x 
  real evalF(const Teuchos::SerialDenseVector<int,double>& x);  	

  /// Evaluate the gradient 
  Teuchos::SerialDenseVector<int,double> evalG();  		

  /// Evaluate the gradient at x 
  Teuchos::SerialDenseVector<int,double> evalG(const Teuchos::SerialDenseVector<int,double>& x);

  /// Evaluate Hessian 
  Teuchos::SerialSymDenseMatrix<int,double> evalH();  		

  /// Evaluate Hessian at x 
  Teuchos::SerialSymDenseMatrix<int,double> evalH(Teuchos::SerialDenseVector<int,double>& x);

  /// Evaluate the function, gradient, and Hessian 
  void eval();  			

  /// Evaluate the constraints at x 
  Teuchos::SerialDenseVector<int,double> evalCF(const Teuchos::SerialDenseVector<int,double>& x); 

  /// Evaluate the constraint gradient at x 
  Teuchos::SerialDenseMatrix<int,double> evalCG(const Teuchos::SerialDenseVector<int,double>& x);  

  /// Evaluate the constraint Hessian at x 
  Teuchos::SerialSymDenseMatrix<int,double> evalCH(Teuchos::SerialDenseVector<int,double>& x);   

  /// Evaluate the constraint Hessian at x 
  OptppArray<Teuchos::SerialSymDenseMatrix<int,double> > evalCH(Teuchos::SerialDenseVector<int,double>& x, int darg);   
  void evalC(const Teuchos::SerialDenseVector<int,double>& x); 

  /// Print status of the nonlinear function to the screen 
  void printState(const char *); 
  /// Print status of the nonlinear function to file 
  void fPrintState(std::ostream *, const char *); 

};

} // namespace OPTPP
#endif

