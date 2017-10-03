#ifndef NonLinearConstraint_h
#define NonLinearConstraint_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.

 P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 ----------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include "NLP.h"
#include "ConstraintBase.h"

/**
 * NonLinearConstraint is a derived class of ConstraintBase.
 * NonLinearConstraint provides common data and functionality
 * to NonLinearEquation and NonLinearInequality.
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date Last modified 02/06/2007
 */

namespace OPTPP {

class NonLinearConstraint: public ConstraintBase{
protected:      
  /// Pointer to an NLP object
  NLP*           nlp_;		
  /// Lower bounds on the nonlinear constraints
  Teuchos::SerialDenseVector<int,double>  lower_;		
  /// Upper bounds on the nonlinear constraints
  Teuchos::SerialDenseVector<int,double>  upper_;		
  /// Value of the nonlinear constraints
  mutable Teuchos::SerialDenseVector<int,double>  cvalue_;		
  /// Nonlinear constraint violation
  mutable Teuchos::SerialDenseVector<int,double>  cviolation_;		
  /// Number of nonlinear constraints
  int   numOfCons_;		
  /// Number of variables
  int   numOfVars_;		
  /// Number of finite lower bounds
  int	nnzl_;			
  /// Number of finite upper bounds
  int	nnzu_;			
  /// Indices of finite bounds
  OptppArray<int> constraintMappingIndices_; 
  /// Standard representation of constraint
  bool  stdForm_;		
  /// Type of constraint 
  Teuchos::SerialDenseVector<int,double> ctype_;		

public:
// Constructors
/**
 * Default Constructor
 */
  NonLinearConstraint();
/**
 * @param nlprob a pointer to NLP 
 * @param numconstraints an integer argument
 */
  NonLinearConstraint(NLP* nlprob, int numconstraints = 1);
/**
 * @param nlprob a pointer to NLP 
 * @param conFlag a bool
 * @param numconstraints an integer argument
 */
  NonLinearConstraint(NLP* nlprob, const bool conFlag, int numconstraints = 1);
/**
 * @param nlprob a pointer to NLP 
 * @param rhs ColumnVector
 * @param numconstraints an integer argument
 */
  NonLinearConstraint(NLP* nlprob, const Teuchos::SerialDenseVector<int,double>& rhs, 
                   int numconstraints = 1); 
/**
 * @param nlprob a pointer to NLP 
 * @param rhs ColumnVector
 * @param conFlag a bool
 * @param numconstraints an integer argument
 */
  NonLinearConstraint(NLP* nlprob, const Teuchos::SerialDenseVector<int,double>& rhs, 
                   const bool conFlag = true, int numconstraints = 1); 
/**
 * @param nlprob a pointer to NLP 
 * @param lower a ColumnVector
 * @param upper a ColumnVector
 * @param numconstraints an integer argument
 */
  NonLinearConstraint(NLP* nlprob, const Teuchos::SerialDenseVector<int,double>& lower, 
                     const Teuchos::SerialDenseVector<int,double>& upper, int numconstraints = 1);

#ifdef DAKOTA_OPTPP
  NonLinearConstraint(NLP* nlprob, const Teuchos::SerialDenseVector<int,double>& lower, 
                     const Teuchos::SerialDenseVector<int,double>& upper, int ne, int ni );
#endif // DAKOTA_OPTPP

 /**
  * Destructor
  */
  virtual ~NonLinearConstraint(){}

 /**
  * @return The number of constraints.
  */
  virtual int getNumOfCons() const {return numOfCons_;}

 /**
  * @return The number of variables.
  */
  virtual int getNumOfVars() const {return numOfVars_;}
        
 /**
  * @return The lower bounds of the constraints.
  */
  virtual Teuchos::SerialDenseVector<int,double> getLower() const {return lower_;}

 /**
  * @return The upper bounds of the constraints.
  */
  virtual Teuchos::SerialDenseVector<int,double> getUpper() const {return upper_;}

 /**
  * @return The value of the nonlinear constraints.
  */
  virtual Teuchos::SerialDenseVector<int,double> getConstraintValue() const {return cvalue_;}

 /**
  * @return  Nonlinear constraint violation.
  */
  virtual Teuchos::SerialDenseVector<int,double> getConstraintViolation() const {return cviolation_;}

 /**
  * @return Indices of constraints with finite bounds.
  */
  virtual OptppArray<int> getConstraintMappingIndices() const 
  		{return constraintMappingIndices_;}

#ifdef DAKOTA_OPTPP
  virtual Teuchos::SerialDenseVector<int,double> getConstraintType() const {return ctype_;}
  virtual Teuchos::SerialDenseVector<int,double> evalResidual(const Teuchos::SerialDenseVector<int,double>& xc) const ;
  virtual void evalCFGH(const Teuchos::SerialDenseVector<int,double>& xc) const ;
  virtual Teuchos::SerialDenseMatrix<int,double> evalGradient(const Teuchos::SerialDenseVector<int,double>& xc) const ;
  virtual Teuchos::SerialSymDenseMatrix<int,double> evalHessian(Teuchos::SerialDenseVector<int,double>& xc) const ;
  virtual OptppArray<Teuchos::SerialSymDenseMatrix<int,double> > evalHessian(Teuchos::SerialDenseVector<int,double>& xc, int darg) const ; 
  virtual bool amIFeasible(const Teuchos::SerialDenseVector<int,double>& xc, double epsilon) const;
#else
  /**
   * Takes one argument and returns a ColumnVector.
   * @param xc a ColumnVector
   * @return The residual of nonlinear constraints evaluated at xc.
   */
  virtual Teuchos::SerialDenseVector<int,double> evalResidual(const Teuchos::SerialDenseVector<int,double>& xc) const = 0;
  virtual void evalCFGH(const Teuchos::SerialDenseVector<int,double>& xc) const = 0;

  /**
   * Takes one argument and returns a Matrix.
   * @param xc a ColumnVector
   * @return The gradient of nonlinear constraints evaluated at xc.
   */
  virtual Teuchos::SerialDenseMatrix<int,double> evalGradient(const Teuchos::SerialDenseVector<int,double>& xc) const = 0;

  /**
   * Takes one argument and returns a SymmetricMatrix
   * @param xc a ColumnVector
   * @return The Hessian of a nonlinear constraint  evaluated at xc 
   */
  virtual Teuchos::SerialSymDenseMatrix<int,double> evalHessian(Teuchos::SerialDenseVector<int,double>& xc) const = 0;

  /**
   * Takes two arguments and returns an array of real SymmetricMatrices.
   * @param xc a ColumnVector
   * @param darg an integer argument
   * @return An array of nonlinear constraint Hessians.
   */
  virtual OptppArray<Teuchos::SerialSymDenseMatrix<int,double> > evalHessian(Teuchos::SerialDenseVector<int,double>& xc, int darg) const  = 0; 

  /**
   * Takes two arguments and returns a bool.
   * @param xc a ColumnVector
   * @param epsilon a real argument.
   * @return The feasibility of nonlinear constraints at xc.
   */
  virtual bool amIFeasible(const Teuchos::SerialDenseVector<int,double>& xc, double epsilon) const = 0;
#endif // DAKOTA_OPTPP

};

} // namespace OPTPP
#endif
