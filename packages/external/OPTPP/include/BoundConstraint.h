#ifndef BoundConstraint_h
#define BoundConstraint_h

/*---------------------------------------------------------------------
  Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
  DE-AC04-94AL85000, there is a non-exclusive license for use of this 
  work by or on behalf of the U.S. Government.
  -----------------------------------------------------------------------*/

#include "ConstraintBase.h"

/**
 * Simple Bounds Class
 * The standard form representation is x >= a.
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date Last modified 02/2007
 */

namespace OPTPP {

class BoundConstraint: public ConstraintBase {

protected:
  /// Number of constraints.
  int          	numOfCons_; 	
  /// Number of variables.
  int          	numOfVars_; 	
  /// Number of finite lower bounds
  int          	nnzl_;		
  /// Number of finite upper bounds
  int          	nnzu_;		

  /// Lower bounds on the variables 
  Teuchos::SerialDenseVector<int,double> 	lower_;   	
  /// Upper bounds on the variables 
  Teuchos::SerialDenseVector<int,double> 	upper_;   	
  /// Value of the variables 
  mutable Teuchos::SerialDenseVector<int,double> 	cvalue_;   	

  /// Indicator of a fixed variable 
  BoolVector   	fixedVar_;   	
  /// Indicator of a free variable 
  BoolVector   	freeVar_;   	

  /// Denotes whether a constraint is written in standard form or not
  const BoolVector   	stdForm_; 
  /// Type of constraint
  Teuchos::SerialDenseVector<int,double>		ctype_; 	
  
  /// Index vector of finite constraints
  OptppArray<int> 	constraintMappingIndices_;  

public:
/**
 * Default Constructor
 * @see BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& lower)
 * @see BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& bound, 
 *                            const BoolVector& bdFlag)
 * @see BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& lower, 
 *                            const Teuchos::SerialDenseVector<int,double>& upper)
 */
  BoundConstraint();

/**
 * @param nc an integer argument 
 * @param lower a Teuchos::SerialDenseVector<int,double>
 * @see BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& bound, 
 *                            const BoolVector& bdFlag)
 * @see BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& lower, 
 *                            const Teuchos::SerialDenseVector<int,double>& upper)
 * @note Assumes all constraints are in the form x >= c
 */
  BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& lower);
  
/**
 * @param nc an integer argument 
 * @param bound a Teuchos::SerialDenseVector<int,double>
 * @param bdFlag a BoolVector
 * @see BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& lower)
 * @see BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& lower, 
 *                            const Teuchos::SerialDenseVector<int,double>& upper)
 * @note User must specify whether the set of constraints is in standard form.  
 */
  BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& bound, const BoolVector& bdFlag);

/**
 * @param nc an integer argument 
 * @param lower a Teuchos::SerialDenseVector<int,double>
 * @param upper a Teuchos::SerialDenseVector<int,double>
 * @see BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& lower)
 * @see BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& bound, 
 *                            const BoolVector& bdFlag)
 * @note Includes lower and upper bounds
 */
  BoundConstraint(int nc, const Teuchos::SerialDenseVector<int,double>& lower, const Teuchos::SerialDenseVector<int,double>& upper);

/**
 * Destructor
 */
  ~BoundConstraint() {;}

/**
 * @return Number of constraints.
 */
  virtual int getNumOfCons() const {return numOfCons_;}

/**
 * @return Number of variables.
 */
  virtual int getNumOfVars() const {return numOfVars_;}

/**
 * @return Lower bounds on the variables.
 */
  virtual Teuchos::SerialDenseVector<int,double>  getLower()  const {return lower_;}

/**
 * Set lower bounds on the variables.
 */
  void setLower(Teuchos::SerialDenseVector<int,double>& x) {lower_ = x;}

/**
 * @return Upper bounds on the variables.
 */
  virtual Teuchos::SerialDenseVector<int,double>  getUpper()  const {return upper_;}

/**
 * Set upper bounds on the variables.
 */
  void  setUpper(Teuchos::SerialDenseVector<int,double>& x) {upper_ = x;}

/**
 * @return Type of constraint.
 */
  virtual Teuchos::SerialDenseVector<int,double>  getConstraintType()  const {return ctype_;}

/**
 * @return Value of constraint, in this case, the current iterate.
 */
  virtual Teuchos::SerialDenseVector<int,double>  getConstraintValue()  const {return cvalue_;}

/**
 * CPJW Placeholder!!
 * @return Constraint violation 
 */
  virtual Teuchos::SerialDenseVector<int,double>  getConstraintViolation()  const {return cvalue_;}

/**
 * @return Indices of constraints with finite bounds 
 */
  OptppArray<int> getConstraintMappingIndices() const 
	  { return constraintMappingIndices_; }  

/**
 * @return Fixed variable indicator. 
 */
  BoolVector getFixedVar() const {return fixedVar_;}

/**
 * @return Free variable indicator. 
 */
  BoolVector getFreeVar()  const {return freeVar_;}

/**
 * @return Standard form representation.
 */
  BoolVector getStdForm()  const {return stdForm_;}

/**
 * Takes two arguments and returns a bool.
 * @param xc a Teuchos::SerialDenseVector<int,double>
 * @param epsilon a real argument
 * @return true - constraints are feasible 
 * @return false - constraints are infeasible 
 */
  virtual bool amIFeasible(const Teuchos::SerialDenseVector<int,double>& xc, double epsilon) const;

/**
 * Takes no arguments and returns a bool.
 * @return true - lower  < upper 
 * @return false - lower > upper 
 */
  bool amIConsistent() const;

// Evaluation Methods

/**
 * Takes one argument and returns a Teuchos::SerialDenseVector<int,double> of reals.
 * @param xc a Teuchos::SerialDenseVector<int,double>
 * @return The residuals of the constraints.
 */
  virtual Teuchos::SerialDenseVector<int,double> evalResidual(const Teuchos::SerialDenseVector<int,double>& xc) const;
  virtual void evalCFGH(const Teuchos::SerialDenseVector<int,double>& xc) const;

private:
/**
 * Takes one argument and returns a real Matrix.
 * @param xc a ColumnVector
 * @return The gradient of the constraints.
 */
  virtual Teuchos::SerialDenseMatrix<int,double> evalGradient(const Teuchos::SerialDenseVector<int,double>& xc) const;
/**
 * Takes one argument and returns a SymmetricMatrix.
 * @param xc a ColumnVector
 * @return The Hessian of the constraints.
 */
  virtual Teuchos::SerialSymDenseMatrix<int,double> evalHessian(Teuchos::SerialDenseVector<int,double>& xc) const;

/**
 * Takes one argument and returns an array of real SymmetricMatrices.
 * @param xc a ColumnVector
 * @param darg an integer argument
 * @return An array of constraint Hessians.
 */
  virtual OptppArray<Teuchos::SerialSymDenseMatrix<int,double> > evalHessian(Teuchos::SerialDenseVector<int,double>& xc, int darg) const;
};

} // namespace OPTPP
#endif
