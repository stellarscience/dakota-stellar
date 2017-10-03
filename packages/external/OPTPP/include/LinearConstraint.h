#ifndef LinearConstraint_h
#define LinearConstraint_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.

 P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 ----------------------------------------------------------------------*/


#include "ConstraintBase.h"

/**
 * LinearConstraint is a derived class of ConstraintBase.
 * LinearConstraint is an abstract class, which
 * provides common data and functionality 
 * to LinearEquation and LinearInequality.
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date Last modified 02/06/2007
 */

namespace OPTPP {

class LinearConstraint: public ConstraintBase {
protected:
  /// Number of constraints	
  int          	numOfCons_;	
  /// Number of variables	
  int          	numOfVars_;	
  /// Number of finite lower bounds	
  int	  	nnzl_;		
  /// Number of finite upper bounds	
  int	  	nnzu_;		
  /// Matrix representation of the constraints
  Teuchos::SerialDenseMatrix<int,double>       A_;		

  /// Matrix-vector product 
  Teuchos::SerialDenseVector<int,double> Ax_;		
  /// Lower bounds on the variables
  Teuchos::SerialDenseVector<int,double> lower_;		
  /// Upper bounds on the variables
  Teuchos::SerialDenseVector<int,double> upper_;		
  /// The value of the linear constraints 
  mutable Teuchos::SerialDenseVector<int,double> cvalue_;		
  /// The constraint violation, zero if all constraints satisfied 
  mutable Teuchos::SerialDenseVector<int,double> cviolation_;		

  /// Indices of finite bounds
  OptppArray<int> constraintMappingIndices_;
  /// Standard form representation of constraints
  bool   	stdForm_;	

public:
// Constructors
/**
 * Default Constructor
 */
  LinearConstraint();
/**
 * @param A a real Matrix
 * @see LinearConstraint(const Matrix& A, const Teuchos::SerialDenseVector<int,double>& b)
 * @see LinearConstraint(const Matrix& A, const Teuchos::SerialDenseVector<int,double>& b, 
 *                                 const bool rowFlag)
 * @see LinearConstraint(const Matrix& A, const Teuchos::SerialDenseVector<int,double>& lower,
 *                            const Teuchos::SerialDenseVector<int,double>& upper)
 */
  LinearConstraint(const Teuchos::SerialDenseMatrix<int,double>& A);
/**
 * @param A a real Matrix
 * @param b a ColumnVector
 * @see LinearConstraint(const Matrix& A);
 * @see LinearConstraint(const Matrix& A, const Teuchos::SerialDenseVector<int,double>& b, 
 *                                 const bool rowFlag)
 * @see LinearConstraint(const Matrix& A, const Teuchos::SerialDenseVector<int,double>& lower,
 *                            const Teuchos::SerialDenseVector<int,double>& upper)
 */
  LinearConstraint(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& b); 
/**
 * @param A a real Matrix
 * @param b a ColumnVector
 * @param rowFlag a bool
 * @see LinearConstraint(const Matrix& A);
 * @see LinearConstraint(const Matrix& A, const Teuchos::SerialDenseVector<int,double>& b)
 * @see LinearConstraint(const Matrix& A, const Teuchos::SerialDenseVector<int,double>& lower,
 *                            const Teuchos::SerialDenseVector<int,double>& upper)
 */
  LinearConstraint(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& b, 
                   const bool rowFlag);
/**
 * @param A a real Matrix
 * @param lower a ColumnVector
 * @param upper a ColumnVector
 * @see LinearConstraint(const Matrix& A);
 * @see LinearConstraint(const Matrix& A, const Teuchos::SerialDenseVector<int,double>& b)
 * @see LinearConstraint(const Matrix& A, const Teuchos::SerialDenseVector<int,double>& b, 
 *                                 const bool rowFlag)
 */
  LinearConstraint(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& lower,
                   const Teuchos::SerialDenseVector<int,double>& upper);

/**
 *  Destructor
 */
  virtual ~LinearConstraint(){}


/**
 * @return The number of constraints.
 */
  virtual int getNumOfCons() const {return numOfCons_;}

/**
 * @return The number of variables.
 */
  virtual int getNumOfVars() const {return numOfVars_;}

/**
 * @return The lower bounds on the variables.
 */
  virtual Teuchos::SerialDenseVector<int,double> getLower() const {return lower_;}

/**
 * @return The upper bounds on the variables.
 */
  virtual Teuchos::SerialDenseVector<int,double> getUpper() const {return upper_;}

/**
 * @return Value of the linear constraints.
 */
  virtual Teuchos::SerialDenseVector<int,double> getConstraintValue() const {return cvalue_;}

/**
 * @return Constraint violation.
 */
  virtual Teuchos::SerialDenseVector<int,double> getConstraintViolation() const {return cviolation_;}

/**
 * @return Indices of constraints with finite bounds 
 */
  OptppArray<int> getConstraintMappingIndices() const 
  		{ return constraintMappingIndices_; }

  
/**
 * Assigns a value to the constraint matrix.
 */
  void setA(Teuchos::SerialDenseMatrix<int,double> & A);

/**
 * Pure Virtual Functions
 */

/**
 * @return Type of constraint - Leqn
 */
  virtual Teuchos::SerialDenseVector<int,double> getConstraintType() const = 0;

/**
 * Takes one argument and returns a ColumnVector.
 * @param xc a ColumnVector
 * @return Matrix-vector product of A and xc.
 */
  virtual Teuchos::SerialDenseVector<int,double> evalAx(const Teuchos::SerialDenseVector<int,double>& xc) const = 0;

/**
 * Takes one argument and returns a ColumnVector.
 * @param xc a ColumnVector
 * @return The residual of the linear equations
 * evaluated at xc.
 */
  virtual Teuchos::SerialDenseVector<int,double> evalResidual(const Teuchos::SerialDenseVector<int,double>& xc) const = 0;
  virtual void evalCFGH(const Teuchos::SerialDenseVector<int,double>& xc) const = 0;

/**
 * Takes one argument and returns a real Matrix.
 * @param xc a ColumnVector
 * @return The gradient of the linear equations
 * evaluated at xc.
 */
  virtual Teuchos::SerialDenseMatrix<int,double> evalGradient(const Teuchos::SerialDenseVector<int,double>& xc) const = 0;

/**
 * Takes two arguments and returns a bool.
 * @param xc a ColumnVector
 * @param epsilon a real argument
 * @return A bool
 */
  virtual bool amIFeasible(const Teuchos::SerialDenseVector<int,double>& xc, double epsilon) const = 0;

private:
/**
 * Takes one argument and returns a SymmetricMatrix
 * @param xc a ColumnVector
 * @return The constraint Hessain evaluated at xc 
 */
  virtual Teuchos::SerialSymDenseMatrix<int,double> evalHessian(Teuchos::SerialDenseVector<int,double>& xc) const ;
/**
 * Takes two arguments and returns an array of real SymmetricMatrices.
 * @param xc a ColumnVector
 * @param darg an integer argument
 * @return An array of constraint Hessians.
 */
  virtual OptppArray<Teuchos::SerialSymDenseMatrix<int,double> > evalHessian(Teuchos::SerialDenseVector<int,double>& xc, int darg) const;

/**
 * Takes one arguments and returns a bool.
 * @param A a Matrix
 * @return A bool
 */
  bool dimMatch(Teuchos::SerialDenseMatrix<int,double>& A);

};

} // namespace OPTPP
#endif
