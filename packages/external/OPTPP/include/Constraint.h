#ifndef Constraint_h
#define Constraint_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.

 P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 ----------------------------------------------------------------------*/


#include "ConstraintBase.h"
#include "OptppSmartPtr.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

/**
 * Constraint is a handle class to ConstraintBase 
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date Last modified 02/2007
 */

namespace OPTPP {

class Constraint{

private:
  /// Smart Pointer to ConstraintBase
  SmartPtr<ConstraintBase> ptr_;  

public:
/**
 * Default Constructor
 * @see Constraint(ConstraintBase* base)
 */
  Constraint();
/**
 * @param base pointer to ConstraintBase 
 * @see Constraint()
 */
  Constraint(ConstraintBase* base);


/**
 * @return Number of constraints.
 */
  int getNumOfCons() const;
/**
 * @return Number of variables.
 */
  int getNumOfVars() const;
/**
 * @return Lower bounds on the variables.
 */
  Teuchos::SerialDenseVector<int,double> getLower() const;
/**
 * @return Upper bounds on the variables.
 */
  Teuchos::SerialDenseVector<int,double> getUpper() const;
/**
 * @return Type of constraint 
 */
  Teuchos::SerialDenseVector<int,double> getConstraintType() const;
/**
 * @return Value of the constraint 
 */
  Teuchos::SerialDenseVector<int,double> getConstraintViolation() const;
/**
 * @return The constraint violation 
 */
  Teuchos::SerialDenseVector<int,double> getConstraintValue() const;
/**
 * @return Indices of constraints with finite bounds 
 */
  OptppArray<int> getConstraintMappingIndices() const;

/**
 * Takes one argument and returns a ColumnVector of reals.
 * @param xcurrent a ColumnVector 
 * @return The residuals of the constraints evaluated at xcurrent. 
 */
  Teuchos::SerialDenseVector<int,double> evalResidual(const Teuchos::SerialDenseVector<int,double>& xcurrent) const;
  void evalCFGH(const Teuchos::SerialDenseVector<int,double>& xcurrent) const;
/**
 * Takes one argument and returns a real Matrix.
 * @param xcurrent a ColumnVector 
 * @return The gradient of the constraints evaluated at xcurrent. 
 */
  Teuchos::SerialDenseMatrix<int,double> evalGradient(const Teuchos::SerialDenseVector<int,double>& xcurrent) const;
/**
 * Takes one argument and returns a SymmetricMatrix
 * @param xcurrent a ColumnVector 
 * @return The constraint Hessian evaluated at xcurrent. 
 */
  Teuchos::SerialSymDenseMatrix<int,double> evalHessian(Teuchos::SerialDenseVector<int,double>& xcurrent) const;
/**
 * Takes two arguments and returns an array of real SymmetricMatrices.
 * @param xcurrent a ColumnVector 
 * @param darg an integer 
 * @return An array of constraint Hessians evaluated at xcurrent. 
 */
  OptppArray<Teuchos::SerialSymDenseMatrix<int,double> > evalHessian(Teuchos::SerialDenseVector<int,double>& xcurrent, int darg) const;

/**
 * Takes two arguments and returns a bool.
 * @param xcurrent a ColumnVector 
 * @param epsilon a real argument 
 * @return A bool  
 */
  bool amIFeasible(const Teuchos::SerialDenseVector<int,double>& xcurrent, double epsilon) const ;

};

} // namespace OPTPP
#endif
