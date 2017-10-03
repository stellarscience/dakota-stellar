#ifndef CompoundConstraint_h
#define CompoundConstraint_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.
 ----------------------------------------------------------------------*/


#include "Constraint.h"
#include "OptppArray.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

/**
 * CompoundConstraint is an array of constraints. 
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date   Last modified 02/2006
 */

namespace OPTPP {

class CompoundConstraint: public ConstraintBase{

private:
  OptppArray<Constraint> constraints_; 	///< Array of constraints
  int numOfSets_;			///< Number of constraint sets
  Teuchos::SerialDenseVector<int,double> lower_;		///< Lower bound on the constraints 
  Teuchos::SerialDenseVector<int,double> upper_;		///< Upper bound on the constraints

public:
 /**
  * Default Constructor
  * @see CompoundConstraint(const Constraint& c1)
  * @see CompoundConstraint(const Constraint& c1, const Constraint& c2)
  * @see CompoundConstraint(const OptppArray<Constraint>& constraints)
  * @see CompoundConstraint(const CompoundConstraint& cc)
  */
  CompoundConstraint(); 
 /**
  * @param c1 a Constraint
  * @see CompoundConstraint(const Constraint& c1, const Constraint& c2)
  * @see CompoundConstraint(const OptppArray<Constraint>& constraints)
  * @see CompoundConstraint(const CompoundConstraint& cc)
  */
  CompoundConstraint(const Constraint& c1); 
 /**
  * @param c1 a Constraint
  * @param c2 a Constraint
  * @see CompoundConstraint(const Constraint& c1)
  * @see CompoundConstraint(const OptppArray<Constraint>& constraints)
  * @see CompoundConstraint(const CompoundConstraint& cc)
  */
  CompoundConstraint(const Constraint& c1, const Constraint& c2); 
 /**
  * @param constraints an array of Constraints 
  * @see CompoundConstraint(const Constraint& c1)
  * @see CompoundConstraint(const Constraint& c1, const Constraint& c2)
  * @see CompoundConstraint(const CompoundConstraint& cc)
  */
  CompoundConstraint(const OptppArray<Constraint>& constraints); 
 /**
  * Copy constructor
  * @param cc a CompoundConstraint
  * @see CompoundConstraint(const Constraint& c1)
  * @see CompoundConstraint(const Constraint& c1, const Constraint& c2)
  * @see CompoundConstraint(const OptppArray<Constraint>& constraints)
  */
  CompoundConstraint(const CompoundConstraint& cc);

 /**
  * Destructor
  */
  ~CompoundConstraint() {;}

  /// Assignment Operator
  CompoundConstraint& operator=(const CompoundConstraint& cc);

  /**
   * @return Integer 
   */
  int compare(const Constraint& c1, const Constraint& c2); 

  /**
   * @return Number of different constraint sets.
   */
  int getNumOfSets() const {return numOfSets_;}

  /**
   * @return Number of nonlinear constraints.
   */
  int getNumOfNLCons() const;

  /**
   * @return Number of total constraints.
   */
  virtual int getNumOfCons() const;

  /**
   * @return Number of variables.
   */
  virtual int getNumOfVars() const;

  /**
   * @return Lower bounds on the constraints 
   */
  virtual Teuchos::SerialDenseVector<int,double> getLower() const;

  /**
   * @return Upper bounds on the constraints 
   */
  virtual Teuchos::SerialDenseVector<int,double> getUpper() const;

  /**
   * @return Constraint type 
   */
  virtual Teuchos::SerialDenseVector<int,double> getConstraintType() const;

  /**
   * @return Value of the entire constraint set 
   */
  virtual Teuchos::SerialDenseVector<int,double> getConstraintValue() const;

  /**
   * @return Violation of the entire constraint set 
   */
  virtual Teuchos::SerialDenseVector<int,double> getConstraintViolation() const;

  /**
   * @return Value of the nonlinear constraints only 
   */
  Teuchos::SerialDenseVector<int,double> getNLConstraintValue() const;

  /**
   * @return Indices of the constraints with finite bounds 
   */
  virtual OptppArray<int> getConstraintMappingIndices() const;

  /**
   * @return Reference to constraint i
   */
  Constraint& operator[]( int i ) { return constraints_[i];}

  /**
   * @return Const reference to constraint i
   */
  const Constraint& operator[](int i) const {return constraints_[i];}

  void computeDistanceToBounds(Teuchos::SerialDenseVector<int,double>& xcurrent, Teuchos::SerialDenseVector<int,double>& d_lower, Teuchos::SerialDenseVector<int,double>& d_upper);

  /**
   * @return Feasible vector with respect to bounds 
   */
  void computeFeasibleBounds(Teuchos::SerialDenseVector<int,double>& xcurrent, double epsilon);

  /**
   * @return Feasible vector with respect to linear and nonlinear inequalities 
   */
  void computeFeasibleInequalities(Teuchos::SerialDenseVector<int,double>& xcurrent, double ftol);

  /**
   * @return Sorted constraints - equations followed by inequalities 
   */
  void insertSort();

  /**
   * @return Output constraint values to the screen 
   */
  void printConstraints();

  /**
   * Takes one argument and returns a ColumnVector
   * @param xcurrent a ColumnVector
   * @return The residual of the constraints.
   */
  virtual Teuchos::SerialDenseVector<int,double> evalResidual(const Teuchos::SerialDenseVector<int,double>& xcurrent ) const ;
  virtual void evalCFGH(const Teuchos::SerialDenseVector<int,double>& xcurrent ) const ;

  /**
   * Takes one argument and returns a real Matrix.
   * @param xcurrent a ColumnVector
   * @return The gradient of the constraints.
   */
  virtual Teuchos::SerialDenseMatrix<int,double> evalGradient(const Teuchos::SerialDenseVector<int,double>& xcurrent ) const ;

  /**
   * Takes two arguments and returns a real SymmetricMatrix
   * @param xcurrent a ColumnVector
   * @param LagMultiplier a ColumnVector
   * @return The Hessian of the constraints multiplied by its assoc. multiplier.
   */
  Teuchos::SerialSymDenseMatrix<int,double> evalHessian(Teuchos::SerialDenseVector<int,double>& xcurrent, 
                            const Teuchos::SerialDenseVector<int,double>& LagMultiplier) const ;
  /**
   * Takes one arguments and returns a Teuchos::SerialSymDenseMatrix<int,double>
   * @param xcurrent a ColumnVector
   * @return The Hessian of the constraints.
   */
  virtual Teuchos::SerialSymDenseMatrix<int,double> evalHessian(Teuchos::SerialDenseVector<int,double>& xcurrent ) const ;

  /**
   * Takes two arguments and returns an array of real SymmetricMatrices
   * @param xcurrent a ColumnVector
   * @param darg an integer argument
   * @return  An array of Hessians.
   */
  virtual OptppArray<Teuchos::SerialSymDenseMatrix<int,double> > evalHessian(Teuchos::SerialDenseVector<int,double>& xcurrent, int darg) const ;


  /**
   * Takes two arguments and returns a bool.
   * @param xcurrent a ColumnVector
   * @param epsilon a real argument
   * @return true - constraints are feasible 
   *
   * @return false - constraints are infeasible
   */
  virtual bool amIFeasible(const Teuchos::SerialDenseVector<int,double>& xcurrent, double epsilon) const;

private:
  /**
   * Sorts an array of constraints.  
   * Equations are followed by inequalities.
   */
  void insertSort(const OptppArray<Constraint>& constraints);
};

} // namespace OPTPP
#endif
