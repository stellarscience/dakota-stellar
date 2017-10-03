#ifndef LinearEquation_h
#define LinearEquation_h


#include "LinearConstraint.h"

/**
 * LinearEquation is a derived class of LinearConstraint.
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date  modified 02/2006
 *
 */

//----------------------------------------------------------------------
// Linear Equations
//----------------------------------------------------------------------

namespace OPTPP {

class LinearEquation: public LinearConstraint {

protected:
  /// Right-hand side of equation
  Teuchos::SerialDenseVector<int,double> b_;		
  /// Type of constraint - Leqn
  Teuchos::SerialDenseVector<int,double> ctype_;		

public:

/**
 * Default Constructor
 * @see LinearEquation(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& rhs);
 */
  LinearEquation();
/**
 * @param A a real Teuchos::SerialDenseMatrix<int,double>
 * @param rhs Teuchos::SerialDenseVector<int,double>
 * @see LinearEquation()
 */
  LinearEquation(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& rhs);

/**
 * Destructor
 */
  virtual ~LinearEquation(){}

/**
 * @return Type of constraint - Leqn
 */
  virtual Teuchos::SerialDenseVector<int,double> getConstraintType() const { return ctype_;};

/**
 * @return The right-hand side of the equation.
 */
  Teuchos::SerialDenseVector<int,double> getB() const {return b_;}

/**
 * Takes one argument and returns a ColumnVector.
 * @param xc a ColumnVector
 * @return Matrix-vector product of A and xc.
 */
  virtual Teuchos::SerialDenseVector<int,double> evalAx(const Teuchos::SerialDenseVector<int,double>& xc) const;

  /**
   * Takes one argument and returns a ColumnVector.
   * @param xc a ColumnVector
   * @return The residual of the linear equations
   * evaluated at xc.
   */
  virtual Teuchos::SerialDenseVector<int,double> evalResidual(const Teuchos::SerialDenseVector<int,double>& xc) const;
  virtual void evalCFGH(const Teuchos::SerialDenseVector<int,double>& xc) const;

  /**
   * Takes one argument and returns a real Matrix.
   * @param xc a ColumnVector
   * @return The gradient of the linear equations
   * evaluated at xc.
   */
  virtual Teuchos::SerialDenseMatrix<int,double> evalGradient(const Teuchos::SerialDenseVector<int,double>& xc) const;

  /**
   * Takes two arguments and returns a bool.
   * @param xc a ColumnVector
   * @param epsilon a real argument.
   * @return The feasibility of the linear equations at xc.
   */
  virtual bool amIFeasible(const Teuchos::SerialDenseVector<int,double>& xc, double epsilon) const;
};

} // namespace OPTPP
#endif
