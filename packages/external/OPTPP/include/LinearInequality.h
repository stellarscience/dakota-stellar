#ifndef LinearInequality_h
#define LinearInequality_h

#ifndef LinearConstraint_h
#include "LinearConstraint.h"
#endif

/**
 * LinearInequality is a derived class of LinearConstraint. 
 * Linear Inequalities - Ax >= b
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date Last modified 02/2006
 *
 */

namespace OPTPP {

class LinearInequality: public LinearConstraint {
protected:
  /// type of constraints
  Teuchos::SerialDenseVector<int,double> 		ctype_; 	

public:
/**
 * Default Constructor
 */
  LinearInequality();

/**
 * @param A a real Matrix
 * @param rhs ColumnVector
 * @see LinearInequality(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& rhs,
 *                                 const bool rowFlag)
 * @see LinearInequality(const Teuchos::SerialDenseMatrix<int,double>& A, 
 *   const Teuchos::SerialDenseVector<int,double>& lower, const Teuchos::SerialDenseVector<int,double>& upper)
 */
  LinearInequality(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& rhs);
/**
 * @param A a real Matrix
 * @param rhs ColumnVector
 * @param rowFlag a bool
 * @see LinearInequality(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& rhs);
 * @see LinearInequality(const Teuchos::SerialDenseMatrix<int,double>& A, 
 *   const Teuchos::SerialDenseVector<int,double>& lower, const Teuchos::SerialDenseVector<int,double>& upper)
 */
  LinearInequality(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& rhs, 
                   const bool rowFlag);
/**
 * @param A a real Matrix
 * @param lower a ColumnVector
 * @param upper a ColumnVector
 * @see LinearInequality(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& rhs);
 * @see LinearInequality(const Teuchos::SerialDenseMatrix<int,double>& A, 
 *   const Teuchos::SerialDenseVector<int,double>& rhs, const bool rowFlag)
 */
  LinearInequality(const Teuchos::SerialDenseMatrix<int,double>& A, const Teuchos::SerialDenseVector<int,double>& lower, 
                   const Teuchos::SerialDenseVector<int,double>& upper);

  /**
   * Destructor
   */
  virtual ~LinearInequality(){}

  /**
   * @return Standard form representation of the constraints.
   */
  bool getStdForm() const {return stdForm_;}

  /**
   * @return Lineq 
   */
  virtual Teuchos::SerialDenseVector<int,double> getConstraintType() const {return ctype_;}

  /**
   * Takes one argument and returns a ColumnVector.
   * @param xc a ColumnVector
   * @return Matrix-vector product of A and xc.
   */
  virtual Teuchos::SerialDenseVector<int,double> evalAx(const Teuchos::SerialDenseVector<int,double>& xc) const;

  /**
   * Takes one argument and returns a ColumnVector.
   * @param xc a ColumnVector
   * @return The residual of the linear inequality constraints 
   * evaluated at xc.
   */
  virtual Teuchos::SerialDenseVector<int,double> evalResidual(const Teuchos::SerialDenseVector<int,double>& xc) const;
  virtual void evalCFGH(const Teuchos::SerialDenseVector<int,double>& xc) const;

  /**
   * Takes one argument and returns a real Matrix.
   * @param xc a ColumnVector
   * @return The gradient of the linear inequality constraints
   * evaluated at xc.
   */
  virtual Teuchos::SerialDenseMatrix<int,double> evalGradient(const Teuchos::SerialDenseVector<int,double>& xc) const;

  /**
   * Takes two arguments and returns a bool.
   * @param xc a ColumnVector
   * @param epsilon a real argument.
   * @return The feasibility of the constraints at xc.
   */
  virtual bool amIFeasible(const Teuchos::SerialDenseVector<int,double>& xc, double epsilon) const;
};

} // namespace OPTPP
#endif

