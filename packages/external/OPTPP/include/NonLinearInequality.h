#ifndef NonLinearInequality_h
#define NonLinearInequality_h

#include "NonLinearConstraint.h"

/**
 * NonLinearInequality is a derived class of NonLinearConstraint.
 * Standard Form  g(x) >= 0 
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date   Last modified 02/2006
 *
 */

namespace OPTPP {

class NonLinearInequality: public NonLinearConstraint{
protected:
 /// Type of constraint - NLineq
 Teuchos::SerialDenseVector<int,double> ctype_;	
 /// Denotes whether we have a 1-sided constraint
 const bool oneSided_;	

public:
/**
 * Default Constructor
 */
  NonLinearInequality();
/**
 * @param nlprob a pointer to an NLP object
 * @param numconstraints an integer argument
 * @note Assumes the right-hand side = 0 and the inequality is in standard form
 */
  NonLinearInequality(NLP* nlprob, int numconstraints = 1);
/**
 * @param nlprob a pointer to an NLP object
 * @param rhs Teuchos::SerialDenseVector<int,double>
 * @param numconstraints an integer argument
 * @note Assumes nonzero right-hand side and the inequality is in standard form
 */
  NonLinearInequality(NLP* nlprob, const Teuchos::SerialDenseVector<int,double>& rhs, 
                   int numconstraints = 1);
/**
 * @param nlprob a pointer to an NLP object
 * @param flag a bool argument
 * @param numconstraints an integer argument
 * @note Assumes the right-hand side = 0 
 */
  NonLinearInequality(NLP* nlprob, const bool flag, int numconstraints = 1);

/**
 * @param nlprob a pointer to an NLP object
 * @param rhs Teuchos::SerialDenseVector<int,double>
 * @param numconstraints an integer argument
 * @param flag a bool argument
 * @note Assumes nonzero right-hand side 
 */
  NonLinearInequality(NLP* nlprob, const Teuchos::SerialDenseVector<int,double>& rhs, const bool flag, 
                   int numconstraints = 1);
/**
 * @param nlprob a pointer to an NLP object
 * @param lower Teuchos::SerialDenseVector<int,double>
 * @param upper Teuchos::SerialDenseVector<int,double>
 * @param numconstraints an integer argument
 */
  NonLinearInequality(NLP* nlprob, const Teuchos::SerialDenseVector<int,double>& lower, 
                   const Teuchos::SerialDenseVector<int,double>& upper, int numconstraints = 1);

/**
 * Destructor
 */
  virtual ~NonLinearInequality(){}

/**
 * @return Type of constraint - NLineq
 */
  Teuchos::SerialDenseVector<int,double> getConstraintType() const {return ctype_;}

  /**
   * Takes one argument and returns a ColumnVector.
   * @param xc a ColumnVector
   * @return The residual of the nonlinear inequalities evaluated at xc.
   */
  Teuchos::SerialDenseVector<int,double> evalResidual(const Teuchos::SerialDenseVector<int,double> & xc) const;
  void evalCFGH(const Teuchos::SerialDenseVector<int,double> & xc) const;

  /**
   * Takes one argument and returns a Matrix.
   * @param xc a ColumnVector
   * @return The gradient of the nonlinear inequalities evaluated at xc.
   */
  Teuchos::SerialDenseMatrix<int,double> evalGradient(const Teuchos::SerialDenseVector<int,double>& xc) const;

  /**
   * Takes one argument and returns a SymmetricMatrix.
   * @param xc a ColumnVector
   * @return The Hessian of the nonlinear inequalities evaluated at xc.
   */
  Teuchos::SerialSymDenseMatrix<int,double> evalHessian(Teuchos::SerialDenseVector<int,double>& xc) const;

  /**
   * Takes two arguments and returns an array of real SymmetricMatrices.
   * @param xc a ColumnVector
   * @param darg an integer argument
   * @return An array of constraint Hessians.
   */
  OptppArray<Teuchos::SerialSymDenseMatrix<int,double> > evalHessian(Teuchos::SerialDenseVector<int,double>& xc, int darg)const;

 /**
  * Takes two arguments and returns a bool.
  * @param xc a ColumnVector
  * @param epsilon a real argument.
  * @return The feasibility of the nonlinear inequalities at xc.
  */
  bool amIFeasible(const Teuchos::SerialDenseVector<int,double>& xc, double epsilon) const;

};

} // namespace OPTPP
#endif

