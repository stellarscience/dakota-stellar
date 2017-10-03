//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003
//------------------------------------------------------------------------

#include "Constraint.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

namespace OPTPP {

// Constructors 
Constraint::Constraint(): ptr_(0){;}
Constraint::Constraint(ConstraintBase* base): ptr_(base){;}

int Constraint::getNumOfCons() const
{
  int result = ptr_->getNumOfCons();
  return result;
}

int Constraint::getNumOfVars() const
{
  int result = ptr_->getNumOfVars();
  return result;
}

SerialDenseVector<int,double> Constraint::getLower() const
{
  SerialDenseVector<int,double> result(ptr_->getLower().length());
    result = ptr_->getLower();
  return result;
}

SerialDenseVector<int,double> Constraint::getUpper() const
{
  SerialDenseVector<int,double> result(ptr_->getUpper().length());
    result = ptr_->getUpper();
  return result;
}

SerialDenseVector<int,double> Constraint::getConstraintType() const
{
  SerialDenseVector<int,double> result(ptr_->getConstraintType().length());
    result = ptr_->getConstraintType();
  return result;
}

SerialDenseVector<int,double> Constraint::getConstraintValue() const
{
  SerialDenseVector<int,double> result(ptr_->getConstraintValue().length());
    result = ptr_->getConstraintValue();
  return result;
}

SerialDenseVector<int,double> Constraint::getConstraintViolation() const
{
  SerialDenseVector<int,double> result(ptr_->getConstraintViolation().length());
    result = ptr_->getConstraintViolation();
  return result;
}

OptppArray<int> Constraint::getConstraintMappingIndices() const
{
  OptppArray<int> result = ptr_->getConstraintMappingIndices();
  return result;
}

void Constraint::evalCFGH(const SerialDenseVector<int,double> & xcurrent) const
{
   ptr_->evalCFGH(xcurrent);
}

SerialDenseVector<int,double> Constraint::evalResidual(const SerialDenseVector<int,double>& xcurrent) const 
{
 
  SerialDenseVector<int,double> result(ptr_->evalResidual(xcurrent));
  //  SerialDenseVector<int,double> result(ptr_->evalResidual(xcurrent).length());
  //   result = ptr_->evalResidual(xcurrent);
   
   return result;
}

SerialDenseMatrix<int,double> Constraint::evalGradient(const SerialDenseVector<int,double>& xcurrent) const 
{
  SerialDenseMatrix<int,double> result(ptr_->evalGradient(xcurrent));
  //  SerialDenseMatrix<int,double> result(ptr_->evalGradient(xcurrent).numRows(),ptr_->evalGradient(xcurrent).numCols());
  //   result = ptr_->evalGradient(xcurrent);
   return result;
}

SerialSymDenseMatrix<int,double> Constraint::evalHessian(SerialDenseVector<int,double>& xcurrent) const 
{
  SerialSymDenseMatrix<int,double> result(ptr_->evalHessian(xcurrent));
  //  SerialSymDenseMatrix<int,double> result(ptr_->evalHessian(xcurrent).numRows());
  //   result = ptr_->evalHessian(xcurrent);
 
   return result;
}

OptppArray<SerialSymDenseMatrix<int,double> > Constraint::evalHessian(SerialDenseVector<int,double>& xcurrent, int darg) const 
{
   OptppArray<SerialSymDenseMatrix<int,double> > result;
   result = ptr_->evalHessian(xcurrent,darg);

   return result;
}

bool Constraint::amIFeasible(const SerialDenseVector<int,double>& xcurrent,double epsilon) const 
{
   bool result;
   result =  ptr_->amIFeasible(xcurrent,epsilon);
   return result;
}

} // namespace OPTPP
