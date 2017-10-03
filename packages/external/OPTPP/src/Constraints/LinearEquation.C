#include "LinearEquation.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003 
//------------------------------------------------------------------------

namespace OPTPP {

// Constructors 
LinearEquation::LinearEquation():
      LinearConstraint(), b_(0), ctype_(0){;}

LinearEquation::LinearEquation(const SerialDenseMatrix<int,double>& A, const SerialDenseVector<int,double>& rhs):
      LinearConstraint(A, rhs), b_(rhs), ctype_(A.numRows())
      {ctype_.resize(numOfCons_); ctype_ = Leqn;}

// Functions for computing various quantities 
SerialDenseVector<int,double> LinearEquation::evalAx(const SerialDenseVector<int,double>& xc) const 
{ 
      int i, index;
      SerialDenseVector<int,double> Ax(numOfCons_);
      SerialDenseMatrix<int,double> tmp(numOfCons_, numOfVars_);
     //  for( i = 1; i <= numOfCons_; i++){
//          index = constraintMappingIndices_[i-1];
// 	 tmp.Row(i) = A_.Row(index);
//       }
//       Ax = tmp*xc;
      for(i=0; i<numOfCons_; i++)
	{index = constraintMappingIndices_[i];
	  for(int j=0; j<numOfVars_; j++)
	    {tmp(i,j) = A_(index,j);}
	}
      Ax.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS, 1.0, tmp, xc, 0.0);
      return Ax;
}

void LinearEquation::evalCFGH(const SerialDenseVector<int,double> & xc) const
{
  return;
}

SerialDenseVector<int,double> LinearEquation::evalResidual(const SerialDenseVector<int,double>& xc) const 
{ 
  
      int i, index;
      // cvalue_         = A_*xc;
      cvalue_.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS, 1.0, A_, xc, 0.0);
      SerialDenseVector<int,double> Ax(evalAx(xc));
      //      SerialDenseVector<int,double> Ax(evalAx(xc).length());
      //	Ax = evalAx(xc);
      SerialDenseVector<int,double> residual(numOfCons_);

      for( i = 0; i < numOfCons_; i++){
         index = constraintMappingIndices_[i];
         residual(i) = Ax(i) - b_(index);
      }
      return residual;
}

SerialDenseMatrix<int,double> LinearEquation::evalGradient(const SerialDenseVector<int,double>& xc) const 
{ 
      int i, index;
      SerialDenseMatrix<int,double> tmp(numOfCons_, numOfVars_);
      SerialDenseMatrix<int,double> tmptrans(numOfVars_, numOfCons_);

//       for( i = 1; i <= numOfCons_; i++){
//          index = constraintMappingIndices_[i-1];
//          tmp.Row(i) = A_.Row(index);
//       }
//       return tmp.t();
      for(i=0;i<numOfCons_;i++)
	{index = constraintMappingIndices_[i];
	for(int j=0;j<numOfVars_;j++)
	  {tmp(i,j) = A_(index,j);}
	}
      for(i=0; i<numOfVars_; i++)
	for(int j=0; j<numOfCons_; j++)
	  {tmptrans(i,j) = tmp(j,i);}
      return tmptrans;

}


bool LinearEquation::amIFeasible(const SerialDenseVector<int,double> & xc, double epsilon) const
{
     int i;
     bool feasible = true;
     SerialDenseVector<int,double> residual(evalResidual(xc));
     //     SerialDenseVector<int,double> residual(evalResidual(xc).length());
     //       residual = evalResidual(xc);
     for(i = 0; i < numOfCons_; i++){
        if( (residual(i) > epsilon) || (residual(i) < -epsilon) ){
           feasible = false;
           break;
        }
     }
     return feasible;
}

} // namespace OPTPP
