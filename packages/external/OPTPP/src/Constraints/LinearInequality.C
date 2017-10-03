//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003 
//------------------------------------------------------------------------

#include "LinearInequality.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

namespace OPTPP {

// Constructors
LinearInequality::LinearInequality():
      LinearConstraint( ), ctype_(0) {}

LinearInequality::LinearInequality(const SerialDenseMatrix<int,double>& A, const SerialDenseVector<int,double>& rhs):
      LinearConstraint(A,rhs,true), ctype_(A.numRows())
      {ctype_.resize(numOfCons_); ctype_ = Lineq;}

LinearInequality::LinearInequality(const SerialDenseMatrix<int,double>& A, const SerialDenseVector<int,double>& rhs, 
                                   const bool rowFlag):
      LinearConstraint(A,rhs,rowFlag), ctype_(A.numRows())
      {ctype_.resize(numOfCons_); ctype_ = Lineq;}

LinearInequality::LinearInequality(const SerialDenseMatrix<int,double>& A, const SerialDenseVector<int,double>& lower, 
                                   const SerialDenseVector<int,double>& upper):
      LinearConstraint(A,lower,upper), ctype_(2*A.numRows())
      {ctype_.resize(numOfCons_); ctype_ = Lineq;}

// Evaluation Methods
SerialDenseVector<int,double> LinearInequality::evalAx(const SerialDenseVector<int,double>& xc) const 
{

      int i, index, nnz = nnzl_ + nnzu_;
      SerialDenseVector<int,double> Ax(numOfCons_);
      SerialDenseMatrix<int,double> tmp(numOfCons_, numOfVars_);
     //  for( i = 1; i <= nnzl_; i++){
//          index = constraintMappingIndices_[i-1];
// 	 tmp.Row(i) = A_.Row(index);
//       }
//       for( i = nnzl_+1; i <= nnz; i++){
//          index = constraintMappingIndices_[i-1];
// 	 tmp.Row(i) = -A_.Row(index);
//       }

      for(i=0;i<nnzl_;i++)
	for(int j=0; j<numOfVars_; j++)
	  {index = constraintMappingIndices_[i];
	    tmp(i,j) = A_(index,j);}
      for(i=nnzl_;i<nnz;i++)
	for(int j=0; j<numOfVars_; j++)
	  {index = constraintMappingIndices_[i];
	    tmp(i,j)=-A_(index,j);}

      //Ax = tmp*xc;
      Ax.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS, 1.0, tmp,xc,0.0); 
      return Ax;
}

void LinearInequality::evalCFGH(const SerialDenseVector<int,double> & xc) const
{
  return;
}

SerialDenseVector<int,double> LinearInequality::evalResidual(const SerialDenseVector<int,double>& xc) const 
{
      int i, index, nnz = nnzl_ + nnzu_;
      // cvalue_               = A_*xc;
      cvalue_.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS, 1.0, A_, xc, 0.0);
      SerialDenseVector<int,double> residual = evalAx(xc);
      for( i = 0; i < nnzl_; i++){
         index = constraintMappingIndices_[i];
         residual(i)-= lower_(index); 
      }
      for( i = nnzl_; i < nnz; i++){
         index = constraintMappingIndices_[i];
         residual(i)+= upper_(index); 
      }
      return residual;
}

SerialDenseMatrix<int,double> LinearInequality::evalGradient(const SerialDenseVector<int,double>& xc) const 
{
      int i, index, nnz = nnzl_ + nnzu_;
      SerialDenseMatrix<int,double> tmp(numOfCons_, numOfVars_); 
      SerialDenseMatrix<int,double> tmptrans(numOfVars_,numOfCons_);

      for(i=0;i<nnzl_;i++)
	for(int j=0; j<numOfVars_; j++)
	  {index = constraintMappingIndices_[i];
	    tmp(i,j) = A_(index,j);}
      for(i=nnzl_;i<nnz;i++)
	for(int j=0; j<numOfVars_; j++)
	  {index = constraintMappingIndices_[i];
	    tmp(i,j)=-A_(index,j);}

//       for( i = 1; i <= nnzl_; i++){
//          index = constraintMappingIndices_[i-1];
//          tmp.Row(i) = A_.Row(index); 
//       }
//       for( i = nnzl_+1; i <= nnz; i++){
//          index = constraintMappingIndices_[i-1];
//          tmp.Row(i) = -A_.Row(index); 
//       }
//       return tmp.t();
      for(i=0; i<numOfVars_; i++)
	for(int j=0; j<numOfCons_; j++)
	  {tmptrans(i,j) = tmp(j,i);}
      return tmptrans;

}

bool LinearInequality::amIFeasible(const SerialDenseVector<int,double>& xc, double epsilon) const
{
     bool feasible = true;
     int i, index;
     SerialDenseVector<int,double> residual(evalResidual(xc));
     //     SerialDenseVector<int,double> residual(evalResidual(xc).length());
     //       residual = evalResidual(xc);

     for(i = 0; i < numOfCons_; i++){
       index = constraintMappingIndices_[i];
       if( residual(i) < -epsilon ){
          cviolation_(index) = residual(i);
          feasible = false;
       }
     }
     return feasible;
}

} // namespace OPTPP
