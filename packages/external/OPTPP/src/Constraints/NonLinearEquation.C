#include "NonLinearEquation.h"

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
NonLinearEquation::NonLinearEquation():
   NonLinearConstraint(), b_(0), ctype_(0){} 

NonLinearEquation::NonLinearEquation(NLP* nlprob, int numconstraints ):
   NonLinearConstraint(nlprob,numconstraints), b_(numconstraints), 
   ctype_(numconstraints) 
   {b_ = 0.0; ctype_.resize(numOfCons_); ctype_ = NLeqn; } 

NonLinearEquation::NonLinearEquation(NLP* nlprob, const SerialDenseVector<int,double>& b,
                                 int numconstraints ):
   NonLinearConstraint(nlprob, b, numconstraints), b_(b),
   ctype_(nlprob->getDim()) 
   {ctype_.resize(numOfCons_); ctype_ = NLeqn;} 

// Functions for computing various quantities 
void NonLinearEquation::evalCFGH(const SerialDenseVector<int,double> & xc) const
{
  nlp_->evalC(xc);
}

SerialDenseVector<int,double> NonLinearEquation::evalResidual(const SerialDenseVector<int,double>& xc) const 
{ 
  
      int i, index;
      SerialDenseVector<int,double> resid( numOfCons_);
      cvalue_ = nlp_->evalCF(xc);
     
      // 08/06/2001 PJW - Okay to only loop over nnzl.
      // We assume nnzl = total constraints and nnzu = 0
      // for the equality constrained case.
      for( i = 0; i < nnzl_; i++){
          index = constraintMappingIndices_[i];
	  resid(i) = cvalue_(index) - b_(index); 
      }
      return resid;
}

SerialDenseMatrix<int,double> NonLinearEquation::evalGradient(const SerialDenseVector<int,double>& xc) const 
{ 
      int j, index;
      SerialDenseMatrix<int,double> grad(numOfVars_, numOfCons_);
      SerialDenseMatrix<int,double> constraint_grad(nlp_->evalCG(xc));
      //      SerialDenseMatrix<int,double> constraint_grad(nlp_->evalCG(xc).numRows(),nlp_->evalCG(xc).numCols());
      //	constraint_grad = nlp_->evalCG(xc);
     
      for( j = 0; j < nnzl_; j++){
          index = constraintMappingIndices_[j];
	  for(int i=0;i<numOfVars_;i++)
	    grad(i,j) = constraint_grad(i,index);
	    // grad.Column(j) = constraint_grad.Column(index);
      }
      return grad;
}

SerialSymDenseMatrix<int,double> NonLinearEquation::evalHessian(SerialDenseVector<int,double>& xc) const 
{ 
  SerialSymDenseMatrix<int,double> hess(numOfCons_), constraint_hess(nlp_->evalCH(xc));
  //  SerialSymDenseMatrix<int,double> hess(nlp_->evalCH(xc).numRows()), constraint_hess(nlp_->evalCH(xc).numRows());
  //     constraint_hess = nlp_->evalCH(xc);
     
     hess = constraint_hess;

     return hess;
}

OptppArray<SerialSymDenseMatrix<int,double> > NonLinearEquation::evalHessian(SerialDenseVector<int,double>& xc,
							  int darg) const 
{ 
     int i, index;
     OptppArray<SerialSymDenseMatrix<int,double> > hess(numOfCons_);
     OptppArray<SerialSymDenseMatrix<int,double> > constraint_hess = nlp_->evalCH(xc,darg);
     
     for( i = 0; i < nnzl_; i++){
         index = constraintMappingIndices_[i];
	 
	 hess[i] = constraint_hess[index];
     }
      
     return hess;
}

bool NonLinearEquation::amIFeasible(const SerialDenseVector<int,double>& xc, double epsilon) const
{
     int i;
     bool feasible = true;
     SerialDenseVector<int,double> residual(evalResidual(xc));
     //     SerialDenseVector<int,double> residual(evalResidual(xc).length());
     //       residual = evalResidual(xc);
     for(i = 0; i < numOfCons_; i++){
        if( (residual(i) < -epsilon) || (residual(i) > epsilon) ){
           feasible = false;
           break;
        }
     }
     return feasible;
}

} // namespace OPTPP
