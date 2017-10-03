#include "NonLinearInequality.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003
// 
// Standard form g(x) >= 0
//------------------------------------------------------------------------

namespace OPTPP {

// Constructors
// PJW: The +/- infinity bound checked is done in the NonLinearConstraint
// constructor.  If any bound is infinite, then the vector ctype_
// must be resized accordingly.
NonLinearInequality::NonLinearInequality():
     NonLinearConstraint(), ctype_(0), oneSided_(true){}

NonLinearInequality::NonLinearInequality(NLP* nlprob, int numconstraints):
     NonLinearConstraint(nlprob,true,numconstraints), 
     ctype_( numconstraints ), oneSided_(true)
     {ctype_.resize(numOfCons_); ctype_ = NLineq;}

NonLinearInequality::NonLinearInequality(NLP* nlprob, const bool flag, 
                         int numconstraints):
     NonLinearConstraint(nlprob,flag,numconstraints), 
     ctype_( numconstraints ), oneSided_(true)
     {ctype_.resize(numOfCons_); ctype_ = NLineq;}

NonLinearInequality::NonLinearInequality(NLP* nlprob, const SerialDenseVector<int,double>& b, 
                         int numconstraints):
     NonLinearConstraint(nlprob,b,true,numconstraints), 
     ctype_( numconstraints ), oneSided_(true)
     {ctype_.resize(numOfCons_); ctype_ = NLineq;}

NonLinearInequality::NonLinearInequality(NLP* nlprob, const SerialDenseVector<int,double>& b, 
                                     const bool flag, int numconstraints):
     NonLinearConstraint(nlprob,b,flag,numconstraints), 
     ctype_( numconstraints ), oneSided_(true)
     {ctype_.resize(numOfCons_); ctype_ = NLineq;}

NonLinearInequality::NonLinearInequality(NLP* nlprob, 
                                     const SerialDenseVector<int,double>& lower, 
                                     const SerialDenseVector<int,double>& upper,
                                     int numconstraints):
     NonLinearConstraint(nlprob, lower, upper,numconstraints), 
     ctype_(2*numconstraints ), oneSided_(false)
     {ctype_.resize(numOfCons_); ctype_ = NLineq;}


// Methods
void NonLinearInequality::evalCFGH(const SerialDenseVector<int,double> & xc) const
{
  nlp_->evalC(xc);
}

SerialDenseVector<int,double> NonLinearInequality::evalResidual(const SerialDenseVector<int,double> & xc) const 
{
  
      int i, index;
      SerialDenseVector<int,double> resid( numOfCons_);
      cvalue_ = nlp_->evalCF(xc);
     
      for( i = 0; i < nnzl_; i++){
          index = constraintMappingIndices_[i];
	  resid(i) = cvalue_(index) - lower_(index); 
      }
      for( i = nnzl_; i < numOfCons_; i++){
          index = constraintMappingIndices_[i];
	  resid(i) = upper_(index) - cvalue_(index); 
      }
      return resid;
}

SerialDenseMatrix<int,double> NonLinearInequality::evalGradient(const SerialDenseVector<int,double> & xc) const 
{
      int j, index;
      SerialDenseMatrix<int,double> grad(numOfVars_, numOfCons_);
      SerialDenseMatrix<int,double> constraint_grad(nlp_->evalCG(xc));
      //      SerialDenseMatrix<int,double> constraint_grad(nlp_->evalCG(xc).numRows(),nlp_->evalCG(xc).numCols());
      //	constraint_grad = nlp_->evalCG(xc);
     
      for( j = 0; j <nnzl_; j++){
          index = constraintMappingIndices_[j];
	  for(int i=0; i<numOfVars_; i++)
	    grad(i,j) = constraint_grad(i,index);
	    //grad.Column(j) = constraint_grad.Column(index);
      }
      for( j = nnzl_; j < numOfCons_; j++){
          index = constraintMappingIndices_[j];
	  for(int i=0; i<numOfVars_; i++)
	    {grad(i,j) = -constraint_grad(i,index);
	      //grad.Column(j) = -constraint_grad.Column(index);
	    }
      }
      return grad;
}

SerialSymDenseMatrix<int,double> NonLinearInequality::evalHessian(SerialDenseVector<int,double> & xc) const 
{
      SerialSymDenseMatrix<int,double> hess, constraint_hess, nconstraint_hess;
      constraint_hess = nlp_->evalCH(xc);
     
      if(oneSided_){
        if(stdForm_) 
          return constraint_hess;
        else
          {constraint_hess *=-1;
	    return constraint_hess;}
      }
      else{
         nconstraint_hess = constraint_hess;
	 nconstraint_hess *= -1;
	 // hess = constraint_hess & nconstraint_hess;
         int alpha = constraint_hess.numRows()+nconstraint_hess.numRows();
	 int beta = constraint_hess.numCols();
	 for(int i=0; i<alpha; i++)
	   for(int j=0; j<beta; j++)
	     {if(i<constraint_hess.numRows())
		 {hess(i,j) = constraint_hess(i,j);}
	       else
		 {hess(i,j) = nconstraint_hess(i,j);}
	     }
	 return hess;
      }
}

OptppArray<SerialSymDenseMatrix<int,double> > NonLinearInequality::evalHessian(SerialDenseVector<int,double>& xc,
                                                        int darg)const 
{
      int i, index;
      OptppArray<SerialSymDenseMatrix<int,double> > hess(numOfCons_,numOfCons_);
      OptppArray<SerialSymDenseMatrix<int,double> > constraint_hess = nlp_->evalCH(xc,darg);

      for( i = 0; i < nnzl_; i++){
          index = constraintMappingIndices_[i];
	  hess[i] = constraint_hess[index];
      }
      for( i = nnzl_; i < numOfCons_; i++){
          index = constraintMappingIndices_[i];
	  constraint_hess[index] *= -1;
	  hess[i] = constraint_hess[index];
      }
      return hess;
}


bool NonLinearInequality::amIFeasible(const SerialDenseVector<int,double> & xc, double epsilon) const
{
      int i, index;
      bool feasible = true;
      SerialDenseVector<int,double> residual = evalResidual(xc);
      for(i = 0; i < numOfCons_; i++){
         index = constraintMappingIndices_[i];
         if(residual(i) < -epsilon ){
            cviolation_(index) = residual(i);
            feasible = false;
         }
      }
      return feasible;
}

} // namespace OPTPP
