
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include "NonLinearConstraint.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 10/19/1999
//------------------------------------------------------------------------

namespace OPTPP {

// Constructors
NonLinearConstraint::NonLinearConstraint():
  nlp_(0), lower_(0), upper_(0), cvalue_(0), cviolation_(0),
  numOfCons_(0), numOfVars_(0), nnzl_(0), nnzu_(0), 
  constraintMappingIndices_(0), stdForm_(true){;}

NonLinearConstraint::NonLinearConstraint(NLP* nlprob, int numconstraints):
   nlp_(nlprob), lower_(numconstraints), upper_(numconstraints), 
   cvalue_(numconstraints), cviolation_(numconstraints),
   numOfCons_(numconstraints), numOfVars_(nlprob->getDim()), 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), stdForm_(true)
   { 
     cvalue_ = 1.0e30; cviolation_ = 0.0;
     lower_ = 0.0; upper_ = MAX_BND;
     nnzl_ = numconstraints;
     for(int i = 0; i < numconstraints; i++)
        constraintMappingIndices_.append(i);
   }

NonLinearConstraint::NonLinearConstraint(NLP* nlprob, const bool conFlag, int numconstraints):
   nlp_(nlprob), lower_(numconstraints), upper_(numconstraints),
   cvalue_(numconstraints), cviolation_(numconstraints),
   numOfCons_(numconstraints), numOfVars_(nlprob->getDim()), 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), stdForm_(conFlag)
   {
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      if (stdForm_) { 
         lower_ = 0.0; upper_ = MAX_BND; 
         nnzl_ = numconstraints;
         for(int i = 0; i < numconstraints; i++)
            constraintMappingIndices_.append(i);
      }
      else{ 
	 lower_ = MIN_BND; upper_ = 0.0; 
         nnzu_ = numconstraints;
         for(int i = 0; i < numconstraints; i++)
            constraintMappingIndices_.append(i);
      }
   }

NonLinearConstraint::NonLinearConstraint(NLP* nlprob, const SerialDenseVector<int,double>& rhs, 
                                                 int numconstraints):
   nlp_(nlprob), lower_(rhs), upper_(rhs), 
   cvalue_(numconstraints), cviolation_(numconstraints), 
   numOfCons_(numconstraints), numOfVars_(nlprob->getDim()), 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), stdForm_(true) 
   {
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      for(int i = 0; i < numconstraints; i++){
         if(lower_(i) > -BIG_BND && upper_(i) < BIG_BND ){
            constraintMappingIndices_.append(i);
	    nnzl_++;
         }
      }
      numOfCons_ = nnzl_;
   }


NonLinearConstraint::NonLinearConstraint(NLP* nlprob, const SerialDenseVector<int,double>& rhs,
                                           const bool conFlag, int numconstraints):
   nlp_(nlprob), lower_(numconstraints), upper_(numconstraints), 
   cvalue_(numconstraints), cviolation_(numconstraints),
   numOfCons_(numconstraints), numOfVars_(nlprob->getDim()), 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), stdForm_(conFlag)
   {
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      if (stdForm_) { 
	 lower_ = rhs;
         upper_ = MAX_BND; 
         for(int i = 0; i < numconstraints; i++){
	    if( lower_(i) > -BIG_BND ){
               constraintMappingIndices_.append(i);
               nnzl_++;
            }
         }
      }
      else{ 
	 lower_ = MIN_BND; 
	 upper_ = rhs; 
         for(int i = 0; i < numconstraints; i++){
	    if( upper_(i) < BIG_BND ){
              constraintMappingIndices_.append(i);
              nnzu_++;
            }
         }
      }
     numOfCons_ = nnzl_ + nnzu_;
   }


NonLinearConstraint::NonLinearConstraint(NLP* nlprob, const SerialDenseVector<int,double>& lower,
                                       const SerialDenseVector<int,double>& upper, int numconstraints):
   nlp_(nlprob), lower_(lower), upper_(upper), 
   cvalue_(numconstraints), cviolation_(numconstraints),
   numOfCons_(2*numconstraints), numOfVars_(nlprob->getDim()), 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), stdForm_(true)
   { 
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      for(int i = 0; i < numconstraints; i++){
	 if( lower_(i) > -BIG_BND ){
           constraintMappingIndices_.append(i);
           nnzl_++;
         }
      }
      for(int i = 0; i < numconstraints; i++){
	  if( upper_(i) < BIG_BND ){
            constraintMappingIndices_.append(i);
            nnzu_++;
          }
      }
      numOfCons_ = nnzl_ + nnzu_;
   }

#ifdef DAKOTA_OPTPP
NonLinearConstraint::NonLinearConstraint(NLP* nlprob, const SerialDenseVector<int,double>& lower,
                                       const SerialDenseVector<int,double>& upper, 
				       int ne, int ni):
   nlp_(nlprob), lower_(lower), upper_(upper), cvalue_(ne+ni), cviolation_(ne+ni),
   numOfCons_(ne+ni), numOfVars_(nlprob->getDim()) , 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), 
   stdForm_(true), ctype_(ne +ni) 
   { 
      OptppArray<int> temp;
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      for(int i = 0; i < ne; i++){
	 if( lower_(i) > -BIG_BND ){
           constraintMappingIndices_.append(i);
	   temp.append(NLeqn);
           nnzl_++;
         }
      }
      for(int i = ne; i < ne + ni; i++){
	 if( lower_(i) > -BIG_BND ){
           constraintMappingIndices_.append(i);
	   temp.append(NLineq);
           nnzl_++;
         }
      }
      for(int i = ne; i < ne + ni; i++){
	 if( upper_(i) < BIG_BND ){
           constraintMappingIndices_.append(i);
	   temp.append(NLineq);
           nnzu_++;
         }
      }
      numOfCons_ = nnzl_ + nnzu_;
      ctype_.resize(numOfCons_);
      for(int i = 0; i < numOfCons_; i++)
         ctype_(i) = temp[i];
   }

void NonLinearConstraint::evalCFGH(const SerialDenseVector<int,double> & xc) const
{
  nlp_->evalC(xc);
}

SerialDenseVector<int,double> NonLinearConstraint::evalResidual(const SerialDenseVector<int,double> & xc) const 
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

SerialDenseMatrix<int,double> NonLinearConstraint::evalGradient(const SerialDenseVector<int,double> & xc) const 
{
      int j, index;
      SerialDenseMatrix<int,double> grad(numOfVars_, numOfCons_);
      SerialDenseMatrix<int,double> constraint_grad(nlp_->evalCG(xc));
      //      SerialDenseMatrix<int,double> constraint_grad(nlp_->evalCG(xc).numRows(),nlp_->evalCG(xc).numCols());
      //	constraint_grad = nlp_->evalCG(xc);
     
      for( j = 0; j < nnzl_; j++){
          index = constraintMappingIndices_[j];
	  //grad.Column(j) = constraint_grad.Column(index);
	  for(int k=0; k< numOfVars_;k++)
	    {grad(k,j) = constraint_grad(k,index);}
      }
      for( j = nnzl_; j < numOfCons_; j++){
          index = constraintMappingIndices_[j];
	  //grad.Column(j) = -constraint_grad.Column(index);
	  for(int k=0; k< numOfVars_;k++)
	    {grad(k,j) = -constraint_grad(k,index);}
      }
      return grad;
}

SerialSymDenseMatrix<int,double> NonLinearConstraint::evalHessian(SerialDenseVector<int,double> & xc) const 
{
   // 09/05/01 PJW Dummy routine 
      SerialSymDenseMatrix<int,double> hess, nconstraint_hess;
      SerialSymDenseMatrix<int,double> constraint_hess(nlp_->evalCH(xc));
      //      constraint_hess.reshape(nlp_->evalCH(xc).numRows());
      //      constraint_hess = nlp_->evalCH(xc);
     
      nconstraint_hess.reshape(constraint_hess.numRows());
      nconstraint_hess = constraint_hess;
      nconstraint_hess *= -1;
     
      //hess = constraint_hess & nconstraint_hess;
      for(int i=0;i<constraint_hess.numRows()+nconstraint_hess.numRows();i++)
	{for(int j=0;j<constraint_hess.numCols();j++)
	    {if(i<constraint_hess.numRows())
		{hess(i,j) = constraint_hess(i,j);}
	      else{hess(i,j) = nconstraint_hess(i,j);}
	    }
	}

      return hess;
}

OptppArray<SerialSymDenseMatrix<int,double> > NonLinearConstraint::evalHessian(SerialDenseVector<int,double>& xc,
                                                        int darg)const 
{
      int i, index;
      OptppArray<SerialSymDenseMatrix<int,double> > hess(numOfCons_);
      OptppArray<SerialSymDenseMatrix<int,double> > constraint_hess = nlp_->evalCH(xc,darg);
       
      for( i = 0; i < nnzl_; i++){
          index = constraintMappingIndices_[i];
	  hess[i] = constraint_hess[index];
      }
      for( i = nnzl_; i < numOfCons_; i++){
          index = constraintMappingIndices_[i];
	  hess[i] = constraint_hess[index];
	  hess[i] *= -1;
      }


      return hess;
}


bool NonLinearConstraint::amIFeasible(const SerialDenseVector<int,double> & xc, double epsilon) const
{
      int i, index;
      bool feasible = true;
      SerialDenseVector<int,double> residual(evalResidual(xc));
      //      SerialDenseVector<int,double> residual(evalResidual(xc).length());
      //	residual = evalResidual(xc);
      for(i = 0; i < numOfCons_; i++){
         index = constraintMappingIndices_[i];
         if(residual(i) < -epsilon ){
            cviolation_(index) = residual(i);
            feasible = false;
         }
      }
      return feasible;
}
#endif // DAKOTA_OPTPP 

} // namespace OPTPP
