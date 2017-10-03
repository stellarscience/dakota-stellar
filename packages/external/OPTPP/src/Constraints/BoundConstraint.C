//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003
//------------------------------------------------------------------------

#include "BoundConstraint.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

namespace OPTPP {

// Constructors
BoundConstraint::BoundConstraint(): numOfCons_(0), numOfVars_(0), 
   nnzl_(0), nnzu_(0), lower_(0), upper_(0), cvalue_(0), 
   fixedVar_(0), freeVar_(0), stdForm_(0), 
   ctype_(0), constraintMappingIndices_(0){}

BoundConstraint::BoundConstraint(int nc, const SerialDenseVector<int,double>& lower):
   numOfCons_(0), numOfVars_(nc), nnzl_(0), nnzu_(0), 
   lower_(nc), upper_(nc), cvalue_(nc), fixedVar_(nc,false), 
   freeVar_(nc,true), stdForm_(nc,true), ctype_(nc), constraintMappingIndices_(0)
   { 
     cvalue_ = 1.0e30;
     lower_  = lower; 
     upper_  = MAX_BND;     
     for(int i = 0; i < nc; i++){
        if (lower_(i) > -BIG_BND){
	    constraintMappingIndices_.append(i);
	    nnzl_++;
        }
     }
     numOfCons_ = nnzl_;
     ctype_.resize(nnzl_);
     ctype_ = Bound;
   }

BoundConstraint::BoundConstraint(int nc, const  SerialDenseVector<int,double>& bound, 
                                 const BoolVector& bdFlag):
   numOfCons_(0), numOfVars_(nc), nnzl_(0), nnzu_(0), 
   lower_(nc), upper_(nc), cvalue_(nc), fixedVar_(nc,false), 
   freeVar_(nc,true), stdForm_(nc,bdFlag), ctype_(nc), constraintMappingIndices_(0)
   { 
     cvalue_ = 1.0e30;
     for(int i = 0; i < nc; i++){
        if ( stdForm_(i+1) ){
            lower_(i)  = bound(i); 
	    upper_(i)  = MAX_BND; 
            if (lower_(i) > -BIG_BND){
	        nnzl_++;
	        constraintMappingIndices_.append(i);
            }
	}
     }
     for(int i = 0; i < nc; i++){
        if ( !stdForm_(i+1) ){
	    lower_(i)  = MIN_BND; 
            upper_(i)  = bound(i); 
            if (upper_(i) < BIG_BND){
	        nnzu_++;
	        constraintMappingIndices_.append(i);
            }
	}
     }
     numOfCons_ = nnzl_ + nnzu_;
     ctype_.resize(numOfCons_);
     ctype_ = Bound;
   }

BoundConstraint::BoundConstraint(int nc, const SerialDenseVector<int,double>& lower,
                                 const SerialDenseVector<int,double>& upper):
   numOfCons_(0), numOfVars_(nc), nnzl_(0), nnzu_(0), 
   lower_(nc), upper_(nc), cvalue_(nc),
   fixedVar_(nc,false), freeVar_(nc,true), stdForm_(nc,true), 
   ctype_(2*nc), constraintMappingIndices_(0)
   { 
     cvalue_ = 1.0e30;
     lower_  = lower; 
     for(int i = 0; i < nc; i++){
        if (lower_(i) > -BIG_BND){
	    nnzl_++;
	    constraintMappingIndices_.append(i);
        }
     }
     upper_  = upper; 
     for(int i = 0; i < nc; i++){
        if (upper_(i) <  BIG_BND){
	    nnzu_++;
	    constraintMappingIndices_.append(i);
        }
     }
     numOfCons_ = nnzl_ + nnzu_;
     ctype_.resize(numOfCons_);
     ctype_ = Bound;
     if(!amIConsistent() ) 
       OptppmathError("Error in Constructor - Lower bound exceeds upper bound");  
   }

void BoundConstraint::evalCFGH(const SerialDenseVector<int,double> & xc) const
{
  return;
}

SerialDenseVector<int,double> BoundConstraint::evalResidual(const SerialDenseVector<int,double>& xc) const 
{ 
  
    int i, index, nnz = nnzl_ + nnzu_;
    SerialDenseVector<int,double> resid(nnz);
    for(i = 0; i < nnzl_; i++){
        index = constraintMappingIndices_[i];
	resid(i) = xc(index) - lower_(index);
    }
    for(i = nnzl_; i < nnz; i++){
        index = constraintMappingIndices_[i];
	resid(i) = upper_(index) - xc(index);
    }
    cvalue_ = xc;
    return resid;
}

SerialDenseMatrix<int,double> BoundConstraint::evalGradient(const SerialDenseVector<int,double>& xc) const 
{ 
    int i, j, nnz = nnzl_+ nnzu_;
    SerialDenseMatrix<int,double> D(numOfVars_, nnz);
    D = 0.0;
    
    for(j = 0; j < nnzl_; j++){
      i = constraintMappingIndices_[j];
      D(i,j) = 1.0;
    }
    for(j = nnzl_; j < nnz; j++){
      i = constraintMappingIndices_[j];
      D(i,j) = -1.0;
    }
    return D;
}

SerialSymDenseMatrix<int,double> BoundConstraint::evalHessian(SerialDenseVector<int,double>& xc) const 
{ 

  SerialSymDenseMatrix<int,double> H(numOfCons_);
    H = 0;

    return H;
}

OptppArray<SerialSymDenseMatrix<int,double> > BoundConstraint::evalHessian(SerialDenseVector<int,double>& xc, 
                                           int darg) const 
{
    OptppArray<SerialSymDenseMatrix<int,double> > Hessian(1);
    SerialSymDenseMatrix<int,double> H(numOfCons_);
    H = 0;
    Hessian[0] = H;

    return Hessian;
}

bool BoundConstraint::amIFeasible(const SerialDenseVector<int,double>& xc, double epsilon) const 
{
    int i;
    bool feasible = true;

    for(i = 0; i < numOfVars_; i++)
      if(xc(i) <  lower_(i) || xc(i) > upper_(i)){
        feasible = false;  
        break;
      }
    return feasible;
}

bool BoundConstraint::amIConsistent() const
{
    int i;
    bool consistent = true;
    for(i = 0; i < numOfVars_; i++)
       if(lower_(i) > upper_(i)){
	  consistent = false;
	  break;
      }
    return consistent;
}

} // namespace OPTPP
