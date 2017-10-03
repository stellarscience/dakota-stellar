//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003 
//------------------------------------------------------------------------

#include "LinearConstraint.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

namespace OPTPP {

// Constructors
LinearConstraint::LinearConstraint():
    numOfCons_(0), numOfVars_(0), nnzl_(0), nnzu_(0),
    A_(0,0), Ax_(0), lower_(0), upper_(0), cvalue_(0), 
    cviolation_(0), constraintMappingIndices_(0), stdForm_(true) {;}

LinearConstraint::LinearConstraint(const SerialDenseMatrix<int,double>& A):
    numOfCons_( A.numRows() ), numOfVars_( A.numCols() ), nnzl_(0), nnzu_(0),
    A_(A), Ax_( A.numRows() ), lower_( A.numRows() ), upper_( A.numRows() ),
    cvalue_( A.numRows() ), cviolation_(0), 
    constraintMappingIndices_(0), stdForm_(true)
    { 
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      lower_ = 0.0; upper_ = MAX_BND;
      for(int i = 0; i < numOfCons_; i++){
 	 constraintMappingIndices_.append(i);
         nnzl_++;
      }
      numOfCons_ = nnzl_;
    }

LinearConstraint::LinearConstraint(const SerialDenseMatrix<int,double>& A, const SerialDenseVector<int,double>& b):
    numOfCons_( A.numRows() ), numOfVars_( A.numCols() ), nnzl_(0), nnzu_(0),
    A_(A), Ax_( A.numRows() ), lower_( b ), upper_( b ), 
    cvalue_( A.numRows() ), cviolation_( A.numRows() ),
    constraintMappingIndices_(0), stdForm_(true)
    {
       cvalue_ = 1.0e30; cviolation_ = 0.0;
       for(int i = 0; i < numOfCons_; i++){
          if(lower_(i) > -BIG_BND ){
             constraintMappingIndices_.append(i);
             nnzl_++;
         }
       }
       numOfCons_ = nnzl_;
     }

LinearConstraint::LinearConstraint(const SerialDenseMatrix<int,double>& A, const SerialDenseVector<int,double>& b,
                                   const bool rowFlag):
    numOfCons_( A.numRows() ), numOfVars_( A.numCols() ), nnzl_(0), nnzu_(0),
    A_(A), Ax_( A.numRows() ), lower_( A.numRows() ), upper_( A.numRows() ), 
    cvalue_( A.numRows()), cviolation_( A.numRows() ),
    constraintMappingIndices_(0), stdForm_(rowFlag)
    {
      int i;

      cvalue_  = 1.0e30; cviolation_ = 0.0;
      if( stdForm_ ){
        lower_ = b;
        upper_ = MAX_BND;
        for(i = 0; i < numOfCons_; i++){
            if (lower_(i) > -BIG_BND){
	        constraintMappingIndices_.append(i);
	        nnzl_++;
            }
	}
      }
      else{
        upper_ = b;
	lower_ = MIN_BND; 
        for(i = 0; i < numOfCons_; i++){
            if (upper_(i) < BIG_BND){
	        constraintMappingIndices_.append(i);
	        nnzu_++;
            }
	}
      }

      numOfCons_ = nnzl_ + nnzu_;
   }


LinearConstraint::LinearConstraint(const SerialDenseMatrix<int,double>& A, const SerialDenseVector<int,double>& lower,
                                   const SerialDenseVector<int,double>& upper):
    numOfCons_( 2*A.numRows() ), numOfVars_( A.numCols() ), nnzl_(0), nnzu_(0),
    A_(A), Ax_( A.numRows() ), lower_( lower ), upper_( upper ),
    cvalue_( A.numRows()), cviolation_( A.numRows()) ,
    constraintMappingIndices_(0), stdForm_(true)
    {
       int i, numconstraints = A.numRows();

       cvalue_  = 1.0e30; cviolation_ = 0.0;
       for(i = 0; i < numconstraints; i++){
           if(lower_(i) > -BIG_BND){
	      constraintMappingIndices_.append(i);
	      nnzl_++;
           }
	}
        for(i = 0; i < numconstraints; i++){
            if(upper_(i) < BIG_BND){
	       constraintMappingIndices_.append(i);
	       nnzu_++;
            }
	}
        numOfCons_ = nnzl_ + nnzu_;
   }

void LinearConstraint::setA(SerialDenseMatrix<int,double>& A)
{
    if( dimMatch(A) ) 
      A_ = A;
    else 
      OptppmathError("Check matrix dimensions.  Error in the setA method. ");
}

bool LinearConstraint::dimMatch(SerialDenseMatrix<int,double>& A)
{
    bool match = true;
    if (numOfCons_ != A.numRows() || numOfVars_ !=  A.numCols() )
		   match = false;
    return match;
}

SerialSymDenseMatrix<int,double> LinearConstraint::evalHessian(SerialDenseVector<int,double>& xc) const
{
  SerialSymDenseMatrix<int,double> H(numOfVars_);
    H = 0;

    return H;
}

OptppArray<SerialSymDenseMatrix<int,double> > LinearConstraint::evalHessian(SerialDenseVector<int,double>& xc, int darg) const
{
    OptppArray<SerialSymDenseMatrix<int,double> >  H(1);
    SerialSymDenseMatrix<int,double> Htmp(numOfVars_);
    Htmp = 0;
    H[0] = Htmp;

    return H;
}

} // namespace OPTPP
