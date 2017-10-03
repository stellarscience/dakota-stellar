//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003
//------------------------------------------------------------------------

#include "CompoundConstraint.h"
#include "ioformat.h"
#include <cstring>
#include <string.h>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;
using std::cout;
using std::strcpy;


namespace OPTPP {

// Constructors
CompoundConstraint::CompoundConstraint(): constraints_(0), 
   numOfSets_(0), lower_(0) , upper_(0){}

CompoundConstraint::CompoundConstraint(const Constraint& c1): 
   constraints_(0), numOfSets_(1)
{ 
   constraints_.append(c1); 
   lower_ = getLower(); 
   upper_ = getUpper(); 
}

CompoundConstraint::CompoundConstraint(const Constraint& c1, 
                                       const Constraint& c2): 
   constraints_(0), numOfSets_(2)
{ 
   constraints_.append(c1);
   constraints_.append(c2);
   insertSort();
   lower_       = getLower();
   upper_       = getUpper();
}

CompoundConstraint::CompoundConstraint(const OptppArray<Constraint>& constraints): 
   constraints_(constraints), numOfSets_(constraints.length())
{  
   insertSort();
   lower_  = getLower();  
   upper_   = getUpper();  
}

CompoundConstraint::CompoundConstraint(const CompoundConstraint& cc):
   constraints_(0), numOfSets_(cc.numOfSets_), lower_(cc.lower_), 
   upper_(cc.upper_) 
{ 
   if(numOfSets_ > 0){
     for(int i = 0; i < numOfSets_; i++)
        constraints_.append(cc[i]);
   }
}

CompoundConstraint& CompoundConstraint::operator=(const CompoundConstraint& cc)
{ 
   if(this != &cc){
      numOfSets_ = cc.numOfSets_;
      lower_     = cc.lower_;
      upper_     = cc.upper_;
      for(int i = 0; i < numOfSets_; i++)
         constraints_.append(cc[i]);
   }                                
   return *this;                       
}

// Accessor Methods 
int CompoundConstraint::getNumOfNLCons() const
{
   int Mi, i, k = 0;
   Constraint test;

   for(i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     SerialDenseVector<int,double> temp(test.getConstraintType().length());
       temp =  test.getConstraintType();
     if (temp(0) == NLeqn || temp(0) == NLineq){
        Mi   = test.getNumOfCons();
        k   += Mi; 
     }
   }
   return k;
}

int CompoundConstraint::getNumOfCons() const
{
   int Mi, i, k = 0;
   Constraint test;

   for(i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     Mi   = test.getNumOfCons();
     k   += Mi; 
   }
   return k;
}

int CompoundConstraint::getNumOfVars() const
{
   int Mi, i, k = 0;
   Constraint test;

   for(i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     Mi   = test.getNumOfVars();
     k   += Mi; 
   }
   if( k!= 0 && k == Mi*numOfSets_ ) 
     return Mi;
   else
     return 0;
}

SerialDenseVector<int,double> CompoundConstraint::getLower() const
{
   SerialDenseVector<int,double> result(getNumOfCons());
   Constraint test;

   int alpha = 0;
   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     int beta = test.getLower().length();
     SerialDenseVector<int,double> temp(beta);
     result.resize(alpha+beta);
     temp =  test.getLower();
     for(int j=alpha; j<alpha+beta; j++)
       {result(j) = temp(j-alpha);}
     alpha = alpha+beta;
   }
   return result;
}

SerialDenseVector<int,double> CompoundConstraint::getUpper() const
{
   SerialDenseVector<int,double> result(getNumOfCons());
   Constraint test;
   int alpha = 0;
   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     int beta = test.getUpper().length();
     SerialDenseVector<int,double> temp(beta);
     result.resize(alpha+beta);
     temp =  test.getUpper();
     for(int j=alpha; j<alpha+beta; j++)
       {result(j) = temp(j-alpha);}
     alpha = alpha+beta;
   }
   return result;

}

SerialDenseVector<int,double> CompoundConstraint::getConstraintType() const
{
   SerialDenseVector<int,double> result(getNumOfCons());
   Constraint test;
   int alpha = 0;
   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     int beta = test.getConstraintType().length();
     SerialDenseVector<int,double> temp(beta);
     result.resize(alpha+beta);
     temp =  test.getConstraintType();
     for(int j=alpha; j<alpha+beta; j++)
       {result(j) = temp(j-alpha);}
     alpha = alpha+beta;
   }
   return result;
}

OptppArray<int> CompoundConstraint::getConstraintMappingIndices() const
{
   OptppArray<int> result;
   Constraint test;

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     OptppArray<int> temp =  test.getConstraintMappingIndices();
     for(int j = 0; j < temp.length(); j++)
        result.append(temp[j]);
   }
   return result;

}

SerialDenseVector<int,double> CompoundConstraint::getConstraintValue() const 
{
   /* 
    * Returns the raw value of all the constraints. 
    */

   Constraint test;
   SerialDenseVector<int,double> temp, value(1);
   int alpha = 0;
   int beta;
   for(int i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     beta = test.getConstraintValue().length();
     temp.resize(beta);
     temp   = test.getConstraintValue();
     value.resize(alpha+beta);
     for(int j=alpha; j<alpha+beta; j++)
       {value(j) = temp(j-alpha);}
     alpha = alpha+beta;

   }
   return value;
}

SerialDenseVector<int,double> CompoundConstraint::getNLConstraintValue() const 
{
   /* 
    * Returns the raw value of the nonlinear equations and inequalities.  
    */

   int j = 0;
   Constraint test;
   SerialDenseVector<int,double> temp, type, value(1), zero(1);

   // If there are no nonlinear constraints, initialize a return value of zero
   zero = 0;
   int alpha = 0;
   int beta;
   for(int i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     type.resize(test.getConstraintType().length());
     type   = test.getConstraintType();

     if(type(0) == NLeqn || type(0) == NLineq){
       beta = test.getConstraintValue().length();
       temp.resize(beta);
        temp   = test.getConstraintValue();
	 value.resize(alpha+beta);
	 for(int k=alpha; k<alpha+beta; k++)
	   {value(k) = temp(k-alpha);}
	 alpha = alpha+beta;

        j++;
     }
   }

   if( j == 0) value = 0.0;

   return value;
}

SerialDenseVector<int,double> CompoundConstraint::getConstraintViolation() const 
{
   /* CPJW - Rethink! 
    * Returns the infeasibility with respect to the constraints.  
    */

   Constraint test;
   SerialDenseVector<int,double> temp, value(1);
   int alpha = 0;
   int beta;
   for(int i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     beta = test.getConstraintViolation().length();
     temp.resize(beta);
     temp   = test.getConstraintViolation();
     value.resize(alpha+beta);
     for(int j=alpha; j<alpha+beta; j++)
       {value(j) = temp(j-alpha);}
     alpha = alpha+beta;
    
   }
   return value;
}

void CompoundConstraint::computeDistanceToBounds(SerialDenseVector<int,double>&xc, SerialDenseVector<int,double>& d_lower, SerialDenseVector<int,double>& d_upper)
{
   /* 
    * Returns a point that is feasible with respect to the bound constraints 
    */

   int i, j, nvars;
   Constraint test;
   SerialDenseVector<int,double> type, lower, upper;

   for(i = 0; i < numOfSets_; i++) {
     test   = constraints_[i];
     type.resize(test.getConstraintType().length());
     type   = test.getConstraintType();

     if(type(0) == Bound){
       nvars = test.getNumOfVars();
       lower.resize(test.getLower().length());
       lower = test.getLower();
       upper.resize(test.getUpper().length());
       upper = test.getUpper();

       for(j = 0; j < nvars; j++){
         d_lower(j) = xc(j) - lower(j);
	 d_upper(j) = upper(j) - xc(j);
       }
     }
   }
}

void CompoundConstraint::computeFeasibleBounds(SerialDenseVector<int,double>&xc, double epsilon) 
{
   /* 
    * Returns a point that is feasible with respect to the bound constraints 
    */

   int i, j, nvars;
   Constraint test;
   SerialDenseVector<int,double> type, lower, upper;

   for(i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     type.resize(test.getConstraintType().length());
     type   = test.getConstraintType();

     if(type(0) == Bound){
       nvars = test.getNumOfVars();
       lower.resize(test.getLower().length());
       lower = test.getLower();
       upper.resize(test.getUpper().length());
       upper = test.getUpper();

       for(j = 0; j < nvars-1; j++){
         if( xc(j) < lower(j) || xc(j) > upper(j)){
           if(lower(j) > -BIG_BND && upper(j) == MAX_BND)
              xc(j) = lower(j) + epsilon;
           else if(upper(j) < BIG_BND && lower(j) == MIN_BND)
              xc(j) = upper(j) + epsilon;
           else
              xc(j) = (lower(j) + upper(j))/2.0 + epsilon;
         }
       }
     }
   }
}

void CompoundConstraint::computeFeasibleInequalities(SerialDenseVector<int,double>&xc, double ftol) 
{
   /* 
    * Returns a point that is feasible with respect to general inequalities 
    */

   int i, j, ncons;
   double alpha = 0.5;
   Constraint test;
   SerialDenseMatrix<int,double> grad_c;
   SerialDenseVector<int,double> g, gTg, type, v, g_tmp;

   for(i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     type.resize(test.getConstraintType().length());
     type   = test.getConstraintType();

     if(type(0) == Lineq || type(0) == NLineq){
       if(!test.amIFeasible(xc,ftol)){
	 v.resize(test.getConstraintViolation().length());
         v      = test.getConstraintViolation();
	 //	 grad_c.reshape(test.evalGradient(xc).numRows(),test.evalGradient(xc).numCols());
	 grad_c.reshape(xc.length(),v.length());
         grad_c = test.evalGradient(xc);
         if(type(0) == Lineq || type(0) == NLineq){
           ncons = v.length();
           gTg.resize(ncons);
           OptppArray<int> indices = test.getConstraintMappingIndices();
           for(j = 0; j <ncons-1; j++ ){
             if( std::abs(v(j)) > alpha ) {
	       int alpha = grad_c.numRows();
	       g.resize(alpha);
	       for(int i=0; i<alpha; i++)
		 { g(i) = grad_c(i,indices[j]);}
	       //g = grad_c(j);
               gTg(j) = g.dot(g);
	       g_tmp.resize(alpha);
	       g_tmp = g;
	       g_tmp *= (-v(j)/gTg(j));
	       xc += g_tmp;
             }
           }
         } 
       }
     }
   }
}


void CompoundConstraint::printConstraints() 
{
   /* 
    * Prints the raw value of the constraints.  
    */

   int i, j, index, ncons, nvars;
   Constraint test;
   SerialDenseVector<int,double> lower, upper, type, value;
   OptppArray<int> mapping;
   char s[2];

   for(i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     type.resize(test.getConstraintType().length());
     type   = test.getConstraintType();

     value.resize(test.getConstraintValue().length());
     value   = test.getConstraintValue();
     lower.resize(test.getLower().length());
     lower   = test.getLower();
     upper.resize(test.getUpper().length());
     upper   = test.getUpper();

     if(type(0) == Bound){
          cout <<"\nBound Constraints: \n";
     }
     else if(type(0) == NLeqn || type(0) == NLineq)
          cout <<"\nNonlinear Constraints: \n";
     else if(type(0) == Leqn || type(0) == Lineq)
          cout <<"\nLinear Constraints: \n";

     if(type(0) != Bound){
          ncons   = test.getNumOfCons();
          mapping = test.getConstraintMappingIndices();
          cout << "Index  Type       Lower   \t Constraint \t Upper \n";
          for (j = 1; j<=ncons; j++) {
               index = mapping[j-1];
               if(type(index-1) == NLeqn ||  type(index-1) == Leqn)   strcpy(s,"E");
               if(type(index-1) == NLineq || type(index-1) == Lineq)  strcpy(s,"I");
               cout << d(index,5) << "\t"   << s << "\t" 
                    << e(lower(index),12,4) << "\t" 
                    << e(value(index),12,4) <<  "\t" 
                    << e(upper(index),12,4) << "\n";
          }
     }
     else{
          nvars = getNumOfVars();
          cout << "Index \t Lower \t\t\t X \t Upper \n";
          for (j = 1; j<=nvars; j++) {
               cout << d(j,5) << "\t"   <<  e(lower(j),12,4) << "\t" 
                    << e(value(j),12,4) <<  "\t" << e(upper(j),12,4)  
                    << "\n";
          }
     }
   }
}

void CompoundConstraint::evalCFGH(const SerialDenseVector<int,double> & xc) const
{
   Constraint test;
   SerialDenseVector<int,double> result(numOfSets_);

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     test.evalCFGH(xc);
   }
}

// Evaluation 
SerialDenseVector<int,double> CompoundConstraint::evalResidual(const SerialDenseVector<int,double>& xc ) const 
{
 
   Constraint test;
   SerialDenseVector<int,double> result(numOfSets_);
   int alpha = 0;
   int beta;
   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     SerialDenseVector<int,double> temp(test.evalResidual(xc));
     beta = temp.length();
     result.resize(alpha+beta);
     //     SerialDenseVector<int,double> temp(beta);
     //     temp =  test.evalResidual(xc);
     for(int j=alpha; j<alpha+beta; j++)
       {result(j) = temp(j-alpha);}
     alpha = alpha+beta;

   }
  
  
   return result;
}

SerialDenseMatrix<int,double> CompoundConstraint::evalGradient(const SerialDenseVector<int,double>& xc ) const 
{
  SerialDenseMatrix<int,double> grad(1,1);
   Constraint test;
   int alpha = 0;
   int beta, chi;
   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     SerialDenseMatrix<int,double> temp(test.evalGradient(xc));
     beta = temp.numRows();
     chi = temp.numCols();
     //     SerialDenseMatrix<int,double> temp(beta,chi);
     //     temp = test.evalGradient(xc);
     grad.reshape(beta,alpha+chi);
     for(int j=0; j<beta; j++)
       for(int k = alpha; k<alpha+chi; k++)
	 {grad(j,k) = temp(j,k-alpha);}
     alpha = alpha + chi;

   }
   return grad;
}

SerialSymDenseMatrix<int,double> CompoundConstraint::evalHessian(SerialDenseVector<int,double>& xc ) const 
{
   
  // Extremely adhoc.  Conceived on 12/07/2000.  Vertical Concatenation
   SerialSymDenseMatrix<int,double> hessian(xc.length());
   hessian = 0;

   return hessian;
}

OptppArray<SerialSymDenseMatrix<int,double> > CompoundConstraint::evalHessian(SerialDenseVector<int,double>& xc, int darg ) const 
{
   
  // Extremely adhoc.  Conceived on 12/07/2000.  Vertical Concatenation
   SerialSymDenseMatrix<int,double> hessianT(xc.length());
   hessianT = 0;
   OptppArray<SerialSymDenseMatrix<int,double> > hessian(1);
   hessian[0] = hessianT;

   return hessian;
}

SerialSymDenseMatrix<int,double> CompoundConstraint::evalHessian(SerialDenseVector<int,double>& xc,
                                                const SerialDenseVector<int,double>& mult) const 
{
   int k, tncons;
   SerialSymDenseMatrix<int,double> hessian(xc.numRows()),temp2(3);
   OptppArray<SerialSymDenseMatrix<int,double> > temp;
   SerialDenseVector<int,double> type;
   Constraint test;

   k       = 0;
   hessian = 0.0;

   for(int i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     type.resize(test.getConstraintType().length());
     type   = test.getConstraintType();
     tncons = test.getNumOfCons();
     k     += tncons;

     if(type(0) == NLeqn || type(0) == NLineq){
        temp   = test.evalHessian(xc, i);
	for(int j = 0; j < temp.length(); j++){ 
	  temp2 = temp[j];
	  temp2 *= mult(j);
	  hessian += temp2;         
	  //hessian += temp[j].scale(mult(j+1));
        }
     }
   }
   return hessian;
}

bool CompoundConstraint::amIFeasible(const SerialDenseVector<int,double>& xc, double epsilon ) const
{
   bool feasible = true;
   SerialDenseVector<int,double> type;
   Constraint test;

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     type.resize(test.getConstraintType().length());
     type = test.getConstraintType();
     if(type(0) == Bound)
          feasible = test.amIFeasible(xc,epsilon);
     if(!feasible){
       //       OptppmathError("The current iterate is infeasible wrt to the bound and
       //                  linear constraints");
       break;
     }
   }
   return feasible;
}

void CompoundConstraint::insertSort(const OptppArray<Constraint>& constraints)
{
   Constraint ctemp;
   int dim = constraints.length();
   OptppArray<Constraint> sorted(dim);
   sorted = constraints;
   int i;

   if(dim == 1)
      constraints_ = sorted;
   else{
      for(int j = 1; j < dim; j++){
         ctemp = sorted[j];
	 i = j - 1;
	 while(i > -1 && compare(sorted[i],ctemp) > 0){
	    sorted[i+1] = sorted[i];
	    i--;
         }
	 sorted[i+1] = ctemp;
      }
      constraints_ = sorted;
  }
}

void CompoundConstraint::insertSort()
{
   Constraint ctemp;
   int dim = constraints_.length();
   int i;

   if(dim > 1){
      for(int j = 1; j < dim; j++){
         ctemp = constraints_[j];
	 i = j - 1;
	 while(i > -1 && compare(constraints_[i],ctemp) > 0){
	    constraints_[i+1] = constraints_[i];
	    i--;
         }
	 constraints_[i+1] = ctemp;
      }
   }
}

int CompoundConstraint::compare(const Constraint& c1, const Constraint& c2)
{
   SerialDenseVector<int,double> ct1 = c1.getConstraintType();
   SerialDenseVector<int,double> ct2 = c2.getConstraintType();
   
   if(ct1(0) < ct2(0))
      return -1;
   else if(ct1(0) > ct2(0))
     return 1;
   else
     return 0;
  
}

} // namespace OPTT
