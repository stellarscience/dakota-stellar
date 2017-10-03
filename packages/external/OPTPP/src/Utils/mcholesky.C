//JWG

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#define WANT_MATH

#include <iostream>
#include <cfloat>

//#include "include.h"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

using Teuchos::SerialSymDenseMatrix;
using Teuchos::SerialDenseMatrix;

namespace OPTPP {
//------------------------------------------------------------------------
//
// Perturbed  Cholesky decomposition 
//
//------------------------------------------------------------------------

static double square(double x) { return x*x; }
SerialDenseMatrix<int,double> PertChol(SerialSymDenseMatrix<int,double>&, double, double&);

SerialDenseMatrix<int,double> MCholesky(SerialSymDenseMatrix<int,double>& S)
{
  //   Tracer trace("MCholesky");
   using std::sqrt;
   using std::fabs;
   int nr = S.numRows();
   SerialDenseMatrix<int,double> L(nr,nr);
   double mcheps = DBL_EPSILON;

   //   Real* s = S.Store(); Real* l = L.Store();

   double maxadd = 0.0;

   int i, j;

   double sqrteps = sqrt(mcheps);
   double maxdiag = 0.0;
   double mindiag = 1.0e10;
   double maxoff  = 0.0;
   for (i=0; i<nr; ++i) {
     maxdiag = max(maxdiag,S(i,i));
     mindiag = min(mindiag,S(i,i));
     if(i != nr)
       {for (j=i; j<=i; ++j) 
	   {maxoff = max(maxoff,S(i,j));
	   }
       }
   }


   double maxposdiag = max(0.0,maxdiag);
   double mu;

   if (mindiag <= sqrteps*maxposdiag) {
     mu = 2.0*(maxposdiag-mindiag)*sqrteps - mindiag;
     maxdiag = maxdiag + mu;
   }
   else mu = 0.0;

   if (maxoff*(1.0 + 2.0*sqrteps) > maxdiag) {
     mu = mu + (maxoff-maxdiag) + 2.0*sqrteps*maxoff;
     maxdiag = maxoff * (1.0 + 2.0*sqrteps);
   }

   if (maxdiag == 0.0) {
     mu = 1.0;
     maxdiag = 1.0;
   }
   if (mu > 0.0) {
     for (i=0; i<nr; ++i) S(i,i) = S(i,i) + mu;
   }

   double maxoffl = sqrt(max(maxdiag,(maxoff/nr)));

   L = PertChol(S,maxoffl,maxadd);
   
   
   if (maxadd > 0.0) {

     double maxev = S(0,0);
     double minev = S(0,0);
     for (i=0; i<nr; ++i) {
       double offrow = 0.0;
       for(j=0; j<=i-1; ++j) offrow += fabs(S(j,i));
       for(j=i+1; j<nr; ++j) offrow += fabs(S(i,j));
       maxev = max(maxev,(S(i,i)+offrow));
       minev = min(minev,(S(i,i)-offrow));
      }
     double sdd = (maxev -minev) * sqrteps - minev;
     sdd = max(sdd,0.0);
     mu = min(maxadd,sdd);
     for (i=0; i<nr; ++i) S(i,i) = S(i,i) + mu;
     
     L = PertChol(S,0.0,maxadd);
   }
       
   return L;
}



SerialDenseMatrix<int,double> PertChol(SerialSymDenseMatrix<int,double>& S, double maxoffl, double& maxadd)
{
  using std::sqrt;
  using std::fabs;
  int i;
  //  Tracer trace("PertChol");
  int nr = S.numRows();
  SerialDenseMatrix<int,double> L(nr,nr);
  double mcheps = DBL_EPSILON;
  
  //  Real* s = S.Store(); Real* l = L.Store();
  double sum;
  double minl2  = 0.0;
  
  int j, k;
  double minl = std::pow(mcheps,.25)*maxoffl;
  
  if (maxoffl == 0.0) {
    double maxdiag = 0.0;
    for (i=0;i<nr;++i) maxdiag = max(maxdiag,fabs(S(i,i)));
    maxoffl = sqrt(maxdiag);
    minl2 = sqrt(mcheps)*maxoffl;
  }
  maxadd = 0.0;
  
  for (j=0; j<nr; j++) {
    sum = 0.0;
    for (i=0; i<=j-1; ++i) {
      sum += L(j,i)*L(j,i); 
    }
    double ljj = S(j,j) - sum;
    
    double  minljj = 0.0;
    
    for (i=j+1; i<nr; ++i) {
      sum = 0.0;
      for(k=0; k<=j-1; ++k) {
	sum += L(i,k) * L(j,k);
      }
      L(i,j) = S(j,i) - sum;
      minljj = max(fabs(L(i,j)),minljj);
    }
    minljj = max((minljj/maxoffl),minl);
    
    if (ljj > square(minljj)) { // Normal Cholesky
      L(j,j) = sqrt(ljj);
    }
    else {//    Modify ljj since it is too small
      if (minljj < minl2) minljj = minl2;
      maxadd = max(maxadd,(square(minljj)-ljj));
      L(j,j) = minljj;
    }
    for (i=j+1; i<nr; ++i){
      L(i,j) = L(i,j) / L(j,j);
    }
    
  }
  return L;
}

} // namespace OPTPP

