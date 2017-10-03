//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#if (defined(__sgi) || defined(__xlc__) || defined(__xlC__))
#define WANT_MATH
#else
#define WANT_STREAM
#define WANT_MATH
#endif


#include "newmat.h"
#include "newmatrc.h"

#include "cblas.h"

using namespace std;


#define REPORT {}
//#define REPORT { static ExeCounter ExeCount(__LINE__,8); ExeCount++; }


//Real ColumnVector::Dot(ColumnVector& y)
//{
//  REPORT                                         // not accessed
//  int n1 = storage; int n2 = y.storage;
//  if (n1 != n2) {
//    cerr << "ColumnVector::Dot: Vectors not of the same size\n";
//    return 0.0;
//  }
//  if (n1<=0) return 0.0;
//  
//  Real* el1=store; Real* el2=y.store;
//  Real sum = ddot(n1,el1,1,el2,1);
//  return sum;
//}


