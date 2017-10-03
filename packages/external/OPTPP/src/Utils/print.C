//Deal with IsZero command
//might have screwed up the indexing for the print statements

//JWG

#if (defined(__sgi) || defined(__xlc__) || defined(__xlC__))
#define WANT_MATH
#else
#define WANT_STREAM
#define WANT_MATH
#endif

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#include "ioformat.h"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

using Teuchos::SerialSymDenseMatrix;
using Teuchos::SerialDenseMatrix;

namespace OPTPP {

class PrintCounter
{
   int count;
   const char* s;
public:
   ~PrintCounter();
   PrintCounter(const char * sx) : count(0), s(sx) {}
   void operator++() { count++; }
};

PrintCounter PCZ("Number of non-zero matrices (should be 1) = ");
PrintCounter PCN("Number of matrices tested                 = ");

PrintCounter::~PrintCounter()
//jcm{ cout << s << count << "\n"; }
{}


void Print(const SerialDenseMatrix<int,double>& X)
{
   ++PCN;
   //   cout << "\nPrint::Matrix type: " << X.Type().Value() << " (";
   cout << X.numRows() << ", ";
   cout << X.numCols() << ")\n\n";
   // if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.numRows(); int nc=X.numCols();
   for (int i=0; i<nr; i++)
   {
      for (int j=0; j<nc; j++)  cout << e(X(i,j),14,6) << "\t";
      cout << "\n";
   }
   cout << flush; ++PCZ;
}


void Print(const SerialSymDenseMatrix<int,double>& X)
{
   ++PCN;
   // cout << "\nMatrix type: " << X.Type().Value() << " (";
   cout << X.numRows() << ", ";
   cout << X.numCols() << ")\n\n";
   //if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.numRows(); int nc=X.numCols();
   int i, j;
   for (i=0; i<nr; i++)
   {
      for (j=0; j<i; j++) cout << e(X(j,i),14,6) << "\t";
      for (j=i; j<nc; j++)  cout << e(X(i,j),14,6) << "\t";
      cout << "\n";
   }
   cout << flush; ++PCZ;
}




void Clean(SerialDenseMatrix<int,double>& A, double c)
{
   int nr = A.numRows(); int nc = A.numCols();
   for (int i=0; i<nr; i++)
   {
      for ( int j=0; j<nc; j++)
      { double a = A(i,j); if ((a < c) && (a > -c)) A(i,j) = 0.0; }
   }
}



void FPrint(ostream *fout, const SerialDenseMatrix<int,double>& X)
{
  //(*fout) << "\nFPrint::Matrix type: " << X.Type().Value() << " (";
   (*fout) << X.numRows() << ", ";
   (*fout) << X.numCols() << ")\n\n";
   // if (X.IsZero()) { (*fout) << "All elements are zero\n" << flush; return; }
    int nr=X.numRows(); int nc=X.numCols();
   for (int i=0; i<nr; i++)
   {
      for (int j=0; j<nc; j++)  (*fout) << e(X(i,j),14,6) << "\t";
      (*fout) << "\n";
   }
   (*fout) << flush; ++PCZ;
}

void FPrint(ostream *fout, const SerialSymDenseMatrix<int,double>& X)
{
   ++PCN;
   //(*fout) << "\nFPrint::Matrix type: " << X.Type().Value() << " (";
   (*fout) << X.numRows() << ", ";
   (*fout) << X.numCols() << ")\n\n";
   // if (X.IsZero()) { (*fout) << "All elements are zero\n" << flush; return; }
   int nr=X.numRows(); int nc=X.numCols();
   int i, j;
   for (i=0; i<nr; i++)
   {
      for (j=0; j<i; j++) (*fout) << e(X(j,i),14,6) << "\t";
      for (j=i; j<nc; j++)  (*fout) << e(X(i,j),14,6) << "\t";
      (*fout) << "\n";
   }
   (*fout) << flush; ++PCZ;
}

} // namespace OPTPP
