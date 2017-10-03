#if defined(SGI) || defined(RS6K)
#define WANT_MATH
#else
#define WANT_STREAM
#define WANT_MATH
#endif

#include "include.h"

#include "newmatap.h"

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void Clean(Matrix&, Real);
void Clean(DiagonalMatrix&, Real);


int main()
{
  Matrix A(5,5);
  double a[25];
  
  DiagonalMatrix SV(5), D(5), D1(5), D2(5); 
  Matrix V;

  SymmetricMatrix H(5); 

  int i, j, n = 5;

  cout << "Enter elements of A: " << "\n";
  for (i=0; i<n; i++)  
    for (j=0; j<n; j++)
      cin >> a[j + n*(i)];

  A << a;
  H << A;

  Print(H);
//  SymmetricMatrix H; H << A;
  
  SVD(A, SV);
  Print(SV);

  SVD(H, D1);
  Print(D1);

  EigenValues(H, D2, V);
  Print(D2);
  Print(V);
}
