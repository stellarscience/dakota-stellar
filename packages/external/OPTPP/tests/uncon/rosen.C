#include "NLP.h"
#include "ioformat.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

using namespace OPTPP;

extern double alpha;
static int fcn_count = 0;
/* Example file to demonstrate the calling sequence to a 
 * simple NLF1 function
 */
void init_rosen (int ndim, SerialDenseVector<int,double>& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  x(0) = -1.2;
  x(1) =  1.0;

}
void rosen(int mode, int n, const SerialDenseVector<int,double>& x, double& fx, SerialDenseVector<int,double>& g, int& result)
{ // Rosenbrock's function
  double f0, f1, x0, x1;
  
  if (n != 2) return;

  x0 = x(0);
  x1 = x(1);
  f0 = (x1 - x0 * x0);
  f1 = 1. - x0;
  
  if (mode & NLPFunction) {
    fx  = 100.* f0*f0 + f1*f1;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(0) = -400.*f0*x0 - 2.*f1; 
    g(1) = 200.*f0;
    result = NLPGradient;
  }
}
void rosen0(int n, const SerialDenseVector<int,double>& x, double& fx, int& result)
{ // Rosenbrock's function
  double f0, f1, x0, x1;

  if (n != 2) return;

  x0 = x(0);
  x1 = x(1);
  f0 = (x1 - x0 * x0);
  f1 = 1. - x0;

  
  fx  = 100.* f0*f0 + f1*f1;

  result = NLPFunction;
}
void rosen2(int mode, int n, const SerialDenseVector<int,double>& x, double& fx, 
	SerialDenseVector<int,double>& g, SerialSymDenseMatrix<int,double>& H, int &result)
// Rosenbrock's function, n = 2 with first and second order derivatives

{ 
  fcn_count++;

  double f0, f1, x0, x1;

  // cout << "\nrosen2: mode = " << mode
  //      << " count = " << fcn_count << "\n";
  // for(int i=1; i<=n; i++) 
  //     cout << "x(" << i << ") = " << x(i) << "\n";

  if (n != 2) return;

  x0 = x(0);
  x1 = x(1);
  f0 = (x1 - x0 * x0);
  f1 = 1. - x0;
  
  if (mode & NLPFunction) {
    fx  = 100.* f0*f0 + f1*f1;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(0) = -400.*f0*x0 - 2.*f1; 
    g(1) = 200.*f0;
    result = NLPGradient;
  }
  
  if (mode & NLPHessian) {
    f0 = (x1 - 3.0*x0*x0);
    H(0,0) = -400.0*f0 + 2.0;
    H(1,0) = -400.0*x0;
    H(1,1) = 200.0;
    result = NLPHessian;
  }
}

void rosen0_least_squares(int n, const SerialDenseVector<int,double>& x, SerialDenseVector<int,double>& fx, int& result)
{ // Rosenbrock's function
  double x0, x1;

  if (n != 2) return;

  x0 = x(0);
  x1 = x(1);
  fx(0) = 10*(x1 - x0 * x0);
  fx(1) = 1. - x0;
  
  result = NLPFunction;
}

void rosen_least_squares(int mode, int n, const SerialDenseVector<int,double>& x, SerialDenseVector<int,double>& fx, 
	       SerialDenseMatrix<int,double>& gx, int& result)
{ // Rosenbrock's function
  double x0, x1;

  if (n != 2) return;

  x0 = x(0);
  x1 = x(1);
  if (mode & NLPFunction) {
    fx(0) = 10*(x1 - x0 * x0);
    fx(1) = 1. - x0;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    gx(0,0) = -20*x0; 
    gx(0,1) =  10; 
    gx(1,0) = -1.0;
    gx(1,1) =  0.0;
    result = NLPGradient;
  }
}
