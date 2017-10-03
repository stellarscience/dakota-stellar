
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstdio>
#else
#include <stdio.h>
#endif

#include "NLF.h"

using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;
using namespace OPTPP;

extern double alpha;
static int fcn_count = 0;
/* Example file to demonstrate the calling sequence to a 
 * simple NLF1 function
 */
void init_ellipsoid (int ndim, ColumnVector& x)
{
  for (int i=1; i<=ndim; i++) x(i) = 0.0;
}

void ellipsoid(int n, const ColumnVector& x, double& fx, int& result)
{
  int i;
  
  fcn_count++;
  cout << "Ellipsoid : count = " << fcn_count << "\n";
  fx = 0.0;
  for (i=1; i<=n; i++) fx += (pow(2.0e0,1.0e0*i)*(x(i)-i)*(x(i)-i));
  result = NLPFunction;
}

void ellipsoid1(int mode, int n, const ColumnVector& x, double& fx, 
	       ColumnVector& g, int& result)
{
  int i;
  
  fcn_count++;
  cout << "Ellipsoid : count = " << fcn_count << "\n";
  fx = 0.0;
  for (i=1; i<=n; i++) {
    fx = fx + pow(2.0e0,1.0e0*i) * (x(i) - i) * (x(i) - i);
    g(i) = pow(2.0e0,1.0e0*i) * 2.0 * (x(i) - i);
  }
  result = NLPFunction | NLPGradient;
}

void ellipsoid2(int mode, int n, const ColumnVector& x, double& fx, 
	       ColumnVector& g, SymmetricMatrix&H, int& result)
{
  int i;
  
  fcn_count++;
  cout << "Ellipsoid : count = " << fcn_count << "\n";
  fx = 0.0;
  for (i=1; i<=n; i++) {
    fx = fx + pow(2.0e0,1.0e0*i) * (x(i) - i) * (x(i) - i);
    g(i) = pow(2.0e0,1.0e0*i) * 2.0 * (x(i) - i);
  }
  H = 0.0;
  for (i=1; i<=n; i++) H(i,i) = pow(2.0e0,1.0e0*i) * 2.0; 

  result = NLPFunction | NLPGradient | NLPHessian;
}

