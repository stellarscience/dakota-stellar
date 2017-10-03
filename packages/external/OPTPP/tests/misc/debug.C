// Example of using a Newton-like optimize class on an FDNLF1 problem class
// This class is used whenever the objective function does not
// have analytic gradients
//

#include <iostream>
#include "OptQNewton.h"
#include "NLF.h"

using NEWMAT::ColumnVector;

using namespace OPTPP;

void init_rosen(int n, ColumnVector& x);
void rosen0(int n, const ColumnVector& x, double& fx, int& result);

int main ()
{
  int n = 2;
  ColumnVector x(n);
  int debug;

  //  Create a Nonlinear problem object

  FDNLF1 nlp(n,rosen0,init_rosen);
  
  //  Initialize and evaluate the function at x

  nlp.initFcn();
  nlp.eval();
  
  //  Build a Quasi-Newton object and optimize 

  OptQNewton objfcn(&nlp);   
  objfcn.setDebug();
  objfcn.optimize();
}
