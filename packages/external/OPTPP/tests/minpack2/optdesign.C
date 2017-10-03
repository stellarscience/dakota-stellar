// Optimal Design with Composite Materials from MINPACK-2 collection

#include <fstream>
#ifdef ANSI_HEADERS
#include <cstdio>
#include <cmath>
#include <cstring>
#else
#include <stdio.h>
#include <math.h>
#include <string.h>
#endif

#include "OptCG.h"
#include "NLF.h"
#include "minpack2.h"

using NEWMAT::ColumnVector;
using namespace OPTPP;

void update_model(int, int, ColumnVector) {}

int main ()
{
  int n = 100;
  ColumnVector x(n);
  
  static char *status_file = {"optdesign.out"};

// Create an NLF1 object and initialize it

  NLF1 nlp(n,optdesign,init_optdesign);
  
  nlp.setX(x);
  nlp.eval();               //  evaluate the function at x
  //  nlp.printState("Initial Guess for Optimal Design Problem");

//  Create a "Tolerances" object and set the tolerances
  
  TOLS tol;
  tol.setDefaultTol();
  tol.setFTol(1.e-9);
  tol.setMaxIter(200);

//  Build an optimization object and optimize 
  OptCG objfcn(&nlp,tol);
  
  objfcn.setOutputFile(status_file,0);
  objfcn.setUpdateModel(update_model);
  objfcn.optimize();
  objfcn.printStatus("Solution");
}
void init_optdesign(int n, ColumnVector& x)
{ /* set x = to the standard guess */

  int nx = (int) sqrt((double) n);
  int ny = n/nx;

  double fx;
  ColumnVector g;
  MINPACK_TASK   task = XStandard;
  double lambda = 0.008;

  dodcfg(nx,ny,x,fx,g,task,lambda);
}
void optdesign(int mode, int n, const ColumnVector& xc, double& fx, 
	       ColumnVector& g, int &result)
{ 

// Mode = 0  Function only
//        1  Gradient
//        2  Hessian
//        3  Function and gradient

  ColumnVector x = xc;
  double lambda = 0.008;
  int nx, ny;

  nx = (int) sqrt((double) n);
  ny = n/nx;

  MINPACK_TASK task;

  if (mode & NLPFunction) {
    task = Function;
    dodcfg(nx,ny,x,fx,g,task,lambda);
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    task = Gradient;
    dodcfg(nx,ny,x,fx,g,task,lambda); 
    result = NLPGradient;
 }
}
