// Steady-State Combustion from MINPACK-2 collection

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
  
  static char *status_file = {"steady.out"};

// Create an NLF1 object and initialize it

  NLF1 nlp(n,steady_comb,init_steady_comb);
  
  nlp.setX(x);
  nlp.eval();               //  Evaluate the function at x
  
  TOLS tol;                 //  Create a "Tolerances" object and 
  tol.setDefaultTol();      //  set the tolerances
  tol.setFTol(1.e-9);
  tol.setMaxIter(200);

  OptCG objfcn(&nlp,tol);   //  Build a Quasi-Newton object and optimize 
  
  objfcn.setOutputFile(status_file,0);
  objfcn.setUpdateModel(update_model);
  objfcn.optimize();
  objfcn.printStatus("Solution");
}
void init_steady_comb(int n, ColumnVector& x)
{
  int nx = (int) sqrt((double) n);
  int ny = n/nx;
  double fx;
  ColumnVector g;
  MINPACK_TASK task = XStandard;
  double lambda = 5.0;

  dsscfg(nx,ny,x,fx,g,task,lambda);
}

void steady_comb(int mode, int n, const ColumnVector& xc, double& fx, ColumnVector& g, int& result)
{ 

// Mode = 0  Function only
//        1  Gradient
//        2  Hessian
//        3  Function and gradient
//        4  XS

  ColumnVector x = xc;
  double lambda = 5.0;
  int nx, ny;

  nx = (int) sqrt((double) n);
  ny = n/nx;

  MINPACK_TASK task;

  if (mode & NLPFunction) {
    task = Function;
    dsscfg(nx,ny,x,fx,g,task,lambda);
    result = NLPFunction;

  }
  if (mode & NLPGradient) {
    task = Gradient;
    dsscfg(nx,ny,x,fx,g,task,lambda); 
    result = NLPGradient;
 }
}
