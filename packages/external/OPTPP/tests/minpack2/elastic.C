// Elastic-Plastic Torsion function from MINPACK-2 collection


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
  int n = 16;
  
  static char *status_file = {"elastic.out"};

// Create an NLF1 object and initialize it

  NLF1 nlp(n,elastic,init_elastic);

  nlp.initFcn();
  nlp.eval();

  nlp.printState("Initial Guess for Elastic-Plastic Torsion Problem");
  

//  Build an optimization object and optimize 

  OptCG objfcn(&nlp);   

  objfcn.setOutputFile(status_file, 0);
  objfcn.setUpdateModel(update_model);
  objfcn.setStepTol(1.e-6);
  objfcn.setFcnTol(1.e-9);
  objfcn.setMaxIter(200);
  objfcn.optimize();
    
  objfcn.printStatus("Solution");
}

void init_elastic(int n, ColumnVector& x)
{ /* set x = to the standard guess */

  int nx = (int) sqrt((double) n);
  int ny = n/nx;
  double fx;
  ColumnVector g;
  MINPACK_TASK task = XUpper;
  double c   = 10.0;

  deptfg(nx,ny,x,fx,g,task,c);
  
}

void elastic(int mode, int n, const ColumnVector& xc, double& fx, ColumnVector& g, int& result)
{ 

// Mode = 0  Function only
//        1  Gradient
//        2  Hessian
//        3  Function and gradient

  ColumnVector x = xc;
  double c   = 10.0;
  int nx, ny;

  nx = (int) sqrt((double) n);
  ny = n/nx;

  MINPACK_TASK task;

  if (mode & NLPFunction) {
    task = Function;
    deptfg(nx,ny,x,fx,g,task,c);
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    task = Gradient;
    deptfg(nx,ny,x,fx,g,task,c); 
    result = NLPGradient;
 }
}
