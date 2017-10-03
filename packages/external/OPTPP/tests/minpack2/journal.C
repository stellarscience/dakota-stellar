// Journal Bearing function from MINPACK-2 collection

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
  
  static char *status_file = {"journal.out"};

// Create an NLF1 object and initialize it

  NLF1 nlp(n,journal,init_journal);
  
  nlp.initFcn();
  nlp.eval(); 
  //nlp.setIsExpensive(true);
  nlp.printState("Initial Guess for Journal Bearing Problem");


//  Create a "Tolerances" object and set the tolerances

  TOLS tol;         
  tol.setDefaultTol();
  tol.setFTol(1.e-9);
  tol.setMaxIter(200);

//  Build an optimization object and optimize 

  OptCG objfcn(&nlp,tol);

  objfcn.setOutputFile(status_file, 1);
  objfcn.setUpdateModel(update_model);
  objfcn.optimize();
  objfcn.printStatus("Solution");
}
void init_journal(int n, ColumnVector& x)
{ /* set x = to the standard guess */

  int nx = (int) sqrt((double) n);
  int ny = n/nx;
  double fx;
  ColumnVector g;

  MINPACK_TASK task = XStandard;
  double ecc  = 0.1;
  double b    = 10.0;

  dpjbfg(nx,ny,x,fx,g,task,ecc,b);

}

void journal(int mode, int n, const ColumnVector& xc, double& fx, ColumnVector& g, int& result)
{ 

// Mode = 0  Function only
//        1  Gradient
//        2  Hessian
//        3  Function and gradient

  ColumnVector x = xc;
  double ecc = 0.1;
  double b   = 10.0;
  int nx, ny;

  nx = (int) sqrt((double) n);
  ny = n/nx;

  MINPACK_TASK task;

 
  if (mode & NLPFunction) {
    task = Function;
    dpjbfg(nx,ny,x,fx,g,task,ecc,b);
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    task = Gradient;
    dpjbfg(nx,ny,x,fx,g,task,ecc,b); 
    result = NLPGradient;
 }
}
