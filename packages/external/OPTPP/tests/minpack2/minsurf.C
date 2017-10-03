// Minimal Surfaces Problem from MINPACK-2 collection

#include <fstream>
#ifdef ANSI_HEADERS
#include <cstdio>
#include <cmath>
#include <string>
#else
#include <stdio.h>
#include <math.h>
#include <string.h>
#endif

#include "ioformat.h"
#include "OptCG.h"
#include "NLF.h"
#include "minpack2.h"

using NEWMAT::ColumnVector;
using namespace OPTPP;

void update_model(int, int, ColumnVector) {}

int main ()
{
  int n = 25;
  
  static char *status_file = {"minsurf.out"};

// Create an NLF1 object and initialize it

  NLF1 nlp(n,minsurf,init_minsurf);
  
  nlp.initFcn();
  nlp.eval(); 

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
void init_minsurf(int n, ColumnVector& x)
{ /* set x = to the standard guess */

  int nx = (int) sqrt((double) n);
  int ny = n/nx;

  double fx;
  ColumnVector g;

  ColumnVector bottom(nx+2), top(nx+2), right(nx+2), left(nx+2);
  MINPACK_TASK   task = XStandard;

  dmsabc(nx,ny,bottom,top,left,right);
  dmsafg(nx,ny,x,fx,g,task,bottom,top,left,right);

}
void minsurf(int mode, int n, const ColumnVector& xc, double& fx, ColumnVector& g, int& result)
{ 

// Mode = 0  Function only
//        1  Gradient
//        2  Hessian
//        3  Function and gradient

  int nx, ny;

  nx = (int) sqrt((double) n);
  ny = n/nx;

  ColumnVector x = xc;
  ColumnVector bottom(nx+2), top(nx+2), right(nx+2), left(nx+2);
  dmsabc(nx,ny,bottom,top,left,right);

  MINPACK_TASK task;

  if (mode & NLPFunction) {
    task = Function;
    dmsafg(nx,ny,x,fx,g,task,bottom,top,left,right);
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    task = Gradient;
    dmsafg(nx,ny,x,fx,g,task,bottom,top,left,right); 
    result = NLPGradient;
 }
}
