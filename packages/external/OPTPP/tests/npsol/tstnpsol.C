//
// Test program for NPSOL objects
//

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#include <fstream>
#ifdef HAVE_STD
#include <cstdio>
#else
#include <stdio.h>
#endif

#include "NLF.h"
#include "OptNPSOL.h"

#include "tstfcn.h"

using NEWMAT::ColumnVector;

using namespace OPTPP;

double const LOW_BOUND =  0.0E0;
double const UP_BOUND  =  1.0E0;

main ()
{
  int          i, n=7;
  ColumnVector x(n);
  
  char *status_file = {"tstnpsol.out"};
  ofstream opt_fp;
  opt_fp.open(status_file);

  int ndim  = n;
  int nclin = 0;
  int ncnln = 0;
  ColumnVector lower(n), upper(n);

  //  Create boundary constraints 

  for(i=1; i<=ndim; i++) {lower(i) = LOW_BOUND;  upper(i) = UP_BOUND;}

  //  Build a optimization object and optimize 

  OptNPSOL objfcn(ndim,nclin,ncnln,illum,init_illum,illum);   
  objfcn.setX(x);
  objfcn.setLower(lower);
  objfcn.setUpper(upper);
  objfcn.setOutput(opt_fp);
  objfcn.optimize();
  objfcn.printStatus("Solution from NPSOL ");
  opt_fp.close();
}

