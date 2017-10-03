// Example of using a PD optimize class on an NLF0 problem class
// This class is used whenever the objective function does not
// have analytic gradients
//

#include <fstream>

#include "OptPDS.h"
#include "NLF.h"
#include "tstfcn.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;

using namespace OPTPP;

extern "C" {
extern void init_twafer(int, double *);
} 

main ()
{
  int i, j;
  int ndim = 4;
  double perturb;
  ColumnVector x(ndim);
  ColumnVector vscale(ndim);
  Matrix init_simplex(ndim,ndim+1);

  char *schemefilename = {"myscheme"};
  char *status_file = {"pds-twaf.out"};
  double powers[4];
  ofstream opt_fp;
  int opt_fd;

  //  Create a Nonlinear problem object

  NLF0 nlp(ndim, twaf, init_twaf);         
  
  //  Initialize and evaluate the function at x

  nlp.initFcn();

  x = nlp.getXc();
  for (i = 0; i < ndim; i++)
    powers[i] = x(i+1);

  /* Create input file needed by TWAFER. */

  init_twafer(ndim, powers);

  nlp.eval();               
  
  //  Build a PDS object and optimize 
  
  OptPDS objfcn(&nlp);      
  opt_fd = objfcn.setOutputFile(status_file);

  objfcn.setFcnTol(1.e-5);
  objfcn.setMaxIter(10);
  objfcn.setSSS(10);

  vscale = x.MaximumAbsoluteValue();
  objfcn.setScale(vscale);

  for (i=1; i <= ndim; i++) {
    for (j=1; j <= ndim+1; j++) {
      init_simplex(i,j) = x(i);
    }
  }

  for (i=1; i<= ndim; i++) {
    perturb = x(i)*.01;
    init_simplex(i,i+1) = x(i) + perturb;
  }

  objfcn.setSimplexType(2);
  objfcn.setSimplex(init_simplex);

  objfcn.setCreateFlag();

  objfcn.setSchemeFileName(schemefilename);

  objfcn.optimize();
  
  objfcn.printStatus("Solution from PDS");

  objfcn.Cleanup();
}
