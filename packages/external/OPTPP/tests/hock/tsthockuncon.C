
#include <iostream>
#include <fstream>

#include "NLF.h"
#include "OptFDNIPS.h"
#include "hockfcns.h"

using NEWMAT::ColumnVector;

using namespace OPTPP;

void update_model(int, int, ColumnVector) {}

int main ()
{
  int n = 2;

  static char *status_file = {"tsthockuncon.out"};
  ofstream opt_fp;
  int opt_fd;


  //  Create a Constrained Nonlinear problem object
  NLF1 nips(n,hsuncon,init_hsuncon);

  //  Build a NIPS object and optimize
  OptNIPS objfcn(&nips, update_model);
  opt_fd = objfcn.setOutputFile(status_file, 0);
  objfcn.setFcnTol(5.0e-8);
  objfcn.setMaxIter(20);
  objfcn.setSearchStrategy(LineSearch);
  objfcn.setMeritFcn(NormFmu);
  objfcn.optimize();
  objfcn.printStatus("Solution from NIPS");
  objfcn.cleanup();


}

