
#include <iostream>
#include <fstream>

#include "NLF.h"
#include "OptFDNIPS.h"

#include "hockfcns.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_RCP.hpp"



using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;

using namespace OPTPP;

void update_model(int, int, SerialDenseVector<int,double>) {}

int main ()
{
  int n = 2;

  static char *status_file = {"tsthock10.out"};

  SerialDenseMatrix<int,double> ATest(2,1);
  SerialDenseVector<int,double> gtp(1),yzm(1);
 //  ATest(0,0) = 5.0;
//   ATest(1,0) = 1.0;
//   yzm(0) = 7.0;
//   std::cout<<"ATest = "<<ATest<<std::endl;
//   std::cout<<"yzm = "<<yzm<<std::endl;
//   gtp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, ATest, yzm, 0.0);
//   std::cout<<"gtp = "<<gtp<<std::endl;


  //  Create a Constrained Nonlinear problem object
  NLF1 nips(n,hs10,init_hs10,create_constraint_hs10);

  //  Build a finite-difference NIPS object and optimize
  OptFDNIPS objfcn(&nips, update_model);
  objfcn.setOutputFile(status_file, 0);
  objfcn.setFcnTol(1.0e-06);
  objfcn.setMaxIter(50);
  objfcn.setSearchStrategy(LineSearch);
  objfcn.setMeritFcn(ArgaezTapia);
  objfcn.optimize();
  objfcn.printStatus("Solution from nips");
  objfcn.cleanup();

#ifdef REG_TEST
  SerialDenseVector<int,double> x_sol = nips.getXc();
  double f_sol = nips.getF();
  ostream* optout = objfcn.getOutputFile();
  if ((0.0 - x_sol(0) <= 1.e-2) && (1.0 - x_sol(1) <= 1.e-2) && (-1.0 - f_sol
								 <=
								 1.e-2))
    *optout << "Hock  10 PASSED" << endl;
  else
    *optout << "Hock  10 FAILED" << endl;
#endif
}
