
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

#ifdef OPTPP_HAVE_MPI
#include "mpi.h"
#endif

#include "GenSet.h"
#include "OptGSS.h"
#include "NLF.h"

#include "tstfcn.h"

using Teuchos::SerialDenseVector;
using std::cerr;
using std::endl;

using namespace OPTPP;

void SetupTestProblem(std::string test_id, USERFCN0 *test_problem, 
		      INITFCN *init_problem);
void update_model(int, int, SerialDenseVector<int,double>) {}

int main (int argc, char* argv[])
{
  int ndim = 2;

  // USERFCN0 test_problem;
  // INITFCN  init_problem;

  // string test_id;

#ifdef OPTPP_HAVE_MPI
  int me;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

  //if (argc != 3) {
  //  cout << "Usage: tstGSS problem_name ndim\n";
  //  exit(1);
  //}

  //test_id = argv[1];
  //ndim    = atoi(argv[2]);

  //ColumnVector x(ndim);
  //ColumnVector vscale(ndim);

  //SetupTestProblem(test_id, &test_problem, &init_problem);

  char status_file[80];
  //strcpy(status_file, test_id.c_str());
  strcpy(status_file,"tstGSS");
#ifdef OPTPP_HAVE_MPI
  sprintf(status_file,"%s.out.%d", status_file, me);
#else
  strcat(status_file,".out");
#endif

  //  Create a Nonlinear problem object
  //NLF0 nlp(ndim, test_problem, init_problem);
  NLF0 nlp(ndim, erosen, init_erosen);
  
  // Create the GeneratingSet object
  GenSetStd gs(ndim);
  
  //  Build an optimization object and optimize 
  OptGSS optobj(&nlp, &gs);   

  optobj.setFullSearch(true);
  // optobj.setUpdateModel(update_model);

  if (!optobj.setOutputFile(status_file, 0))
    cerr << "main: output file open failed" << endl;

  optobj.optimize();

  optobj.printStatus("Final Status:",true);

#ifdef REG_TEST
  ColumnVector x_sol = nlp.getXc();
  double f_sol = nlp.getF();
  ostream* optout = optobj.getOutputFile();
  double x_tol = 1e-1; // high tolerance due to slow convergence
  double f_tol = 1e-2;
  if ((1.0 - x_sol(0) <= x_tol) && (1.0 - x_sol(1) <= x_tol) && (f_sol
								 <=
								 f_tol))
    *optout << "GSS 1 PASSED" << endl;
  else
    *optout << "GSS 1 FAILED" << endl;
#endif

  optobj.cleanup();

#ifdef OPTPP_HAVE_MPI
  MPI_Finalize();
#endif    

}

