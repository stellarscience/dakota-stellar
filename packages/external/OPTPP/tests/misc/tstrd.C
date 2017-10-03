//
// Test program for Readoptinput
//
// 1. Quasi Newton with trust regions on an NLF1
// 2. Quasi Newton with More Line Search on an NLF1
// 3. Quasi Newton with Backtracking Line Search on an NLF1
//

#include <string>
#include <iostream>
#include <fstream>

#include "OptNewton.h"
#include "OptQNewton.h"
#include "OptFDNewton.h"
#include "NLF.h"


using NEWMAT::ColumnVector;

using namespace OPTPP;

void ReadOptInput();

int main ()
{
  int n = 2;
  ColumnVector x(n);
  
  char *status_file = {"tstread.out"};
  ofstream opt_fp;
  opt_fp.open(status_file);

  ReadOptInput();

}

void ReadOptInput() // Read opt.input file if it exists
{
//  NLP1* nlp = nlprob();

/* A VERY simple routine for reading the optimization parameters
 * We should really make this more general, but as a first pass this
 * will have to do.
 * 
 * The input file should be of the form keyword = value
 * where keyword is one of the following
 * 
 * search      = trustregion
 * diff_option = forward
 * max_iter    = 100
 * max_feval   = 1000
 * grad_tol    = 1.e-6
 * fcn_tol     = 1.e-9
 * max_step    = 100.0
 * fcn_accrcy  = 1.e-9
 *
 */
  
  int  max_iter, max_feval;
  real grad_tol,  fcn_tol, max_step, fcn_accrcy;

  char token[80], ignore[80];
//
// Keywords allowed
//
  string keyword;
  string equals("=");
  string cdiff_option("diff_option");
  string cfcn_accrcy("fcn_accrcy");
  string cfcn_tol("fcn_tol");
  string cgrad_tol("grad_tol");
  string cmaxfeval("maxfeval");
  string cmaxiter("maxiter");
  string cmax_step("max_step");
  string csearch("search");

  string diff_option;
  string search;
  SearchStrategy s = TrustRegion;
  
  int keyword_count = 0;

// 
// Default name of input file
//
  char *opt_input  = {"opt.input"};
  int optin_eof;
//
// Open opt.input file and check to see if we succeeded
//

  ifstream optin(opt_input);
  if (!optin.rdbuf()->is_open()) {
    cout << "ReadOptInput: No opt.input file found\n";
    cout << "ReadOptInput: Default values will be used\n";
    return;
  }

  cout << "ReadOptInput: Reading opt.input file\n";

  optin >> token;

  while (!optin.eof()) {

    keyword = token;
    keyword_count++;

    cout << keyword_count << " keyword = " << keyword << "\n";

    if (keyword == cdiff_option) {

      optin >> equals >> token;
      diff_option = token;

    }    
    else if (keyword == cfcn_accrcy) {
      optin >> equals >> fcn_accrcy;
//      nlp->SetFcnAccrcy(fcn_accrcy);
    }    
    else if (keyword == cfcn_tol) {
      optin >> equals >> fcn_tol;
//      SetFcnTol(fcn_tol);
    }    
    else if (keyword == cgrad_tol) {
      optin >> equals >> grad_tol;
//      SetGradTol(grad_tol);
    }    
    else if (keyword == cmaxfeval) {
      optin >> equals >> max_feval;
//      SetMaxFeval(max_feval);
    }    
    else if (keyword == cmaxiter) {
      optin >> equals >> max_iter;
//      SetMaxIter(max_iter);
    }
    else if (keyword == cmax_step) {
      optin >> equals >> max_step;
//      SetMaxStep(max_step);
    }
    else if (keyword == csearch) {
      optin >> equals >> token;
      search = token;
      if ( search == "trustregion")
	s = TrustRegion;
      else if ( search == "linesearch")
	s = LineSearch;
//      SetSearchStrategy(s);
    }
    else {
      cout << "Unrecognized keyword '" << keyword << "'. "
	<< "Skipping the rest of this line\n";
      optin.getline(ignore, sizeof(ignore));
    }
  optin >> token;
  }

  cout << "\n\n======  Summary of input file  ======\n\n";

  cout << csearch      << " = " << search << "\n";
  cout << cdiff_option << " = " << diff_option << "\n";
  cout << cmaxiter     << " = " << max_iter << "\n";
  cout << cmaxfeval    << " = " << max_feval << "\n";
  cout << cgrad_tol    << " = " << grad_tol << "\n";
  cout << cfcn_tol     << " = " << fcn_tol << "\n";
  cout << cmax_step    << " = " << max_step << "\n";
  cout << cfcn_accrcy  << " = " << fcn_accrcy << "\n";

}
