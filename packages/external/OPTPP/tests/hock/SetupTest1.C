//------------------------------------------------------------------------
// Copyright (C) 1993, 1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include <string>
#include <iostream>

#include "NLF.h"
#include "tstfcn.h"

using namespace OPTPP;

void SetupTestProblem(string test_id, USERFCN1 *test_problem, 
		      INITFCN *init_problem)
{
  //
  // Given a character string test_id, set the function pointers
  // test_problem and init_problem 
  // 

  string tstfcn[2];
  USERFCN1 test_problems[2];
  INITFCN  init_problems[2];

  bool found = false;

  tstfcn       [0] = "rosen";
  test_problems[0] = rosen;
  init_problems[0] = init_rosen;

  tstfcn       [1] = "erosen";
  test_problems[1] = erosen1;
  init_problems[1] = init_erosen;

  int maxprob = 3;

  //
  // Loop over possible test problems
  //
  int i = 0;
  while (!found & i < maxprob) {
    if (test_id == tstfcn[i]) {
      test_id      = tstfcn[i];
      *test_problem = test_problems[i];
      *init_problem = init_problems[i];
      found = true;
    }
    i++;
  }
  if (!found) {
    cout << "SetupTestProblem: No problem specified. Using rosen as default\n";
    *test_problem = rosen;
    *init_problem = init_rosen;
  }
}

