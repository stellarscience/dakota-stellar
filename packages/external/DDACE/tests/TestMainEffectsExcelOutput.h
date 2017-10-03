#ifndef MainEffectsExcelOutput_h
#define MainEffectsExcelOutput_h

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <new>



#include "suite.h"
#include "Factor.h"
#include "Response.h"
#include "MainEffectsExcelOutput.h"




class TestMainEffectsExcelOutput : public Test
{
  public:
    TestMainEffectsExcelOutput();
    ~TestMainEffectsExcelOutput();
   void run();
   void testComputeExcelOutput();
};

#endif
