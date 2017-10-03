#ifndef TestMainEffectsConverter_h
#define TestMainEffectsConverter_h

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
#include "MainEffectsConverter.h"

class TestMainEffectsConverter : public Test
{
  public:
    TestMainEffectsConverter();
    ~TestMainEffectsConverter();

    void run();
    void testConstructor();
    
    void testConvert();
    void testConvertTableOfDoublesToArray();
    void testConvertAllDoublesToCountingNumbers();
    void testSliceOutOneInputVarAndOneOutputVar();
};

#endif
