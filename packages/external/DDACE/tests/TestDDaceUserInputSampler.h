#ifndef TestDDaceUserInputSampler_h
#define TestDDaceUserInputSampler_h

#include "Distribution.h"
#include "UniformDistribution.h"
#include "DDaceUserInputSampler.h"
#include "suite.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

class TestDDaceUserInputSampler : public Test
{
  public:
    TestDDaceUserInputSampler();
    ~TestDDaceUserInputSampler();

    void run();
    void testDDaceUserInputSampler();  
    void testGetSamples();
    void testClone();
    void testPrint();
    void testTypeName();
    void testGetParameter();

  private:
    std::vector< std::vector<double> >    test_data;
};

#endif
