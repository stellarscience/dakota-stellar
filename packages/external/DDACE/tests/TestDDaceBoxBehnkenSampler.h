#ifndef TestDDaceBoxBehnkenSampler_h
#define TestDDaceBoxBehnkenSampler_h

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "Distribution.h"
#include "UniformDistribution.h"
#include "DDaceBoxBehnkenSampler.h"
#include "suite.h"

class TestDDaceBoxBehnkenSampler : public Test
{
  public:
    TestDDaceBoxBehnkenSampler();
    ~TestDDaceBoxBehnkenSampler();

    void run();
    void testDDaceBoxBehnkenSampler();
    void testGetSamples();
    void testClone();
    void testPrint();
    void testTypeName();
  
  private:
    std::vector<Distribution>        dists;
    std::vector< std::vector<double> >     test_data;
};

#endif
