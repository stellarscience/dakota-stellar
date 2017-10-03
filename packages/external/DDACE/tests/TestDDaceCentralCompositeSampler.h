#ifndef TestDDaceCentralCompositeSampler_h
#define TestDDaceCentralCompositeSampler_h

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "Distribution.h"
#include "UniformDistribution.h"
#include "DDaceCentralCompositeSampler.h"
#include "suite.h"

class TestDDaceCentralCompositeSampler : public Test
{
  public:
    TestDDaceCentralCompositeSampler();
    ~TestDDaceCentralCompositeSampler();

    void run();
    void testDDaceCentralCompositeSampler();
    void testGetSamples();
    void testClone();
    void testPrint();
    void testTypeName();
  
  private:
    std::vector<Distribution>        dists;
    std::vector< std::vector<double> >     test_data;
};

#endif
