#ifndef TestDDaceRandomSampler_h
#define TestDDaceRandomSampler_h

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "Distribution.h"
#include "UniformDistribution.h"
#include "DDaceRandomSampler.h"
#include "suite.h"

class TestDDaceRandomSampler : public Test
{
  public:
    TestDDaceRandomSampler();
    ~TestDDaceRandomSampler();

    void run();
    void testDDaceRandomSamplerWithDist(); 
    void testDDaceRandomSamplerWithoutDist(); 
    void testGetSamples();
    void testClone();
    void testPrint();
    void testTypeName();
    void testGetParameter();

  private:
    std::vector<Distribution>       dists;
    std::vector< std::vector<double> >    test_data;
    std::vector<double>             lb;
    std::vector<double>             ub;
    int                       seed;
};

#endif
