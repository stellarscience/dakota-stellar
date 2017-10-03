#ifndef TestDDaceFactorialSampler_h
#define TestDDaceFactorialSampler_h

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "Distribution.h"
#include "UniformDistribution.h"
#include "DDaceFactorialSampler.h"
#include "suite.h"

class TestDDaceFactorialSampler : public Test
{
  public:
    TestDDaceFactorialSampler();
    ~TestDDaceFactorialSampler();

    void run();
    void testDDaceFactorialSampler4(); 
    void testDDaceFactorialSampler2(); 
    void testGetSamplesWithoutNoise();
    void testGetSamplesWithNoise();
    void testClone();
    void testPrint();
    void testTypeName();
    void testGetParameter();

  private:
    std::vector<Distribution>       dists;
    std::vector< std::vector<double> >    test_data;
    std::vector< std::vector<double> >    test_data_wn;
    std::vector<double>             lb;
    std::vector<double>             ub;
    int                       seed;
};

#endif
