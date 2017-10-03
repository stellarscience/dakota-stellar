#ifndef TestDDaceLHSampler_h
#define TestDDaceLHSampler_h

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "Distribution.h"
#include "UniformDistribution.h"
#include "DDaceLHSampler.h"
#include "suite.h"

class TestDDaceLHSampler : public Test
{
  public:
    TestDDaceLHSampler();
    ~TestDDaceLHSampler();

    void run();
    void testDDaceLHSamplerWithDist(); 
    void testDDaceLHSamplerWithoutDist(); 
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
