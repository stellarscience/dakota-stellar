#ifndef TestDDaceOASampler_h
#define TestDDaceOASampler_h

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
//#include "DDace.h"
#include "Distribution.h"
#include "UniformDistribution.h"
#include "DDaceOASampler.h"
#include "suite.h"

class TestDDaceOASampler : public Test
{
  public:
    TestDDaceOASampler();
    ~TestDDaceOASampler();

    void run();
    void testDDaceOASamplerWithDist(); 
    void testDDaceOASamplerWithoutDist(); 
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
