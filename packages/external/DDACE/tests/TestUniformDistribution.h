#ifndef TestUniformDistribution_h
#define TestUniformDistribution_h

#include "Distribution.h"
#include "UniformDistribution.h"
#include "UniformDistribution.h"
#include "suite.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

class TestUniformDistribution : public Test
{
  public:
    TestUniformDistribution();
    ~TestUniformDistribution();

    void run();
    void testUniformDistributionNoBounds(); 
    void testUniformDistributionWithBounds(); 
    void testClone();
    void testGetDeviateNoProb();
    void testGetDeviateWithProb();
    void testGetCDF();
    void testLowerBound();
    void testUpperBound();
    void testMean();
    void testStdDev();
    void testPrint();
    void testPrintAttributes();
    void testTypeName();

  private:
    double  lb;
    double  ub;
    int     seed;
    double  deviates[2];
};

#endif
