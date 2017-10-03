#ifndef TestDistribution_h
#define TestDistribution_h

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "Distribution.h"
#include "UniformDistribution.h"
#include "suite.h"

class TestDistribution : public Test
{
  public:
    TestDistribution();
    ~TestDistribution();

    void run();
    void testUsePseudoRandom();
    void testDistributionDefault(); 
    void testDistributionArg(); 
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
    double  deviate;    
};

#endif
