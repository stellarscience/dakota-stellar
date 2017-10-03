#ifndef TestNormalDistribution_h
#define TestNormalDistribution_h

//#include "DDace.h"
#include "Distribution.h"
#include "NormalDistribution.h"
#include "suite.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

class TestNormalDistribution : public Test
{
  public:
    TestNormalDistribution();
    ~TestNormalDistribution();

    void run();
    void testNormalDistributionSigmaWithNumDeviations(); 
    void testNormalDistributionSigmaWithoutNumDeviations();
    void testNormalDistributionBoundsWithoutNumDeviations();
    void testNormalDistributionBoundsWithNumDeviations();
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
    Mean*           mean;
    StdDeviation*   stdDev;
    std::vector<double>   data;

    double          lb;
    double          ub;
    double          lb_est;
    double          ub_est;
    double          sigma;
    int             seed;
    int             numDevs;
    double          deviates[4];
};

#endif
