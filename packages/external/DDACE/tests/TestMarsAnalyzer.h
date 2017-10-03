#ifndef TestMarsAnalyzer_h
#define TestMarsAnalyzer_h

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "Distribution.h"
#include "UniformDistribution.h"
#include "NormalDistribution.h"
#include "Mean.h"
#include "DDaceSampler.h"
#include "DDaceLHSampler.h"
#include "FuncApprox.h"
#include "Mars.h"
#include "suite.h"

class TestMarsAnalyzer : public Test
{
  public:
    TestMarsAnalyzer();
    ~TestMarsAnalyzer();

    void run();
    void testMarsOutput();
    bool close(double a, double b);

  private:
   std::vector<Distribution>       dists;
    std::vector<DDaceSamplePoint>   data;
    std::vector<double>	      y;
    std::vector<double>	      settings;
    int                       seed;

};

#endif
