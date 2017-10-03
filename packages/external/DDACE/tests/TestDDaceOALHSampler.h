#ifndef TestDDaceOALHSampler_h
#define TestDDaceOALHSampler_h

#include "Distribution.h"
#include "UniformDistribution.h"
#include "DDaceOALHSampler.h"
#include "suite.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

///////////////////////////
// Written By: J Cramp
// Date: July 2004
///////////////////////////

class TestDDaceOALHSampler : public Test
{
  public:
    TestDDaceOALHSampler();
    ~TestDDaceOALHSampler();

    void run();
    void testDDaceOALHSampler(); 
    void testGetSamples();
    void testClone();
    void testPrint();
    void testTypeName();
    void testGetParameter();
    void testGetDesign();
    void testGetOA();

  private:
    std::vector< std::vector<int> >       test_oa;
    std::vector< std::vector<int> >       test_design;
    std::vector< std::vector<double> >    test_data;
    double                    lb;
    double                    ub;
    int                       seed;
};

#endif
