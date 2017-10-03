#ifndef TestDDaceSampler_h
#define TestDDaceSampler_h

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "DDaceSampler.h"
#include "DDaceSamplePoint.h"
#include "suite.h"


class TestDDaceSampler : public Test {
 public:
  TestDDaceSampler();
  ~TestDDaceSampler();

  void run();
  void testDDaceSampler();
  void testDDaceSamplerBase();
  void testGetSamples();
  void testPrint();
  void testTypeName();
  void testNSamples();
  void testNInputs();
  void testGetParameter();
  void testDist();
  void testLowerBounds();
  void testUpperBounds();
  void testNoise();

 private:
  std::vector< std::vector<double> >  data_d;   // data as Array< Array<double> >
  std::vector<DDaceSamplePoint> data_sp;  // data as Array<DDaceSamplePoint>
  std::vector<double> lb;                 // lower bounds
  std::vector<double> ub;                 // upper bounds
};

#endif
