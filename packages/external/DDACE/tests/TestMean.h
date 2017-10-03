#ifndef TestMean_h
#define TestMean_h

#include <cstdio>
#include <cstdlib>
#include "SmartPtr.h"
#include <iostream>
#include <string>
#include "Mean.h"
#include "suite.h"


class TestMean : public Test {
 public:
  TestMean();
  ~TestMean();

  void run();
  void testMean();
  void testMeanDouble();
  void testMeanArray();
};

#endif
