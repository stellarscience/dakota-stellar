#ifndef TestStdDeviation_h
#define TestStdDeviation_h

#include <cstdio>
#include <cstdlib>
#include "SmartPtr.h"
#include <iostream>
#include <string>
#include "StdDeviation.h"
#include "arrcmp.h"
#include "suite.h"


class TestStdDeviation : public Test {
 public:
  TestStdDeviation();
  ~TestStdDeviation();

  void run();
  void testStdDeviationDouble();
  void testStdDeviationArray();
};

#endif
