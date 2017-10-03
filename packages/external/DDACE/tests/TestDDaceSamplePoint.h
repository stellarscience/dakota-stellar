#ifndef TestDDaceSamplePoint_h
#define TestDDaceSamplePoint_h

#include <cstdio>
#include <cstdlib>
#include "SmartPtr.h"
#include <iostream>
#include <fstream>

#include <string>
#include "DDaceSamplePoint.h"
#include "suite.h"
#include <math.h>


class TestDDaceSamplePoint : public Test {
 public:
  TestDDaceSamplePoint();
  ~TestDDaceSamplePoint();

  void run();
  void testDDaceSamplePoint();
  void testDDaceSamplePointDouble();
  void testDDaceSamplePointString();
  void testGetDataType();
  void testSetDataType();
  void testIndex();
  void testLength();
  void testOperator();
  void testGetDoubleValue();
  void testGetStringValue();
  void testParameters();
  void testGetDoubleValues();
  void testGetStringValues();
  void testPrint();
};

#endif

