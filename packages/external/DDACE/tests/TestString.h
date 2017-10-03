#ifndef TestString_h
#define TestString_h

#include <cstdio>
#include <cstdlib>
#include "SmartPtr.h"
#include <iostream>
#include <string>
#include "String.h"
#include "Array.h"
#include "suite.h"


class TestString : public Test {
 public:
  TestString();
  ~TestString();

  void run();
  void testStringConstructors();
  void testNewStringConstructors();
  void testAccessors();
  void testOperators();
  void testSearch();
  void testMiscFuncs();
};

#endif
