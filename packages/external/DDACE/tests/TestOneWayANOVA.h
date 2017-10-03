#ifndef TestOneWayANOVA_h
#define TestOneWayANOVA_h

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "suite.h"
#include "Factor.h"
#include "Response.h"
#include "OneWayANOVA.h"

class TestOneWayANOVA : public Test
{
  public:
    TestOneWayANOVA();
    ~TestOneWayANOVA();

    void run();
    void testConstructor();
    void testGetANOVATables();

  private:
     
};

#endif
