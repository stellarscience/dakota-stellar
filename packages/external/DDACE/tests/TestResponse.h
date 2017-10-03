#ifndef TestResponse_h
#define TestResponse_h

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "suite.h"
#include "Response.h"


class TestResponse : public Test
{
  public:
    TestResponse();
    ~TestResponse();

    void run();
    void testNullConstructor();
    void testConstructorWithVectorOfDoubles();
    void testCopyConstructor();
    void testGetAveragePop();
    void testGetSumOfSquaresPop();
    void testGetVariancePop();
    void testgetNumOfObservations();
    void testBracketNotation();
    void testGetSumPop();

  protected:
       DDaceMainEffects::Response response;
     
};

#endif
