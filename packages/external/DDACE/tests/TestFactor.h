#ifndef TestFactor_h
#define TestFactor_h

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "suite.h"
#include "Factor.h"
#include "Response.h"



class TestFactor : public Test
{
  public:
    TestFactor();
    ~TestFactor();
    


    void run();
    void testConstructor();
    void testSumOfSquaresWithinGroups();
    void testSumOfSquaresBetweenGroups();
    void testVarianceBetweenGroups();
    void testVarianceWithinGroups();
    void testFdata();
    void testGetNumberOfLevels();
    void testGetNumberOfObservations();
    void testDoFWithin();
    void testDoFBetween();
    void testGetAllLevelAverages();
    void testGetLevelSum();
    void testGetLevelAverage();
    void testGetLevelSumOfSquares();
    

  protected:
     
};

#endif
