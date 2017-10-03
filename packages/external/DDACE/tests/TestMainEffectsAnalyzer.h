#ifndef TestMainEffectsAnalyzer_h
#define TestMainEffectsAnalyzer_h

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "ColumnHeader.h"
#include "suite.h"
#include "ColumnHeader.h"
#include "MainEffectsAnalyzer3.h"

class TestMainEffectsAnalyzer : public Test
{
  public:
    TestMainEffectsAnalyzer();
    ~TestMainEffectsAnalyzer();

    void run();
    void testAllMainEffects();

  private:
	std::vector< std::vector < DataValue > > 	data;
	std::vector< ColumnHeader > 		columnHeaders;
	MainEffectsAnalyzer3					x;
    	int                       				seed;
};

#endif
