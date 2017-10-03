#ifndef TEST_STAT_H
#define TEST_STAT_H

#include "suite.h"
#include "Statistics.h"
#include "arrcmp.h"

class TestStatistics : public Test
{
  public:
	TestStatistics();
	~TestStatistics();

	void run();
	void testAll();
};

#endif
