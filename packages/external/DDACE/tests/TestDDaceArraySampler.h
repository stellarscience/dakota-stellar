#ifndef TestDDaceArraySampler_h
#define TestDDaceArraySampler_h
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
//#include "DDace.h"
#include "DDaceSamplePoint.h"
#include "DDaceArraySampler.h"
#include "suite.h"


class TestDDaceArraySampler : public Test {
	public:
	    TestDDaceArraySampler();
	    ~TestDDaceArraySampler();
            void run();
	    void testDDaceSamplePoint();

};

#endif
