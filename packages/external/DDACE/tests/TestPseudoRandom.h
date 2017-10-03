#ifndef TESTPSEUDORANDOM_H_
#define TESTPSEUDORANDOM_H_


#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>

#include "PseudoRandomTestsOnly.h"
#include "suite.h"



class TestPseudoRandom : public Test {
	public:
	    TestPseudoRandom();
	    ~TestPseudoRandom();
        void run();
        void testPseudoRandom();
	    void testGetPseudoRandom();
	    void testSetSeed();   
};

#endif 

