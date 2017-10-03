#include "TestPseudoRandom.h"

TestPseudoRandom::TestPseudoRandom(){;}

TestPseudoRandom::~TestPseudoRandom(){;}


void TestPseudoRandom::run(){
	testPseudoRandom();
	testSetSeed();
	testGetPseudoRandom();
}


void TestPseudoRandom::testPseudoRandom(){
	PseudoRandomTestsOnly x;
	double number = x.getPseudoRandom();
	_test(number == 0.0);
}

void TestPseudoRandom::testGetPseudoRandom(){
	
	/* first number should be 0 */
	PseudoRandomTestsOnly x;
	double number = x.getPseudoRandom();
	_test(number == 0.0);
	
	/* 2nd number should be 0.001 */
	number = x.getPseudoRandom();
	_test(fabs(number-0.001) < 0.0001);
	
	/* 3rd number should be 0.002 */
	number = x.getPseudoRandom();
	_test(fabs(number-0.002) < 0.0001);
	
	/* 4th number should be 0.003 */
	number = x.getPseudoRandom();
	_test(fabs(number-0.003) < 0.0001);
	
	/* new number should be 0.998 */
	x.setSeed(998);
	number = x.getPseudoRandom();
	_test(fabs(number-0.998) < 0.0001);	
	
	/* 2nd number should be 0.999 */
	number = x.getPseudoRandom();
	_test(fabs(number-0.999) < 0.0001);		
	
	/* 3rd number should be 0.0 */
	number = x.getPseudoRandom();
	_test(fabs(number-0.0) < 0.0001);		
	
	/* 4th number should be 0.001 */
	number = x.getPseudoRandom();
	_test(fabs(number-0.001) < 0.0001);		
	
	/* 5th number should be 0.002 */
	number = x.getPseudoRandom();
	_test(fabs(number-0.002) < 0.0001);		
	
}


void TestPseudoRandom::testSetSeed() {
	
	
	/* first number should be 0 */
	PseudoRandomTestsOnly x;
	double number = x.getPseudoRandom();
	_test(number == 0.0);
	
	/* 2nd number should be 0.001 */
	number = x.getPseudoRandom();
	_test(fabs(number-0.001) < 0.0001);
	
	
	/* new number should be 0.500 */
	x.setSeed(500);
	number = x.getPseudoRandom();
	_test(fabs(number - 0.500) < 0.0001);
	
	/* 2nd number should be 0.501 */
	number = x.getPseudoRandom();
	_test(fabs(number-.501) < 0.0001);		
	
	/* new number should be 0.999 */
	x.setSeed(999);
	number = x.getPseudoRandom();
	_test(fabs(number - 0.999) < 0.0001);
	
	/* 2nd number should be 0.000 */
	number = x.getPseudoRandom();
	_test(fabs(number - 0.0) < 0.0001);		
}	
