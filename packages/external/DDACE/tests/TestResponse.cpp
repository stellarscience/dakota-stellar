#include "TestResponse.h"
#include "arrcmp.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

TestResponse::TestResponse() {
	std::vector<double> data(10);
	data[0] = 11.06;
	data[1] = 15.4;
	data[2] = 11.81;
	data[3] = 20.48;
	data[4] = 20.32;
	data[5] = 10.96;
	data[6] = 15.86;
	data[7] = 10.8;
	data[8] = 19.73;
	data[9] = 20.57;
	
	this->response = DDaceMainEffects::Response(data);
}


TestResponse::~TestResponse() {
}



void TestResponse::run()
{
    testNullConstructor();
    testConstructorWithVectorOfDoubles();
    testCopyConstructor();
    testGetAveragePop();
    testGetSumOfSquaresPop();
    testGetVariancePop();
    testgetNumOfObservations();
    testBracketNotation();  
    testGetSumPop();   
}


void TestResponse::testNullConstructor() {
	DDaceMainEffects::Response response;
	_test(response.getNumOfObservations()==0);
}


void TestResponse::testConstructorWithVectorOfDoubles(){

	_test(fabs(this->response.getAveragePop()-15.699) < 0.000000001);
	_test(this->response.getNumOfObservations()==10);
	_test(fabs(this->response.getSumOfSquaresPop()-167.405) < 0.001);
	_test(fabs(this->response.getVariancePop()-18.6006) < 0.0001);
	
}


void TestResponse::testCopyConstructor() {
	DDaceMainEffects::Response response(this->response);
	_test(fabs(response.getAveragePop()-15.699) < 0.000000001);
	_test(response.getNumOfObservations()==10);
	_test(fabs(response.getSumOfSquaresPop()-167.405) < 0.001);
	_test(fabs(response.getVariancePop()-18.6006) < 0.0001);
	
}

void TestResponse::testGetSumPop(){

	_test(fabs(this->response.getSumPop()-156.99) < 0.000000001);
}


void TestResponse::testGetAveragePop(){

	_test(fabs(this->response.getAveragePop()-15.699) < 0.000000001);
}


void TestResponse::testGetSumOfSquaresPop() {

	_test(fabs(this->response.getSumOfSquaresPop()-167.405) < 0.001);
}


void TestResponse::testGetVariancePop(){
	_test(fabs(this->response.getVariancePop()-18.6006) < 0.0001);	
}


void TestResponse::testgetNumOfObservations(){
	_test(this->response.getNumOfObservations()==10);
}

void TestResponse::testBracketNotation(){
	_test(fabs(this->response[0] - 11.06) < 0.01);
    _test(fabs(this->response[1] - 15.4) < 0.01);
    _test(fabs(this->response[2] - 11.81) < 0.01);
	_test(fabs(this->response[3] - 20.48) < 0.01);
	_test(fabs(this->response[4] - 20.32) < 0.01);
	_test(fabs(this->response[5] - 10.96) < 0.01);
	_test(fabs(this->response[6] - 15.86) < 0.01);
	_test(fabs(this->response[7] - 10.8) < 0.01);
	_test(fabs(this->response[8] - 19.73) < 0.01);
	_test(fabs(this->response[9] - 20.57) < 0.01);
}





