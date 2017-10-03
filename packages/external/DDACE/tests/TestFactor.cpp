#include "TestFactor.h"
#include "arrcmp.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;
    
	TestFactor::TestFactor(){ }
	
	TestFactor::~TestFactor(){ }
	
	void TestFactor::run() {
		testConstructor();
	    testSumOfSquaresWithinGroups();
	    testSumOfSquaresBetweenGroups();
	    testVarianceBetweenGroups();
	    testVarianceWithinGroups();
	    testFdata();
	    testGetNumberOfLevels();
	    testGetNumberOfObservations();
	    testDoFWithin();
	    testDoFBetween();
	    testGetAllLevelAverages();	
	    testGetLevelSum();
	    testGetLevelAverage();
	    testGetLevelSumOfSquares();
	}
	
	void TestFactor::testConstructor(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;   
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	        DDaceMainEffects::Factor(dataInput, 2, response);    	
		std::vector<int> factors = factor.getFactors();
		_test(factors.size()==10);
		for (int i=0; i<10; i++)
		    _test(factors[i] == dataInput[i]);
		response = factor.getResponse();
		_test(response.getNumOfObservations()==10);
		for (int i=0; i<10; i++)
		    _test(response[i] == dataOutput[i]);
		    
	
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;   
	    response = DDaceMainEffects::Response(dataOutput);
	    factor = DDaceMainEffects::Factor(dataInput, 15, response);    	
		factors = factor.getFactors();
		_test(factors.size()==10);
		for (int i=0; i<10; i++)
		    _test(factors[i] == dataInput[i]);
		response = factor.getResponse();
		_test(response.getNumOfObservations()==10);
		for (int i=0; i<10; i++)
		    _test(response[i] == dataOutput[i]);	    
		    
		
	}
	
	
	void TestFactor::testGetLevelSumOfSquares(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);    		      
		_test(fabs(factor.getLevelSumOfSquares(0)-27.3789) < 0.0001);
        _test(fabs(factor.getLevelSumOfSquares(1)-0.4281) < 0.0001);
		
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
	    factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);    		
		_test(fabs(factor.getLevelSumOfSquares(0)-27.3789) < 0.0001);
        _test(fabs(factor.getLevelSumOfSquares(1)-0.4281) < 0.0001);	}
	
	
	void TestFactor::testGetLevelAverage(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);    		      
		_test(fabs(factor.getLevelAverage(0)-12.6483) < 0.0001);
        _test(fabs(factor.getLevelAverage(1)-20.275) < 0.0001);
		
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
	    factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);    		
		_test(fabs(factor.getLevelAverage(0)-12.6483) < 0.0001);
        _test(fabs(factor.getLevelAverage(1)-20.275) < 0.0001);
	}
	
	
	void TestFactor::testGetLevelSum(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);    		      
		_test(fabs(factor.getLevelSum(0)-75.89) < 0.0001);
        _test(fabs(factor.getLevelSum(1)-81.1) < 0.0001);
		
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
	    factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);    		
		_test(fabs(factor.getLevelSum(0)-75.89) < 0.0001);
        _test(fabs(factor.getLevelSum(1)-81.1) < 0.0001);
	}

	void TestFactor::testGetNumberOfObservations(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);    		
		_test(factor.getNumberOfObservations()==10);
		
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
	    factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);    		
		_test(factor.getNumberOfObservations()==10);		
	}
	
	void TestFactor::testGetNumberOfLevels(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);    		
	    _test(factor.getNumberOfLevels()==2);
	    
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
        factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);    		
	    _test(factor.getNumberOfLevels()==2);	    
	}	
	
	void TestFactor::testDoFBetween(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);	    
	    _test(factor.doFBetween()==1);
	    
	    
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
        factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);	    
	    _test(factor.doFBetween()==1);
	    
	}
		
	
    /*-------------------------------------------------------------------*/
    /*-------------------------------------------------------------------*/
    
	
	void TestFactor::testGetAllLevelAverages(){
				
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);	    
	    vector<double> levelAverages = factor.getAllLevelAverages();
	    _test(levelAverages.size()==2);    	
	    _test(fabs(levelAverages[0]-12.6483)<0.0001);
	    _test(fabs(levelAverages[1]-20.275)<0.001);

	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
       factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);	    
	    levelAverages = factor.getAllLevelAverages();	    
	    _test(levelAverages.size()==2);    	
	    _test(fabs(levelAverages[0]-12.6483)<0.0001);
	    _test(fabs(levelAverages[1]-20.275)<0.001);
	    	    
	    
	}	
	


    /*-------------------------------------------------------------------*/
    /*-------------------------------------------------------------------*/
    
	
	void TestFactor::testSumOfSquaresWithinGroups(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);	    		
		_test(fabs(factor.sumOfSquaresWithinGroups()-27.807)<0.001);
		
		dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
        factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);	    		
		_test(fabs(factor.sumOfSquaresWithinGroups()-27.807)<0.001);
		
	}

	
	
    /*-------------------------------------------------------------------*/
    /*-------------------------------------------------------------------*/
    
	
	void TestFactor::testSumOfSquaresBetweenGroups(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);
		_test(fabs(factor.sumOfSquaresBetweenGroups()-139.599)<0.001);
		
		
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;
	    response = DDaceMainEffects::Response(dataOutput);
        factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);
		_test(fabs(factor.sumOfSquaresBetweenGroups()-139.599)<0.001);		
	}
		
	
	
	
	
	

	
    /*-------------------------------------------------------------------*/
    /*-------------------------------------------------------------------*/
    
	
	
	void TestFactor::testVarianceBetweenGroups(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);    		
		_test(fabs(factor.varianceBetweenGroups()-139.599)<0.001);
		
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
        factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);    		
		_test(fabs(factor.varianceBetweenGroups()-139.599)<0.001);		
		
	}
	
	
    /*-------------------------------------------------------------------*/
    /*-------------------------------------------------------------------*/
    
	
	
	void TestFactor::testVarianceWithinGroups(){
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);    		
		_test(fabs(factor.varianceWithinGroups()-3.47587)<0.00001);	
		
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
        factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);    		
		_test(fabs(factor.varianceWithinGroups()-3.47587)<0.00001);		
	}
	
	

	
	
    /*-------------------------------------------------------------------*/
    /*-------------------------------------------------------------------*/
    
	
	void TestFactor::testFdata(){
		
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);    		
	    _test(fabs(factor.Fdata()-40.1621)<0.0001);
	    
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
        factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);    		
	    _test(fabs(factor.Fdata()-40.1621)<0.0001);	    
		
	}
	
	
    /*-------------------------------------------------------------------*/
    /*-------------------------------------------------------------------*/
    
	
	
	void TestFactor::testDoFWithin(){
		
		std::vector<int> dataInput(10);
	    std::vector<double> dataOutput(10);
	    dataInput[0] = 0;  dataOutput[0]=11.06;
	    dataInput[1] = 0;  dataOutput[1]=15.4;
	    dataInput[2] = 0;  dataOutput[2]=11.81;
	    dataInput[3] = 1;  dataOutput[3]=20.48;
	    dataInput[4] = 1;  dataOutput[4]=20.32;	    
	    dataInput[5] = 0;  dataOutput[5]=10.96;
	    dataInput[6] = 0;  dataOutput[6]=15.86;
	    dataInput[7] = 0;  dataOutput[7]=10.8;
	    dataInput[8] = 1;  dataOutput[8]=19.73;
	    dataInput[9] = 1;  dataOutput[9]=20.57;	    	
	    DDaceMainEffects::Response response(dataOutput);
	    DDaceMainEffects::Factor factor = 
	       DDaceMainEffects::Factor(dataInput, 2, response);    	
	    _test(factor.doFWithin()==8);
	    
	    dataInput[0] = 10;  dataOutput[0]=11.06;
	    dataInput[1] = 10;  dataOutput[1]=15.4;
	    dataInput[2] = 10;  dataOutput[2]=11.81;
	    dataInput[3] = 12;  dataOutput[3]=20.48;
	    dataInput[4] = 12;  dataOutput[4]=20.32;	    
	    dataInput[5] = 10;  dataOutput[5]=10.96;
	    dataInput[6] = 10;  dataOutput[6]=15.86;
	    dataInput[7] = 10;  dataOutput[7]=10.8;
	    dataInput[8] = 12;  dataOutput[8]=19.73;
	    dataInput[9] = 12;  dataOutput[9]=20.57;	    	
	    response = DDaceMainEffects::Response(dataOutput);
        factor = 
	       DDaceMainEffects::Factor(dataInput, 15, response);    	
	    _test(factor.doFWithin()==8);	    
	}
	
	
	







