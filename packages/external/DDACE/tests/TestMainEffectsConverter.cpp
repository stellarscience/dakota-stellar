#include "TestMainEffectsConverter.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <limits>

using namespace std;


TestMainEffectsConverter::
    TestMainEffectsConverter(){
}

TestMainEffectsConverter::
    ~TestMainEffectsConverter(){
}

void TestMainEffectsConverter::run() {
	
	testConstructor();    
	testConvert();
    testConvertTableOfDoublesToArray();
	testConvertAllDoublesToCountingNumbers();
	testSliceOutOneInputVarAndOneOutputVar();
}


void TestMainEffectsConverter::
    testConvertTableOfDoublesToArray() {
    	
    	MainEffectsConverter x;
    	
    	/* empty input data */
    	std::vector<std::vector<double> > vectorInputDataPoints;
        ValueAndRowIndexAndColumnIndex * array = 
    	    x.convertTableOfDoublesToArray
    	    (vectorInputDataPoints);
    	//Can't test for an empty array
    	delete[] array;
    	
    	/* one data point */
    	vectorInputDataPoints.clear();
    	std::vector<double> row1;
    	row1.push_back(5);
    	vectorInputDataPoints.push_back(row1);
        array = 
    	    x.convertTableOfDoublesToArray
    	    (vectorInputDataPoints);
    	_test(array[0].value==5);    
    	_test(array[0].indexRow==0);
    	_test(array[0].indexColumn==0);
        delete[] array;
    	
    	/* one row, two columns */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(5);    	
    	row1.push_back(6);
    	vectorInputDataPoints.push_back(row1);
        array = 
    	    x.convertTableOfDoublesToArray
    	    (vectorInputDataPoints);
    	    
    	_test(array[0].value==5);    	
    	_test(array[0].indexRow==0);
    	_test(array[0].indexColumn==0);
    	
    	_test(array[1].value==6);  
    	_test(array[1].indexRow==0);
    	_test(array[1].indexColumn==1);
    	    	
    	delete[] array;
    	
    	

    	/* one row, three columns */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(5);    	
    	row1.push_back(6);
    	row1.push_back(7);
    	vectorInputDataPoints.push_back(row1);
        array = 
    	    x.convertTableOfDoublesToArray
    	    (vectorInputDataPoints);
    	    
    	_test(array[0].value==5);    	
    	_test(array[0].indexRow==0);
    	_test(array[0].indexColumn==0);
    	    	    	
    	_test(array[1].value==6);  
    	_test(array[1].indexRow==0);
    	_test(array[1].indexColumn==1);
    	    	
    	_test(array[2].value==7);      	
    	_test(array[2].indexRow==0);
    	_test(array[2].indexColumn==2);
    	
    	
    	delete[] array;
    	
    	
    	/* two rows, one column */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(5);    	
        std::vector<double>row2;
        row2.push_back(6);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
        array = 
    	    x.convertTableOfDoublesToArray
    	    (vectorInputDataPoints);
    	    
    	_test(array[0].value==5);    	    	
    	_test(array[0].indexRow==0);
    	_test(array[0].indexColumn==0);
    	
    	_test(array[1].value==6);     
    	_test(array[1].indexRow==1);
    	_test(array[1].indexColumn==0);
    	    	 
    	delete[] array;
    	
    	
    	/* three rows in ZYX order, one column */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	std::vector<double>row3;
    	row1.push_back(7);    	
    	row2.push_back(6);
        row3.push_back(5);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
    	vectorInputDataPoints.push_back(row3);    	
        array = 
    	    x.convertTableOfDoublesToArray
    	    (vectorInputDataPoints);
    	    
    	_test(array[0].value==7);    	
    	_test(array[0].indexRow==0);
    	_test(array[0].indexColumn==0);
    	    	
    	_test(array[1].value==6);      	  		    	
    	_test(array[1].indexRow==1);
    	_test(array[1].indexColumn==0);
    	
    	_test(array[2].value==5);     	    	  		    	    	
    	_test(array[2].indexRow==2);
    	_test(array[2].indexColumn==0);
    	
    	delete[] array;
    	
    	
    	/* two rows, two columns in random order */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row1.push_back(5);    	
    	row1.push_back(7);
    	row2.push_back(5);
        row2.push_back(6);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
        array = 
    	    x.convertTableOfDoublesToArray
    	    (vectorInputDataPoints);
    	    
    	_test(array[0].value==5);    	
    	_test(array[0].indexRow==0);
    	_test(array[0].indexColumn==0);
    	    	
    	_test(array[1].value==7);      	  		    	
    	_test(array[1].indexRow==0);
    	_test(array[1].indexColumn==1);
    	
    	_test(array[2].value==5);      	  		    	    	
    	_test(array[2].indexRow==1);
    	_test(array[2].indexColumn==0);
    	
    	_test(array[3].value==6);      
    	_test(array[3].indexRow==1);
    	_test(array[3].indexColumn==1);
    	
    	delete[] array;
    	
    	/* three rows, three columns in random order */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row3.clear();
    	row1.push_back(5);    	
    	row1.push_back(7);
    	row1.push_back(6);
    	row2.push_back(5);
        row2.push_back(6);
        row2.push_back(5);
    	row3.push_back(7);
        row3.push_back(6);
        row3.push_back(7);        
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
        vectorInputDataPoints.push_back(row3);
        array = 
    	    x.convertTableOfDoublesToArray
    	    (vectorInputDataPoints);
    	    
    	_test(array[0].value==5);
    	_test(array[0].indexRow==0);
    	_test(array[0].indexColumn==0);
    	
    	_test(array[1].value==7);
    	_test(array[1].indexRow==0);
    	_test(array[1].indexColumn==1);
    	
        _test(array[2].value==6);
    	_test(array[2].indexRow==0);
    	_test(array[2].indexColumn==2);
        
    	_test(array[3].value==5);
    	_test(array[3].indexRow==1);
    	_test(array[3].indexColumn==0);
    	
    	_test(array[4].value==6);
    	_test(array[4].indexRow==1);
    	_test(array[4].indexColumn==1);
    	
        _test(array[5].value==5);
    	_test(array[5].indexRow==1);
    	_test(array[5].indexColumn==2);
        
    	_test(array[6].value==7);
    	_test(array[6].indexRow==2);
    	_test(array[6].indexColumn==0);
    	
    	_test(array[7].value==6);
    	_test(array[7].indexRow==2);
    	_test(array[7].indexColumn==1);
    	
        _test(array[8].value==7);
    	_test(array[8].indexRow==2);
    	_test(array[8].indexColumn==2);
        
        delete[] array;
}


void TestMainEffectsConverter::
    testConvertAllDoublesToCountingNumbers(){
    	
    	MainEffectsConverter x;
    	
    	/* empty input data */
    	std::vector<std::vector<double> > vectorInputDataPoints;
        VectorCountingNumbersAndCount vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        std::vector<std::vector<int> > vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        int count = vectorAndCount.count;
        _test(count==0);    	    

    	
    	/* one data point */
    	vectorInputDataPoints.clear();
    	std::vector<double> row1;
    	row1.push_back(5);
    	vectorInputDataPoints.push_back(row1);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==1);    	    
        _test(vectorCountingNumbers[0][0]==0);
    	
    	/* one row, two columns */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(5);    	
    	row1.push_back(6);
    	vectorInputDataPoints.push_back(row1);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==2);    	    
        _test(vectorCountingNumbers[0][0]==0);
        _test(vectorCountingNumbers[0][1]==1);
        
            	
    	/* one row, two columns in ZYX order */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(6);    	
    	row1.push_back(5);
    	vectorInputDataPoints.push_back(row1);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==2);    	    
        _test(vectorCountingNumbers[0][0]==1);
        _test(vectorCountingNumbers[0][1]==0);    	

    	/* one row, three columns */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(5);    	
    	row1.push_back(6);
    	row1.push_back(7);
    	vectorInputDataPoints.push_back(row1);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==0);
        _test(vectorCountingNumbers[0][1]==1);
        _test(vectorCountingNumbers[0][2]==2);
            	
    	/* one row, three columns in ZYX order */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(7);    	
    	row1.push_back(6);
    	row1.push_back(5);
    	vectorInputDataPoints.push_back(row1);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==2);
        _test(vectorCountingNumbers[0][1]==1);
        _test(vectorCountingNumbers[0][2]==0);
            	
    	/* one row, three columns in random order */    	
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(5);    	
    	row1.push_back(7);
    	row1.push_back(6);
    	vectorInputDataPoints.push_back(row1);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==0);
        _test(vectorCountingNumbers[0][1]==2);
        _test(vectorCountingNumbers[0][2]==1);
            	
    	/* one row, three columns in random order */    	
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(6);    	
    	row1.push_back(5);
    	row1.push_back(7);
    	vectorInputDataPoints.push_back(row1);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==1);
        _test(vectorCountingNumbers[0][1]==0);
        _test(vectorCountingNumbers[0][2]==2);
            	
    	/* one row, three columns in random order */    	
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(7);    	
    	row1.push_back(5);
    	row1.push_back(6);
    	vectorInputDataPoints.push_back(row1);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==2);
        _test(vectorCountingNumbers[0][1]==0);
        _test(vectorCountingNumbers[0][2]==1);    	
    	
    	/* two rows, one column */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row1.push_back(5);    	
        std::vector<double>row2;
        row2.push_back(6);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==2);    	    
        _test(vectorCountingNumbers[0][0]==0);
        _test(vectorCountingNumbers[1][0]==1);
            	
    	/* two rows in ZYX order, one column */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row1.push_back(6);    	
        row2.push_back(5);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==2);    	    
        _test(vectorCountingNumbers[0][0]==1);
        _test(vectorCountingNumbers[1][0]==0);
            	
    	/* three rows, one column */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row1.push_back(5);    	
    	row2.push_back(6);
        std::vector<double>row3;
        row3.push_back(7);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
    	vectorInputDataPoints.push_back(row3);    	
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==0);
        _test(vectorCountingNumbers[1][0]==1);
        _test(vectorCountingNumbers[2][0]==2);
        
            	
    	/* three rows in ZYX order, one column */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row3.clear();
    	row1.push_back(7);    	
    	row2.push_back(6);
        row3.push_back(5);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
    	vectorInputDataPoints.push_back(row3);    	
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==2);
        _test(vectorCountingNumbers[1][0]==1);
        _test(vectorCountingNumbers[2][0]==0);
            	
    	/* three rows in random order, one column */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row3.clear();
    	row1.push_back(5);    	
    	row2.push_back(7);
        row3.push_back(6);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
    	vectorInputDataPoints.push_back(row3);    	
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==0);
        _test(vectorCountingNumbers[1][0]==2);
        _test(vectorCountingNumbers[2][0]==1);
            	
    	/* three rows in ZYX order, one column */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row3.clear();
    	row1.push_back(6);    	
    	row2.push_back(5);
        row3.push_back(7);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
    	vectorInputDataPoints.push_back(row3);    	
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==1);
        _test(vectorCountingNumbers[1][0]==0);
        _test(vectorCountingNumbers[2][0]==2);
            	
    	/* three rows in ZYX order, one column */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row3.clear();
    	row1.push_back(6);    	
    	row2.push_back(7);
        row3.push_back(5);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
    	vectorInputDataPoints.push_back(row3);    	
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==1);
        _test(vectorCountingNumbers[1][0]==2);
        _test(vectorCountingNumbers[2][0]==0);
            	
    	/* three rows in ZYX order, one column */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row3.clear();
    	row1.push_back(7);    	
    	row2.push_back(5);
        row3.push_back(6);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
    	vectorInputDataPoints.push_back(row3);    	
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==2);
        _test(vectorCountingNumbers[1][0]==0);
        _test(vectorCountingNumbers[2][0]==1);
            	
    	/* two rows, two columns in random order */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row1.push_back(5);    	
    	row1.push_back(7);
    	row2.push_back(5);
        row2.push_back(6);
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==0);
        _test(vectorCountingNumbers[0][1]==2);
        _test(vectorCountingNumbers[1][0]==0);
        _test(vectorCountingNumbers[1][1]==1);
            	
    	/* three rows, three columns in random order */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row3.clear();
    	row1.push_back(5);    	
    	row1.push_back(7);
    	row1.push_back(6);
    	row2.push_back(5);
        row2.push_back(6);
        row2.push_back(5);
    	row3.push_back(7);
        row3.push_back(6);
        row3.push_back(7);        
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
        vectorInputDataPoints.push_back(row3);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==0);
        _test(vectorCountingNumbers[0][1]==2);
        _test(vectorCountingNumbers[0][2]==1);
        _test(vectorCountingNumbers[1][0]==0);
        _test(vectorCountingNumbers[1][1]==1);    	
        _test(vectorCountingNumbers[1][2]==0);        
        _test(vectorCountingNumbers[2][0]==2);
        _test(vectorCountingNumbers[2][1]==1);    	
        _test(vectorCountingNumbers[2][2]==2);        
        
        
    	/* three rows, three columns, with noise in random order */
    	vectorInputDataPoints.clear();
    	row1.clear();
    	row2.clear();
    	row3.clear();
    	row1.push_back(5.0001);    	
    	row1.push_back(6.9999);
    	row1.push_back(6);
    	row2.push_back(4.9999);
        row2.push_back(6.0001);
        row2.push_back(5);
    	row3.push_back(7);
        row3.push_back(5.9999);
        row3.push_back(7.0001);        
    	vectorInputDataPoints.push_back(row1);
    	vectorInputDataPoints.push_back(row2);
        vectorInputDataPoints.push_back(row3);
        vectorAndCount = 
    	    x.convertAllDoublesToCountingNumbers
    	    (vectorInputDataPoints);
        vectorCountingNumbers =
            vectorAndCount.vectorCountingNumbers;
        count = vectorAndCount.count;
        _test(count==3);    	    
        _test(vectorCountingNumbers[0][0]==0);
        _test(vectorCountingNumbers[0][1]==2);
        _test(vectorCountingNumbers[0][2]==1);
        _test(vectorCountingNumbers[1][0]==0);
        _test(vectorCountingNumbers[1][1]==1);    	
        _test(vectorCountingNumbers[1][2]==0);        
        _test(vectorCountingNumbers[2][0]==2);
        _test(vectorCountingNumbers[2][1]==1);    	
        _test(vectorCountingNumbers[2][2]==2);            
    	
    	      	
 }





void TestMainEffectsConverter::
    testSliceOutOneInputVarAndOneOutputVar(){
    	
    	MainEffectsConverter x;
    	std::vector< std::vector<int> > vectorInputDataPoints(0);
    	std::vector< std::vector<double> > vectorOutputDataPoints(0);
    	int indexOfInputVariable = 0;
    	int indexOfOutputVariable = 0;
    	int numberOfValuesAvailableForOneInputVar = 0;
    	std::vector<double> oneRowOfDoubles(0);
    	std::vector<int> oneRowOfInts(0);
    	std::vector<double> row2OfDoubles(0);
    	std::vector<int> row2OfInts(0);    	
    	DDaceMainEffects::Factor factor;
    	std::vector<int> factors; 
    	DDaceMainEffects::Response responses;
    	
    	
    	/* empty data sets */
    	vectorInputDataPoints.clear();
    	vectorOutputDataPoints.clear();
    	indexOfInputVariable = 0;
    	indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 0;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 0);
    	_test(factor.getNumberOfLevels()==0);
    	_test(responses.getNumOfObservations()==0);
    	
    	/* left data set is empty */
    	vectorInputDataPoints.clear();
    	vectorOutputDataPoints.clear();
    	oneRowOfDoubles.push_back(100.0);
    	vectorOutputDataPoints.push_back(oneRowOfDoubles);
    	indexOfInputVariable = 0;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 1;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 0);
    	_test(factor.getNumberOfLevels()==0);
    	_test(responses.getNumOfObservations()==0);
    	
    	/* right data set is empty */
    	vectorInputDataPoints.clear();
    	oneRowOfInts.push_back(0);    	
    	vectorInputDataPoints.push_back(oneRowOfInts);
    	vectorOutputDataPoints.clear();
    	indexOfInputVariable = 0;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 1;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 0);
    	_test(factor.getNumberOfLevels()==0);
    	_test(responses.getNumOfObservations()==0);

    	
    	/* one data point */
    	vectorInputDataPoints.clear();
    	vectorInputDataPoints.push_back(oneRowOfInts);    	
    	vectorOutputDataPoints.clear();
    	vectorOutputDataPoints.push_back(oneRowOfDoubles);    	
    	indexOfInputVariable = 0;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 1;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 1);
    	_test(factors[0] == 0);
    	_test(factor.getNumberOfLevels()==1);
    	_test(responses.getNumOfObservations()==1);
    	_test(responses[0] == 100.0);
    	
    	/* one row, 2 columns */
    	vectorInputDataPoints.clear();
    	oneRowOfInts.clear();
    	oneRowOfInts.push_back(0);
    	oneRowOfInts.push_back(1);
    	vectorInputDataPoints.push_back(oneRowOfInts);    	
    	vectorOutputDataPoints.clear();
    	oneRowOfDoubles.clear();
    	oneRowOfDoubles.push_back(100.0);
    	oneRowOfDoubles.push_back(101.0);
    	vectorOutputDataPoints.push_back(oneRowOfDoubles);    	
    	
    	indexOfInputVariable = 0;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 2;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 1);
    	_test(factors[0] == 0);
    	_test(factor.getNumberOfLevels()==1);
    	_test(responses.getNumOfObservations()==1);
    	_test(responses[0] == 100.0);    	
    	
  	    indexOfInputVariable = 0;
        indexOfOutputVariable = 1;
    	numberOfValuesAvailableForOneInputVar = 2;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 1);
    	_test(factors[0] == 0);
    	_test(factor.getNumberOfLevels()==1);
    	_test(responses.getNumOfObservations()==1);
    	_test(responses[0] == 101.0);    	
    	
    	indexOfInputVariable = 1;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 2;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 1);
    	_test(factors[0] == 1);
    	_test(factor.getNumberOfLevels()==1);
    	_test(responses.getNumOfObservations()==1);
    	_test(responses[0] == 100.0);    	    	
    	
    	indexOfInputVariable = 1;
        indexOfOutputVariable = 1;
    	numberOfValuesAvailableForOneInputVar = 2;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 1);
    	_test(factors[0] == 1);
    	_test(factor.getNumberOfLevels()==1);
    	_test(responses.getNumOfObservations()==1);
    	_test(responses[0] == 101.0);  
    	
    	/* two rows, 1 column */
    	vectorInputDataPoints.clear();
    	oneRowOfInts.clear();
    	oneRowOfInts.push_back(0);
    	vectorInputDataPoints.push_back(oneRowOfInts);   
    	row2OfInts.clear();
    	row2OfInts.push_back(1);
    	vectorInputDataPoints.push_back(row2OfInts);   
    	vectorOutputDataPoints.clear();
    	oneRowOfDoubles.clear();
    	oneRowOfDoubles.push_back(100.0);
    	vectorOutputDataPoints.push_back(oneRowOfDoubles);    	
    	row2OfDoubles.clear();
    	row2OfDoubles.push_back(101.0);
    	vectorOutputDataPoints.push_back(row2OfDoubles);       	
    	
    	indexOfInputVariable = 0;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 2;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 2);
    	_test(factors[0] == 0);
    	_test(factors[1] == 1);
    	_test(factor.getNumberOfLevels()==2);
    	_test(responses.getNumOfObservations()==2);
    	_test(responses[0] == 100.0);    	
    	_test(responses[1] == 101.0); 

    	/* two rows, 1 column */
    	vectorInputDataPoints.clear();
    	oneRowOfInts.clear();
    	oneRowOfInts.push_back(1);
    	vectorInputDataPoints.push_back(oneRowOfInts);   
    	row2OfInts.clear();
    	row2OfInts.push_back(0);
    	vectorInputDataPoints.push_back(row2OfInts);   
    	vectorOutputDataPoints.clear();
    	oneRowOfDoubles.clear();
    	oneRowOfDoubles.push_back(101.0);
    	vectorOutputDataPoints.push_back(oneRowOfDoubles);    	
    	row2OfDoubles.clear();
    	row2OfDoubles.push_back(100.0);
    	vectorOutputDataPoints.push_back(row2OfDoubles);       	
    	
    	indexOfInputVariable = 0;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 2;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 2);
    	_test(factors[0] == 1);
    	_test(factors[1] == 0);
    	_test(factor.getNumberOfLevels()==2);
    	_test(responses.getNumOfObservations()==2);
    	_test(responses[0] == 101.0);    	
    	_test(responses[1] == 100.0);  
    	
    	/* two rows, 2 columns */
    	vectorInputDataPoints.clear();
    	oneRowOfInts.clear();
    	oneRowOfInts.push_back(0);
    	oneRowOfInts.push_back(1);
    	vectorInputDataPoints.push_back(oneRowOfInts);    	
    	row2OfInts.clear();
    	row2OfInts.push_back(2);
    	row2OfInts.push_back(3);
    	vectorInputDataPoints.push_back(row2OfInts);       	
    	vectorOutputDataPoints.clear();
    	oneRowOfDoubles.clear();
    	oneRowOfDoubles.push_back(100.0);
    	oneRowOfDoubles.push_back(101.0);
    	vectorOutputDataPoints.push_back(oneRowOfDoubles);    	
    	row2OfDoubles.clear();
    	row2OfDoubles.push_back(102.0);
        row2OfDoubles.push_back(103.0);
    	vectorOutputDataPoints.push_back(row2OfDoubles);       	
    	
    	
    	indexOfInputVariable = 0;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 4;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 2);
    	_test(factors[0] == 0);
    	_test(factors[1] == 2);
    	_test(factor.getNumberOfLevels()==2);
    	_test(responses.getNumOfObservations()==2);
    	_test(responses[0] == 100.0);    	 	
    	_test(responses[1] == 102.0);    	 	   	  	  	
    	
    	indexOfInputVariable = 0;
        indexOfOutputVariable = 1;
    	numberOfValuesAvailableForOneInputVar = 4;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 2);
    	_test(factors[0] == 0);
    	_test(factors[1] == 2);
    	_test(factor.getNumberOfLevels()==2);
    	_test(responses.getNumOfObservations()==2);
    	_test(responses[0] == 101.0);    	 	
    	_test(responses[1] == 103.0);      	
    	
    	indexOfInputVariable = 1;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 4;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 2);
    	_test(factors[0] == 1);
    	_test(factors[1] == 3);
    	_test(factor.getNumberOfLevels()==2);
    	_test(responses.getNumOfObservations()==2);
    	_test(responses[0] == 100.0);    	 	
    	_test(responses[1] == 102.0);      	
    	
    	indexOfInputVariable = 1;
        indexOfOutputVariable = 1;
    	numberOfValuesAvailableForOneInputVar = 4;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 2);
    	_test(factors[0] == 1);
    	_test(factors[1] == 3);
    	_test(factor.getNumberOfLevels()==2);
    	_test(responses.getNumOfObservations()==2);
    	_test(responses[0] == 101.0);    	 	
    	_test(responses[1] == 103.0);      
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	/* two rows, 2 columns */
    	vectorInputDataPoints.clear();
    	oneRowOfInts.clear();
    	oneRowOfInts.push_back(3);
    	oneRowOfInts.push_back(2);
    	vectorInputDataPoints.push_back(oneRowOfInts);    	
    	row2OfInts.clear();
    	row2OfInts.push_back(1);
    	row2OfInts.push_back(0);
    	vectorInputDataPoints.push_back(row2OfInts);       	
    	vectorOutputDataPoints.clear();
    	oneRowOfDoubles.clear();
    	oneRowOfDoubles.push_back(103.0);
    	oneRowOfDoubles.push_back(102.0);
    	vectorOutputDataPoints.push_back(oneRowOfDoubles);    	
    	row2OfDoubles.clear();
    	row2OfDoubles.push_back(101.0);
        row2OfDoubles.push_back(100.0);
    	vectorOutputDataPoints.push_back(row2OfDoubles);       	
    	
    	
    	indexOfInputVariable = 0;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 4;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 2);
    	_test(factors[0] == 3);
    	_test(factors[1] == 1);
    	_test(factor.getNumberOfLevels()==2);
    	_test(responses.getNumOfObservations()==2);
    	_test(responses[0] == 103.0);    	 	
    	_test(responses[1] == 101.0);    	 	   	  	  	
    	
    	indexOfInputVariable = 0;
        indexOfOutputVariable = 1;
    	numberOfValuesAvailableForOneInputVar = 4;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 2);
    	_test(factors[0] == 3);
    	_test(factors[1] == 1);
    	_test(factor.getNumberOfLevels()==2);
    	_test(responses.getNumOfObservations()==2);
    	_test(responses[0] == 102.0);    	 	
    	_test(responses[1] == 100.0);      	
    	
    	indexOfInputVariable = 1;
        indexOfOutputVariable = 0;
    	numberOfValuesAvailableForOneInputVar = 4;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 2);
    	_test(factors[0] == 2);
    	_test(factors[1] == 0);
    	_test(factor.getNumberOfLevels()==2);
    	_test(responses.getNumOfObservations()==2);
    	_test(responses[0] == 103.0);    	 	
    	_test(responses[1] == 101.0);      	
    	
    	indexOfInputVariable = 1;
        indexOfOutputVariable = 1;
    	numberOfValuesAvailableForOneInputVar = 4;
    	factor =
    	    x.sliceOutOneInputVarAndOneOutputVar
    	    (vectorInputDataPoints,
    	     vectorOutputDataPoints,
    	     indexOfInputVariable,
    	     indexOfOutputVariable,
    	     numberOfValuesAvailableForOneInputVar);
    	factors = factor.getFactors();
    	responses = factor.getResponse();
    	_test(factors.size() == 2);
    	_test(factors[0] == 2);
    	_test(factors[1] == 0);
    	_test(factor.getNumberOfLevels()==2);
    	_test(responses.getNumOfObservations()==2);
    	_test(responses[0] == 102.0);    	 	
    	_test(responses[1] == 100.0);     		
  	
}

void TestMainEffectsConverter::testConvert(){
	
	std::vector<std::vector<double> > vectorInputDataPoints(0);
	std::vector<std::vector<double> > vectorOutputDataPoints(0);
	std::vector<DDaceMainEffects::Factor>  vectorFactors(0);
	DDaceMainEffects::Factor factor;	
    std::vector<double> row1OfDoubles(0);
    std::vector<double> row2OfDoubles(0);
    std::vector<double> row3OfDoubles(0);    
    std::vector<double> row4OfDoubles(0);    
    std::vector<int> factors; 
    DDaceMainEffects::Response responses;
    
	
	
	MainEffectsConverter x;
	
	/* empty data sets */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 0);
    
    /* left side empty */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	vectorOutputDataPoints.push_back(row1OfDoubles);
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 0);

    /* right side empty */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 0);    
    
    /* one data point */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(100.0);
	vectorOutputDataPoints.push_back(row2OfDoubles);
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 1);     
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 0);
    _test(responses[0] = 100.0);
    
    /* one row, two columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	row1OfDoubles.push_back(6.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(100.0);
	row2OfDoubles.push_back(101.0);
	vectorOutputDataPoints.push_back(row2OfDoubles);
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 4);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 0);
    _test(responses[0] == 100.0);    
    
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 0);
    _test(responses[0] == 101.0);       
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 1);
    _test(responses[0] == 100.0);           
    
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 1);
    _test(responses[0] == 101.0);  
    



    /* one row, three columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	row1OfDoubles.push_back(6.0);
	row1OfDoubles.push_back(7.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(100.0);
	row2OfDoubles.push_back(101.0);
	row2OfDoubles.push_back(102.0);
	vectorOutputDataPoints.push_back(row2OfDoubles);
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 9);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 0);
    _test(responses[0] == 100.0);    
    
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 0);
    _test(responses[0] == 101.0);       
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 0);
    _test(responses[0] == 102.0);           
    
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 1);
    _test(responses[0] == 100.0);  
    
    factor = vectorFactors[4];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 1);
    _test(responses[0] == 101.0);      
    
    factor = vectorFactors[5];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 1);
    _test(responses[0] == 102.0);      
    
    factor = vectorFactors[6];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 2);
    _test(responses[0] == 100.0);      
    
    factor = vectorFactors[7];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 2);
    _test(responses[0] == 101.0);          
    
    factor = vectorFactors[8];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 1);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==1);
    _test(factors[0] == 2);
    _test(responses[0] == 102.0);          
    
    /* two rows, one column */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(6.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(100.0);
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(101.0);
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 1);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 1);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 101.0);    
    
    /* two rows, one column */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(6.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(5.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(101.0);
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(100.0);
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 1);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 0);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 100.0);     
    
    
    /* two rows, one column */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(5.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(101.0);
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(100.0);
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 1);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 0);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 100.0);         
    
    /* two rows, two columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	row1OfDoubles.push_back(6.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(7.0);
	row2OfDoubles.push_back(8.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(100.0);
	row3OfDoubles.push_back(101.0);	
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(102.0);
	row4OfDoubles.push_back(103.0);	
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 4);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 2);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 2);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);                 	
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 3);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 3);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);      
    
    /* two rows, two columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(8.0);
	row1OfDoubles.push_back(7.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(6.0);
	row2OfDoubles.push_back(5.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(103.0);
	row3OfDoubles.push_back(102.0);	
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(101.0);
	row4OfDoubles.push_back(100.0);	
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 4);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 3);
    _test(factors[1] == 1);    
    _test(responses[0] == 103.0);
    _test(responses[1] == 101.0);       
    
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 3);
    _test(factors[1] == 1);    
    _test(responses[0] == 102.0);
    _test(responses[1] == 100.0);                 	
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 2);
    _test(factors[1] == 0);    
    _test(responses[0] == 103.0);
    _test(responses[1] == 101.0);       
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 2);
    _test(factors[1] == 0);    
    _test(responses[0] == 102.0);
    _test(responses[1] == 100.0);   
    
    
    /* two rows, two columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	row1OfDoubles.push_back(5.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(7.0);
	row2OfDoubles.push_back(8.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(100.0);
	row3OfDoubles.push_back(101.0);	
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(102.0);
	row4OfDoubles.push_back(103.0);	
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 4);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 1);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 1);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);                 	
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 2);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 2);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);   
    

    /* two rows, two columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	row1OfDoubles.push_back(6.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(5.0);
	row2OfDoubles.push_back(8.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(100.0);
	row3OfDoubles.push_back(101.0);	
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(102.0);
	row4OfDoubles.push_back(103.0);	
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 4);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 0);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 0);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);                 	
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 2);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 2);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);     
    
    
    
    /* two rows, two columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	row1OfDoubles.push_back(6.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(7.0);
	row2OfDoubles.push_back(5.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(100.0);
	row3OfDoubles.push_back(101.0);	
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(102.0);
	row4OfDoubles.push_back(103.0);	
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 4);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 2);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);   
    
        
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    
    _test(factors[0] == 0);
    _test(factors[1] == 2);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);                	
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 0);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 0);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);     
     
    /* two rows, two columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	row1OfDoubles.push_back(6.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(6.0);
	row2OfDoubles.push_back(8.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(100.0);
	row3OfDoubles.push_back(101.0);	
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(102.0);
	row4OfDoubles.push_back(103.0);	
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 4);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 1);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 1);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);                 	
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 2);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 2);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);      
    
    
    /* two rows, two columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	row1OfDoubles.push_back(6.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(7.0);
	row2OfDoubles.push_back(6.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(100.0);
	row3OfDoubles.push_back(101.0);	
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(102.0);
	row4OfDoubles.push_back(103.0);	
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 4);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 2);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 2);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);                 	
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 1);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==1);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 1);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0); 
    
    
    /* two rows, two columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(5.0);
	row1OfDoubles.push_back(6.0);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(7.0);
	row2OfDoubles.push_back(7.0);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(100.0);
	row3OfDoubles.push_back(101.0);	
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(102.0);
	row4OfDoubles.push_back(103.0);	
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 4);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 2);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 2);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);                 	
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 2);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 2);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0); 
    
    
    /* two rows, two columns */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row1OfDoubles.push_back(4.9999);
	row1OfDoubles.push_back(6.0001);
	vectorInputDataPoints.push_back(row1OfDoubles);
	row2OfDoubles.clear();
	row2OfDoubles.push_back(5.9999);
	row2OfDoubles.push_back(5.0001);
	vectorInputDataPoints.push_back(row2OfDoubles);
	row3OfDoubles.clear();
	row3OfDoubles.push_back(100.0);
	row3OfDoubles.push_back(101.0);	
	vectorOutputDataPoints.push_back(row3OfDoubles);
	row4OfDoubles.clear();
	row4OfDoubles.push_back(102.0);
	row4OfDoubles.push_back(103.0);	
	vectorOutputDataPoints.push_back(row4OfDoubles);	
	vectorFactors = x.convert(vectorInputDataPoints, vectorOutputDataPoints);
    _test(vectorFactors.size() == 4);     
    
    factor = vectorFactors[0];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 1);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[1];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 0);
    _test(factors[1] == 1);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);                 	
    
    factor = vectorFactors[2];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 0);    
    _test(responses[0] == 100.0);
    _test(responses[1] == 102.0);       
    
    factor = vectorFactors[3];
    factors = factor.getFactors();
    responses = factor.getResponse();
    _test(factors.size() == 2);   
    _test(factor.getNumberOfLevels()==2);
    _test(responses.getNumOfObservations()==2);
    _test(factors[0] == 1);
    _test(factors[1] == 0);    
    _test(responses[0] == 101.0);
    _test(responses[1] == 103.0);      
         
         
    
                    	
    
}


void TestMainEffectsConverter::testConstructor(){
}






