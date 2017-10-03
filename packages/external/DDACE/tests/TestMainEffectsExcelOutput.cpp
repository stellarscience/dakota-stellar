#include "TestMainEffectsExcelOutput.h"
#include <cstdlib>
#include <cstdio>
std::ostringstream ss3;
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <limits>

using namespace std;

TestMainEffectsExcelOutput::TestMainEffectsExcelOutput(){ }

TestMainEffectsExcelOutput::~TestMainEffectsExcelOutput(){ }

void TestMainEffectsExcelOutput::run() {
	testComputeExcelOutput();

}



void TestMainEffectsExcelOutput::
    testComputeExcelOutput() {
	
	
	std::vector<std::vector<double> > vectorInputDataPoints(0);
	std::vector<std::vector<double> > vectorOutputDataPoints(0);
	std::vector<DDaceMainEffects::Factor>  vectorFactors(0);
	DDaceMainEffects::Factor factor;	
    std::vector<double> row1OfDoubles(0);
    std::vector<double> row2OfDoubles(0);
    std::vector<double> row3OfDoubles(0);    
    std::vector<double> row4OfDoubles(0);    
    std::vector<double> row5OfDoubles(0);
    std::vector<double> row6OfDoubles(0);
    std::vector<double> row7OfDoubles(0);    
    std::vector<double> row8OfDoubles(0);    
    std::vector<double> row9OfDoubles(0);
    std::vector<double> row10OfDoubles(0);
    std::vector<double> row11OfDoubles(0);    
    std::vector<double> row12OfDoubles(0);    
    std::vector<double> row13OfDoubles(0);
    std::vector<double> row14OfDoubles(0);
    std::vector<double> row15OfDoubles(0);    
    std::vector<double> row16OfDoubles(0);    
    std::vector<double> row17OfDoubles(0);
    std::vector<double> row18OfDoubles(0);    
    std::vector<double> row19OfDoubles(0);    
    std::vector<double> row20OfDoubles(0);
    std::string s;
    
    std::vector<int> factors; 
    DDaceMainEffects::Response responses;
    
	
	
	MainEffectsExcelOutput x;
	
    /* The NaNs are throwing floating point exceptions */
//	/* empty data sets */
//	vectorInputDataPoints.clear();
//	vectorOutputDataPoints.clear();
//	s = x.computeExcelOutput
//	    (vectorInputDataPoints, 
//	     vectorOutputDataPoints);
//	_test(s == "");
//	
//	/* one data point */
//	vectorInputDataPoints.clear();
//	vectorOutputDataPoints.clear();
//	row1OfDoubles.clear();
//	row2OfDoubles.clear();
//	row1OfDoubles.push_back(5.0);
//	vectorInputDataPoints.push_back(row1OfDoubles);
//	row2OfDoubles.push_back(100.0);
//	vectorOutputDataPoints.push_back(row2OfDoubles);	
//	s = x.computeExcelOutput
//	    (vectorInputDataPoints, 
//	     vectorOutputDataPoints);
//	          
//    
//    ss1 << "in(0),out(0),nObservations,sumOfAllObservations,";
//    ss1 << "avgOfAllObservation,sumOfSquaresOfAllObservations,";
//    ss1 << "degreesOfFreedomOfAllObservations,varianceOfAllObservations,";
//    ss1 << "sum,average,sumOfSquares,variance,sumOfSquaresBetweenGroups,";
//    ss1 << "degreesOfFreedomBetweenGroups,varianceBetweenGroups,";
//    ss1 << "sumOfSquaresWithinGroups,degreesOfFreedomWithinGroups,";
//    ss1 << "varianceWithinGroups,F";
//    ss1 << "\n";
//    ss1 << "F,R,1,100,100,0,0,0,100,100,0,0,0,0,nan,0,0,nan,nan";
//    ss1 << "\n";    
//    _test(s==ss1.str()); //NO NO NO!!! Not everybody handles nan the same way
//    ss1.str("");

	
	
	
    /* The NaNs are throwing floating point exceptions */	
//	/* one row, 2 columns */
//	vectorInputDataPoints.clear();
//	vectorOutputDataPoints.clear();
//	row1OfDoubles.clear();
//	row2OfDoubles.clear();
//	row1OfDoubles.push_back(5.0);
//	row1OfDoubles.push_back(6.0);
//	vectorInputDataPoints.push_back(row1OfDoubles);
//	row2OfDoubles.push_back(100.0);
//	row2OfDoubles.push_back(101.0);	
//	vectorOutputDataPoints.push_back(row2OfDoubles);	
//	s = x.computeExcelOutput
//	    (vectorInputDataPoints, 
//	     vectorOutputDataPoints);
//	     
//	 
//                                   
//    ss2 << "in(0),in(1),out(0),out(1),nObservations,sumOfAllObservations,";
//    ss2 << "avgOfAllObservation,sumOfSquaresOfAllObservations,";
//    ss2 << "degreesOfFreedomOfAllObservations,varianceOfAllObservations,";
//    ss2 << "sum,average,sumOfSquares,variance,sumOfSquaresBetweenGroups,";
//    ss2 << "degreesOfFreedomBetweenGroups,varianceBetweenGroups,";
//    ss2 << "sumOfSquaresWithinGroups,degreesOfFreedomWithinGroups,";
//    ss2 << "varianceWithinGroups,F";
//    ss2 << "\n";
//    ss2 << "F,,R,,1,100,100,0,0,0,100,100,0,0,0,0,nan,0,0,nan,nan";
//    ss2 << "\n";
//    ss2 << "F,,,R,1,101,101,0,0,0,101,101,0,0,0,0,nan,0,0,nan,nan";
//    ss2 << "\n";    
//    ss2 << "F,R,,1,100,100,0,0,0,100,100,0,0,0,0,nan,0,0,nan,nan";
//    ss2 << "\n";    
//    ss2 << "F,,R,1,101,101,0,0,0,101,101,0,0,0,0,nan,0,0,nan,nan";
//    ss2 << "\n";    
//    
//    //_test(s==ss2.str());  //NO NO NO!!! Not everybody handles nan the same way
//    ss2.str("");
	
	
	/* 10 rows, 1 column */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row2OfDoubles.clear();	
	row3OfDoubles.clear();
	row4OfDoubles.clear();		
	row5OfDoubles.clear();
	row6OfDoubles.clear();		
	row7OfDoubles.clear();
	row8OfDoubles.clear();	
	row9OfDoubles.clear();
	row10OfDoubles.clear();		
	row11OfDoubles.clear();
	row12OfDoubles.clear();	
	row13OfDoubles.clear();
	row14OfDoubles.clear();		
	row15OfDoubles.clear();
	row16OfDoubles.clear();		
	row17OfDoubles.clear();
	row18OfDoubles.clear();	
	row19OfDoubles.clear();
	row20OfDoubles.clear();
	

    row1OfDoubles.push_back(0);  row11OfDoubles.push_back(11.06);
    row2OfDoubles.push_back(0);  row12OfDoubles.push_back(15.4);
    row3OfDoubles.push_back(0);  row13OfDoubles.push_back(11.81);
    row4OfDoubles.push_back(1);  row14OfDoubles.push_back(20.48);
    row5OfDoubles.push_back(1);  row15OfDoubles.push_back(20.32);    
    row6OfDoubles.push_back(0);  row16OfDoubles.push_back(10.96);
    row7OfDoubles.push_back(0);  row17OfDoubles.push_back(15.86);
    row8OfDoubles.push_back(0);  row18OfDoubles.push_back(10.8);
    row9OfDoubles.push_back(1);  row19OfDoubles.push_back(19.73);
    row10OfDoubles.push_back(1); row20OfDoubles.push_back(20.57);	
    
	vectorInputDataPoints.push_back(row1OfDoubles);
	vectorInputDataPoints.push_back(row2OfDoubles);
	vectorInputDataPoints.push_back(row3OfDoubles);
	vectorInputDataPoints.push_back(row4OfDoubles);
	vectorInputDataPoints.push_back(row5OfDoubles);			
	vectorInputDataPoints.push_back(row6OfDoubles);
	vectorInputDataPoints.push_back(row7OfDoubles);
	vectorInputDataPoints.push_back(row8OfDoubles);
	vectorInputDataPoints.push_back(row9OfDoubles);
	vectorInputDataPoints.push_back(row10OfDoubles);	
	
	vectorOutputDataPoints.push_back(row11OfDoubles);
	vectorOutputDataPoints.push_back(row12OfDoubles);
	vectorOutputDataPoints.push_back(row13OfDoubles);
	vectorOutputDataPoints.push_back(row14OfDoubles);
	vectorOutputDataPoints.push_back(row15OfDoubles);				
	vectorOutputDataPoints.push_back(row16OfDoubles);
	vectorOutputDataPoints.push_back(row17OfDoubles);
	vectorOutputDataPoints.push_back(row18OfDoubles);
	vectorOutputDataPoints.push_back(row19OfDoubles);
	vectorOutputDataPoints.push_back(row20OfDoubles);	
	
	s = x.computeExcelOutput (vectorInputDataPoints, 
	     vectorOutputDataPoints);	     


    ss3 << "in(0),out(0),nObservations,sumOfAllObservations,";
    ss3 << "avgOfAllObservation,sumOfSquaresOfAllObservations,";
    ss3 << "degreesOfFreedomOfAllObservations,varianceOfAllObservations,";
    ss3 << "sum,average,sumOfSquares,variance,sumOfSquaresBetweenGroups,";
    ss3 << "degreesOfFreedomBetweenGroups,varianceBetweenGroups,";
    ss3 << "sumOfSquaresWithinGroups,degreesOfFreedomWithinGroups,";
    ss3 << "varianceWithinGroups,F";
    ss3 << "\n";
    ss3 << "F,R,10,156.99,15.699,167.405,9,18.6006,75.89,12.6483,27.3789,5.47578,139.599,1,139.599,27.807,8,3.47587,40.1621";
    ss3 << "\n";
    ss3 << "F,R,,,,,,,81.1,20.275,0.4281,0.1427,,,,,,,";
    ss3 << "\n";    

    cout << ss3.str() << endl;
    
    _test(s==ss3.str());        
    //ss3.str("");    
    
    
    
	/* 1 column of input, 2 columns of output */
	vectorInputDataPoints.clear();
	vectorOutputDataPoints.clear();
	row1OfDoubles.clear();
	row2OfDoubles.clear();	
	row3OfDoubles.clear();
	row4OfDoubles.clear();		
	row5OfDoubles.clear();
	row6OfDoubles.clear();		
	row7OfDoubles.clear();
	row8OfDoubles.clear();	
	row9OfDoubles.clear();
	row10OfDoubles.clear();		
	row11OfDoubles.clear();
	row12OfDoubles.clear();	
	row13OfDoubles.clear();
	row14OfDoubles.clear();		
	row15OfDoubles.clear();
	row16OfDoubles.clear();		
	row17OfDoubles.clear();
	row18OfDoubles.clear();	
	row19OfDoubles.clear();
	row20OfDoubles.clear();
	

    row1OfDoubles.push_back(0);  row11OfDoubles.push_back(11.06);
    row2OfDoubles.push_back(0);  row12OfDoubles.push_back(15.4);
    row3OfDoubles.push_back(0);  row13OfDoubles.push_back(11.81);
    row4OfDoubles.push_back(1);  row14OfDoubles.push_back(20.48);
    row5OfDoubles.push_back(1);  row15OfDoubles.push_back(20.32);    
    row6OfDoubles.push_back(0);  row16OfDoubles.push_back(10.96);
    row7OfDoubles.push_back(0);  row17OfDoubles.push_back(15.86);
    row8OfDoubles.push_back(0);  row18OfDoubles.push_back(10.8);
    row9OfDoubles.push_back(1);  row19OfDoubles.push_back(19.73);
    row10OfDoubles.push_back(1); row20OfDoubles.push_back(20.57);	
    

    row11OfDoubles.push_back(158.9);
    row12OfDoubles.push_back(191.58);
    row13OfDoubles.push_back(158.7);
    row14OfDoubles.push_back(166.48);
    row15OfDoubles.push_back(165.28);    
    row16OfDoubles.push_back(157.02);
    row17OfDoubles.push_back(191.24);
    row18OfDoubles.push_back(152.09);
    row19OfDoubles.push_back(165.5);
    row20OfDoubles.push_back(167.43);	
    
    
	vectorInputDataPoints.push_back(row1OfDoubles);
	vectorInputDataPoints.push_back(row2OfDoubles);
	vectorInputDataPoints.push_back(row3OfDoubles);
	vectorInputDataPoints.push_back(row4OfDoubles);
	vectorInputDataPoints.push_back(row5OfDoubles);			
	vectorInputDataPoints.push_back(row6OfDoubles);
	vectorInputDataPoints.push_back(row7OfDoubles);
	vectorInputDataPoints.push_back(row8OfDoubles);
	vectorInputDataPoints.push_back(row9OfDoubles);
	vectorInputDataPoints.push_back(row10OfDoubles);	
	
	vectorOutputDataPoints.push_back(row11OfDoubles);
	vectorOutputDataPoints.push_back(row12OfDoubles);
	vectorOutputDataPoints.push_back(row13OfDoubles);
	vectorOutputDataPoints.push_back(row14OfDoubles);
	vectorOutputDataPoints.push_back(row15OfDoubles);				
	vectorOutputDataPoints.push_back(row16OfDoubles);
	vectorOutputDataPoints.push_back(row17OfDoubles);
	vectorOutputDataPoints.push_back(row18OfDoubles);
	vectorOutputDataPoints.push_back(row19OfDoubles);
	vectorOutputDataPoints.push_back(row20OfDoubles);	
	
	s = x.computeExcelOutput (vectorInputDataPoints, 
	     vectorOutputDataPoints);	     

    std::ostringstream ss4;

    ss4 << "in(0),out(0),out(1),nObservations,sumOfAllObservations,";
    ss4 << "avgOfAllObservation,sumOfSquaresOfAllObservations,";
    ss4 << "degreesOfFreedomOfAllObservations,varianceOfAllObservations,";
    ss4 << "sum,average,sumOfSquares,variance,sumOfSquaresBetweenGroups,";
    ss4 << "degreesOfFreedomBetweenGroups,varianceBetweenGroups,";
    ss4 << "sumOfSquaresWithinGroups,degreesOfFreedomWithinGroups,";
    ss4 << "varianceWithinGroups,F";
    ss4 << "\n";
    ss4 << "F,R,,10,156.99,15.699,167.405,9,18.6006,75.89,12.6483,27.3789,5.47578,139.599,1,139.599,27.807,8,3.47587,40.1621";
    ss4 << "\n";
    ss4 << "F,R,,,,,,,,81.1,20.275,0.4281,0.1427,,,,,,,";
    ss4 << "\n";
    ss4 << "F,,R,10,1674.22,167.422,1652.05,9,183.561,1009.53,168.255,1638.71,327.742,10.4083,1,10.4083,1641.64,8,205.205,0.0507217";
    ss4 << "\n";
    ss4 << "F,,R,,,,,,,664.69,166.173,2.92468,0.974892,,,,,,,";
    ss4 << "\n";    
    
    
    _test(s==ss4.str());        
    //ss4.str("");      
    
	
}





