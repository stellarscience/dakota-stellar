#include "OneWayANOVA.h"
#include <iostream>

int main(int argc, char** argv)
{
	vector<double> resp;
	resp.push_back(11.06);
	resp.push_back(15.4);
	resp.push_back(11.81);
	resp.push_back(20.48);
	resp.push_back(20.32);
	resp.push_back(10.96);
	resp.push_back(15.86);
	resp.push_back(10.8);
	resp.push_back(19.73);
	resp.push_back(20.57);

	Response response(resp);

	vector<int> facs;
	facs.push_back(0);
	facs.push_back(0);
	facs.push_back(0);
	facs.push_back(1);
	facs.push_back(1);
	facs.push_back(0);
	facs.push_back(0);
	facs.push_back(0);
	facs.push_back(1);
	facs.push_back(1);

	Factor factor(facs,2,response);
	vector<Factor> factors;
	
	factors.push_back(factor);
/**
	std::cout << "Variation Source    SS     DoF     Variance    Fdata" << std::endl;
	std::cout << "Between Groups: " << factor.sumOfSquaresBetweenGroups()
			<< "   " << factor.doFBetween() << "   "
			<< factor.varianceBetweenGroups() << "   " << factor.Fdata()
			<< std::endl;
	std::cout << "Within Groups : " << factor.sumOfSquaresWithinGroups()
			<< "   " << factor.doFWithin() << "   " 
			<< factor.varianceWithinGroups() << "    " << std::endl;
*/
	OneWayANOVA screwTest(factors,response);

	screwTest.printANOVATable(0);
	
	return 0;
}
