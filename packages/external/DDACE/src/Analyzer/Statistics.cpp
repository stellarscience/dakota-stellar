#include "Statistics.h"

double Statistics::sum(std::vector<double> vec)
{
	double ans = 0;
	for(int i = 0; i < (int) vec.size(); i++)
	{
		ans += vec[i];
	}
	return ans;
}

double Statistics::average(std::vector<double> vec)
{
	//BUG FIX:
	if (vec.size()==0) return(0.0);
	
	
	return (sum(vec) / vec.size());
}

double Statistics::sumOfSquares(std::vector<double> vec, double num)
{
	double sum = 0;
	for(int i = 0; i < (int) vec.size(); i++)
	{
		sum += pow((vec[i]-num),2);
	}
	return sum;
}

double Statistics::variance(std::vector<double> vec)
{
	//BUG FIX:
	if (vec.size() <= 1) return(0.0);
	
	
	double mean = average(vec);
	return sumOfSquares(vec,mean) / (vec.size() - 1);
}

double Statistics::stdDeviation(std::vector<double> vec)
{
	return sqrt(variance(vec));
} 
