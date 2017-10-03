#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
#include <iostream>
#include <cmath>


class Statistics
{
  public:
	static double sum(std::vector<double> data);
	static double average(std::vector<double> data);
	static double sumOfSquares(std::vector<double> data, double num);
	static double variance(std::vector<double> data);
	static double stdDeviation(std::vector<double> data);
};

#endif
