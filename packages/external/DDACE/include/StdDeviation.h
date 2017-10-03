#ifndef STDDEV_H
#define STDDEV_H

#include <cmath>
#include "Mean.h"


	/**\ingroup Random 
	 * Wrap a double in a StdDeviation object, to ensure proper ordering
	 * of arguments to distribution constructors. 
	 * @author Leslea Lehoucq
	 */

	class StdDeviation
		{
		public:
			/** create a standard deviation object */
			StdDeviation(double val) : stdDev_(val) {;}
			/** find the standard deviation of a data set */
			StdDeviation(const std::vector<double>& data);

			/** get the double value of a StandardDeviation object */
			double value() const;

		protected:
			double stdDev_;
		};

	inline StdDeviation::StdDeviation(const std::vector<double>& data)
		{
			Mean m(data);
			double tmpSum = 0;
			for(int i = 0; i < (int) data.size(); i++)
				tmpSum += pow((data[i] - m.value()), 2);
			stdDev_ = sqrt(tmpSum/(data.size() - 1));
		}

	inline double StdDeviation::value() const
		{
			return stdDev_;
		}

#endif


