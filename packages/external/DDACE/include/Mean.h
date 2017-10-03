#ifndef MEAN_H
#define MEAN_H

#include <vector>


	/**\ingroup Random 
 * Wrap a double in a Mean object, to ensure proper ordering of
 * arguments to Distribution ctors.
 */
	class Mean
		{
		public:
			/** create a Mean object with value zero */
			Mean() : mean_(0) {;}

			/** create a Mean object */
			Mean(double val) : mean_(val) {;}

			/** find the mean of a data set */
			Mean(const std::vector<double>& data);

			/** return the double value of this Mean object */
			double value() const;

		protected:
			double mean_;
		};

	inline Mean::Mean(const std::vector<double>& data)
		{
			double tmpSum = 0;
			for(int i = 0; i < (int) data.size(); i++)
				tmpSum += data[i];
			mean_ = tmpSum/(data.size());
		}

	inline double Mean::value() const
		{
			return mean_;
		}
#endif
