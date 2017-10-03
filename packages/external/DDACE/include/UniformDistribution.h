#ifndef UNIFORMDIST_H
#define UNIFORMDIST_H

#include "Distribution.h"
#include <vector>
#include <stdexcept>
#include <cmath>


	/**\ingroup Random 
	 * Generates random numbers from a uniform distribution.
	 */

	class UniformDistribution : public DistributionBase
		{
		public:
			/** create a uniform distribution on the unit interval [0,1] */
			UniformDistribution();
  
			/** create a uniform distribution on the interval [lowerBound, upperBound] */
			UniformDistribution(double lowerBound, double upperBound);
  
			virtual ~UniformDistribution(){;}
  
			virtual DistributionBase* clone() const;
  
			/** draws a uniformly-distributed deviate */
			virtual double getDeviate() const;

			/** find x such that CDF(x) = prob */
			virtual double getDeviate(double prob) const;

			/** evaluate the CDF at x */
			virtual double getCDF(double x) const ;
    
			/** find the lower bound of the distribution */
			virtual double lowerBound() const;

			/** find the upper bound of the distribution */
			virtual double upperBound() const;

			/** find the mean of the distribution */
			virtual double mean() const;

			/** find the standard deviation of the distribution */
			virtual double stdDev() const;

			/** print a description */
			virtual void print(std::ostream& os) const;
			virtual void printAttributes(std::ostream& os) const;
  
			/** return the name "UniformDistribution" */
			virtual const std::string& typeName() const {return typeName_;}
		protected:
			double lowerBound_;
			double upperBound_;
			static std::string typeName_;
		};

	std::vector<int> randomIVector(const int leng);
	std::vector<int> randomIVector(const int leng,
					 const UniformDistribution& dev);
	std::vector<int> randomIVector(const int leng, const int seed);



#endif // UNIFORMDIST_H


