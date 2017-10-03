#ifndef NORMALDIST_H
#define NORMALDIST_H

#include "Distribution.h"
#include "Mean.h"
#include "StdDeviation.h"
#include <stdexcept>
#include <cmath>


	/**\ingroup Random 
	 * Generates random numbers from a normal (Gaussian) distribution.
	 */
	class NormalDistribution : public DistributionBase
		{
		public:
	/** Normal(m, sigma^2) constructor, assuming sampling range is
	 * within numDeviations of the mean, m (sampling range
	 * is 2*numDeviations wide). */
			NormalDistribution(const Mean& m, 
					   const StdDeviation& sigma, 
  					   double numDeviations); 

	/** Normal(m, sigma^2) constructor, assuming sampling range is
	 * within two standard deviations of the mean, m (sampling range
	 * is four standard deviations wide). */
			NormalDistribution(const Mean& m, 
					   const StdDeviation& sigma);

	/** Normal(m, sigma^2) constructor; we use lower and upper to 
	 * calculate m, sigma.  Use the assumption that .95 of the
	 * distrbition lies between lower and upper, and that the
	 * mean, m, lies mid-way between lower and upper:
	 *      m = lower + (upper - lower)/2
	 * Using standard Normal (Norm(0,1)) theory, and the transformation
	 * to a standard Normal using Z = (X - m)/sigma where Z~Normal(0,1),
 	 * and X~Normal(m, sigma^2), we can calculate sigma. */
			NormalDistribution(double lower, double upper);
	
	/** Calculate m and sigma from lower and upper, with the assumption
	 * that the distance between lower and upper
	 * is  2*numDeviations.  Therefore 
	 *      m = lower + (upper - lower)/2
	 *  sigma = (upper - lower)/(2*numDeviations). */
			NormalDistribution(double lower, double upper,
					   double numDeviations);

  

			virtual ~NormalDistribution(){;}
  
			virtual DistributionBase* clone() const;

			/** get a normally-distributed deviate */
			virtual double getDeviate() const; 

			/** find x such that CDF(x) = prob */
			virtual double getDeviate(double prob) const; 

			/** evaluate the CDF at x */
			virtual double getCDF(double x) const ;

			/** find the lower bound of the distribution */
			virtual double lowerBound() const { return lower_; }

			/** find the upper bound of the distribution */
			virtual double upperBound() const { return upper_; }

			/** find the mean of the distribution */
			virtual double mean() const { return mean_; }

			/** find the standard deviation of the distribution */
			virtual double  stdDev() const { return stdDev_; }

			/** print a description */
			virtual void print(std::ostream& os) const;
			virtual void printAttributes(std::ostream& os) const;
  
			/** return the name "NormalDistribution" */
			virtual const std::string& typeName() const {return typeName_;}
		protected:
			double getUntruncatedCDF(double x) const ;
			double mean_;
			double stdDev_;
			// Mean mean_;
			// StdDeviation stdDev_;
			double lower_;  // assume these will be needed....
					// Can calculate them

			double upper_;  // from input from any of the 
					// constructors above.
					// normal CDF evaluated at 
					// truncation points
			double pLower_; 
			double pUpper_; 
			static std::string typeName_ ;
		};
#endif




