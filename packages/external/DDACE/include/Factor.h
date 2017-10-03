#ifndef FACTOR_H
#define FACTOR_H

#include <vector>
#include <stdexcept>
#include "SmartPtr.h"
#include "Response.h"

namespace DDaceMainEffects {

class Factor
{
  public:
  
    Factor();


	/**
	 * This is the constructor for a Factor object.
	 * A "Factor" is 2 parallel 1-D vectors.  One of the vectors contains
	 * the values for one input data variable.  The other vector contains
	 * the values for one output data variable.  The two vectors are 
	 * parallel; whenever we do a run, or experiment, we append one value
	 * to the input vector and we append one value to the output vector.
	 * @param vector<int> factors: The values that we acquired for each
	 * and every run from one input data variable (or factor).  The first
	 * element contains the value that we acquired during the first run.
	 * The value must be an integer.
	 * The value must be equal to or greater than zero.
	 * The value must be less than nLevels.
	 * This vector must parallel the response vector.
	 * @param int nLevels: the number of different values that are
	 * available for an input variable.  The numberOfLevels may or may
	 * not be equal to the number of distinct input values.  For example,
	 * if an input variable can take one of the following values
	 * {0,1,2,3,4,5,...,1000} then nLevels is 1001.  If we have 3 input
	 * variables and if the values of the 3 input variables are {0,1,1}
	 * then the number of distinct values is 2.
	 * data points.  If the nLevels is 5, an input data point (or factor)
	 * can be set to either 0, 1, 2, 3, or 4.  nLevels MUST be larger
	 * then that biggest value in factors.
	 * @param resp: The values that we acquired, or computed, for each
	 * and every run for one output data variable (or response).  The first
	 * element contains the value that we acquired, or computed,
	 * during the first run.
	 */
	Factor(std::vector<int> factors, int nLevels, Response response);
	
	
	/**
	 * Get the factors.
	 * This is a 1-D vector of values for one input data variable.
	 * Each time we do a run or experiment, a new value is appended
	 * to this vector.
	 */
	std::vector<int> getFactors();
	
    /**
     * Get the response.
     * This is a 1-D vector of values for one output data variable.
     * Each time we a run or experiment, a new value is appended
     * to this vector.
     */
    Response getResponse();
	
	// calculates the sum of squares within groups as the summation of
	// the sum of squares of each different level.
	double sumOfSquaresWithinGroups();
	
	
	
	// calculates the sum of squares between groups as the difference
	// between the total sum of squares of the population and the sum
	// of squares within groups
	double sumOfSquaresBetweenGroups();
	

	// calculates the variance between groups by dividing the sum of
	// squares between groups by the degrees of freedom between groups
	double varianceBetweenGroups();
	

	// calculates the variance within groups by dividing the sum of
	// squares within groups by the degrees of freedom within groups
	double varianceWithinGroups();
	

	// calculates the Fdata by dividing the variance between sources
	// by the variance within sources.
	/**
	 *	double Fdata()
	 *	returns the number of observations in the experiment.
	**/
	double Fdata();


	/**
	 * The number of distinct input values.
	 * <p>
	 * EXAMPLE:
	 * If factor contains the following input values
	 * {0,0,1,1,1,555} then the number of distinct values is 3.
	 * <p>
	 * The number of distinct input values may or may not
	 * be equal to the number of levels.  The number of levels
	 * is the number of values that are available for an input
	 * variable.  For example, if an input variable can be set to
	 * any of the following values {0,1,2,3,4,5,6,7,...,1000} then
	 * the number of levels is 1001.
	**/
	int getNumberOfLevels() const;

	/**
	 *	int getNumberOfObservations()
	 *	returns the number of observations in the experiment.
	**/
	int getNumberOfObservations() const;

	/**
	 *	int doFWithin()
	 *	returns the number of degrees of freedom within sources
	**/
	int doFWithin() const;

	/**
	 *	int doFBetween()
	 *	returns the number of degrees of freedom between sources
	**/
	int doFBetween() const;

	/**
	 *	Computes the sum of all output values for a single level.
	 *  <p>
	 *  The input parameter is NOT one of the values in the input vector.
	 *  The input parameter is an INDEX.  That is, after duplicates are
	 *  removed, the input data values are sorted.  The smallest distinct
	 *  input data value has an index of 0.  
	 * <p>
	 * EXAMPLE: <br>
	 * Suppose our (input,output) data values look like this: <br>
	 * {(5,1), (15,2), (10,3), (5,3), (10,5)} <br>
	 * We have 3 distinct input values: {5,15,10} <br>
	 * The input values are sorted:  {5,10,15} <br>
	 * The smallest value (5) is mapped to the index value 0. <br>
	 * The 2nd highest value (10) is mapped to the index value 1. <br>
	 * The 3rd highest value (15) is mapped to the index value 2. <br>
	 * The sum of the values for index 0 is (1+3) = 4. <br>
	 * The sum of the y values for index 1 is (3+5) = 8. <br>
	 * The sum of the y values for index 2 is 2. <br>
	 *	@param int indexOfLevels: indentifies one of the levels.  
	 *  The lowest value in the input data has an index value of 0.
	 * @return The sum of all output values for data values that have
	 * a specific input value.
	 */
	double getLevelSum(int indexOfLevel);

	/**
	 *	Computes the average output value of a single level.
	 *  <p>
	 *  The input parameter is NOT one of the values in the input vector.
	 *  The input parameter is an INDEX.  That is, after duplicates are
	 *  removed, the input data values are sorted.  The smallest distinct
	 *  input data value has an index of 0.  
	 * <p>
	 * EXAMPLE: <br>
	 * Suppose our (input,output) data values look like this: <br>
	 * {(5,1), (15,2), (10,3), (5,3), (10,5)} <br>
	 * We have 3 distinct input values: {5,15,10} <br>
	 * The input values are sorted:  {5,10,15} <br>
	 * The smallest value (5) is mapped to the index value 0. <br>
	 * The 2nd highest value (10) is mapped to the index value 1. <br>
	 * The 3rd highest value (15) is mapped to the index value 2. <br>
	 * The average y value for index 0 is (1+3)/2 = 2. <br>
	 * The average y value for index 1 is (3+5)/2 = 4. <br>
	 * The average y value for index 2 is 2. <br>
	 *	@param int indexOfLevels: indentifies one of the levels.  
	 *  The lowest value in the input data has an index value of 0.
	 * @return The average output value for all data values that have
	 * a specific input value.
	 */
	double getLevelAverage(int indexOfLevel);
	
	
	/**
	 *	Computes the sum of squares of a single level.
	 *  <p>
	 *  The input parameter is NOT one of the values in the input vector.
	 *  The input parameter is an INDEX.  That is, after duplicates are
	 *  removed, the input data values are sorted.  The smallest distinct
	 *  input data value has an index of 0.  
	 * <p>
	 * EXAMPLE: <br>
	 * Suppose our (input,output) data values look like this: <br>
	 * {(5,1), (15,2), (10,3), (5,3), (10,5)} <br>
	 * We have 3 distinct input values: {5,15,10} <br>
	 * The input values are sorted:  {5,10,15} <br>
	 * The smallest value (5) is mapped to the index value 0. <br>
	 * The 2nd highest value (10) is mapped to the index value 1. <br>
	 * The 3rd highest value (15) is mapped to the index value 2. <br>
	 * Compute the sum of squares for the y values {1,3} at index 0. <br>
	 * Compute the sum of squares for the y values {3,5} at index 1. <br>
	 * Compute the sum of squares for the y values {2} at index 2. <br>
	 *	@param int indexOfLevels: indentifies one of the levels.  
	 *  The lowest value in the input data has an index value of 0.
	 * @return The sum of squares for all data values that have
	 * a specific input value.
	 */
	double getLevelSumOfSquares(int indexOfLevel);	

	/**
	 *	Computes the variance of a single level.
	 *  <p>
	 *  The input parameter is NOT one of the values in the input vector.
	 *  The input parameter is an INDEX.  That is, after duplicates are
	 *  removed, the input data values are sorted.  The smallest distinct
	 *  input data value has an index of 0.  
	 * <p>
	 * EXAMPLE: <br>
	 * Suppose our (input,output) data values look like this: <br>
	 * {(5,1), (15,2), (10,3), (5,3), (10,5)} <br>
	 * We have 3 distinct input values: {5,15,10} <br>
	 * The input values are sorted:  {5,10,15} <br>
	 * The smallest value (5) is mapped to the index value 0. <br>
	 * The 2nd highest value (10) is mapped to the index value 1. <br>
	 * The 3rd highest value (15) is mapped to the index value 2. <br>
	 * Compute the variance for the y values {1,3} at index 0. <br>
	 * Compute the variance for the y values {3,5} at index 1. <br>
	 * Compute the variance for the y values {2} at index 2. <br>
	 *	@param int indexOfLevels: indentifies one of the levels.  
	 *  The lowest value in the input data has an index value of 0.
	 * @return The variance for all data values that have
	 * a specific input value.
	 */
	double getLevelVariance(int indexOfLevel);	

	/**
	 *	Computes the average output value for each and every level.
	 *  <p>
	 * EXAMPLE: <br>
	 * Suppose our (input,output) data values look like this: <br>
	 * {(5,1), (15,2), (10,3), (5,3), (10,5)} <br>
	 * We have 3 distinct input values: {5,15,10} <br>
	 * The input values are sorted:  {5,10,15} <br>
	 * The average y value for the smallest input value (5) is (1+3)/2 = 2. <br>
	 * The average y value for next highest input value (10) is (3+5)/2 = 4.<br>
	 * The average y value for the last input value (15) is 2. <br>
	 * The average output values are returned as the vector {2, 4, 15} <br>
	 * @return The average output values of groups of (input,output) data
	 * values.  
	 */
	std::vector<double> getAllLevelAverages();

  private:

	/**
	 * Find all (input,output) data pairs for data points that have
	 * a specific input data value.  The output data value from all found
	 * pairs are collected into a vector.  The vector is appended
	 * to the lvlResponses vector.
	 * <p>
	 * Example: <br>
	 * If the (input,output) data pairs are: <br>
	 * {(5,1), (15,2), (10,3), (15,4), (5,5)} <br>
	 * and if the inputDataValue is 15 <br>
	 * then {2,4} is placed inside a newly created vector <br>
	 * That vector is appended to the lvlResponses vector.
	 * @param indexOfInputDataValue The the input data
	 * value that we are looking for.
	 */
	void extractAllOutputDataValuesThatHaveThisInputDataValue
	     (int inputDataValue);

    /*
     * The values that we acquired, or computed, for each
     * and every run for one output data variable (or response).  
     */
     Response response_;


    /**
     * The values that we acquired for each
     * and every run from one input data variable (or factor).  The first
     * element contains the value that we acquired during the first run.
     * The value must be an positive integer and less than nLevels.
     */
     std::vector<int>	factors_;

     /*
      * The number of (input,output) data pairs.
      */
      int nObs_;

      /* 
       * number of distinct values for the factors
       * <p>
       * Example:  if the input data values are {0,0,1,1,1,2}
       * then the number of distinct values is 3.
       */
      int nLevels_;

    /**
     * Groups (input,output) pairs of data points into groups of
     * distinct input values.  The output values are extracted from
     * each group and placed inside the leveResponses_ vector.
     * <p>
     * EXAMPLE:  If the (input,output) pairs of data contain <br>
     * {(10,1), (15,1), (5,2), (15,3), (10,4)} <br>
     * then levelResponses_ would contain: <br>
     * { all output values that have an input value of 5}, <br>
     * { all output values that have an input value of 10}, <br>
     * { all output values that have an input value of 15} <br>
     * levelIndices_ = {{2}, {1,4}, {1,3}}. <br>
     */
     std::vector<Response>	levelResponses_;
};

}//namespace

#endif
