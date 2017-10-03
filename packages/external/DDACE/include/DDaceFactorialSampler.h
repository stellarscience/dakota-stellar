#ifndef DDACEFACTSAMPLER_H
#define DDACEFACTSAMPLER_H

#include "SmartPtr.h"
#include "DDaceSamplePoint.h"
#include "Distribution.h"
#include <algorithm>
#include "UniformDistribution.h"
#include <fstream>
#include <iostream>
#include <cmath>

/**
 * Sampler  class for generating  a Full  Factorial sampling.   This is
 * also called grid sampling in some contexts.
 * <br>
 * Goal:  Generate every possible combination of values for the
 * input factors.  
 * <br>
 * EXAMPLE: <br>
 * Add 1 or 10 drops of red food coloring to pitcher of water <br> 
 * AND <br>
 * Add 1 or 10 drops of green food coloring to the same pitcher <br>
 * AND <br>
 * Add 1 or 10 drops of blue food coloring to the same pitcher <br>
 * THEN <br>
 * 1) The input factors are:  red, green, blue <br>
 * 2) The symbols or levels are:  1Drop, 10Drops <br>
 * 3) The levels or symbols are : 1Drop, 10Drops <br>
 * 4) The samples are: <br>
 *    pitcher with  1 drop  of red,  1 drop  of green,  1 drop  of blue <br>
 *    pitcher with  1 drop  of red,  1 drop  of green, 10 drops of blue <br>
 *    pitcher with  1 drop  of red, 10 drops of green,  1 drop  of blue <br>
 *    pitcher with  1 drop  of red, 10 drops of green, 10 drops of blue <br>
 *    pitcher with 10 drops of red,  1 drop  of green,  1 drop  of blue <br>
 *    pitcher with 10 drops of red,  1 drop  of green, 10 drops of blue <br>
 *    pitcher with 10 drops of red, 10 drops of green,  1 drop  of blue <br>
 *    pitcher with 10 drops of red, 10 drops of green, 10 drop  of blue <br>
 * 5) numberOfSamples = 2SymbolsForRed * 2SymbolsForGreen * 2SymbolsForBlue<br>
 *    numberOfSamples = numberOfSymbols ^ numberOfInputFactors = 2^3 = 8 <br>
 * <br>
 * NOTES: <br>
 * 1) For a 2-symbol, or 2-level, experiment, 
 * the 2 symbols, or 2 levels, are usually represented as 
 * {-1, +1} or as {low,high}.  For example, a typical Full Factorial
 * chart would look something like this:<br>
 * &nbsp;&nbsp;&nbsp;&nbsp;-1 -1 -1 <br>
 * &nbsp;&nbsp;&nbsp;&nbsp;-1 -1 +1 <br>
 * &nbsp;&nbsp;&nbsp;&nbsp;-1 +1 -1 <br>
 * &nbsp;&nbsp;&nbsp;&nbsp;-1 +1 +1 <br>
 * &nbsp;&nbsp;&nbsp;&nbsp;+1 -1 -1 <br>
 * &nbsp;&nbsp;&nbsp;&nbsp;+1 -1 +1 <br>
 * &nbsp;&nbsp;&nbsp;&nbsp;+1 +1 -1 <br>
 * &nbsp;&nbsp;&nbsp;&nbsp;+1 +1 +1 <br>
 * 2) Most Full Factorial experiments use 2 symbols, or 2 levels.  
 * Those two symbols, or two levels,
 * represent the lowest possible value for a factor and the highest
 * possible value for a factor. <br>
 * 3) To get more information, try googling on "full factorial"
 */

#include "DDaceSampler.h"

class DDaceFactorialSampler : public DDaceSamplerBase
{
public:


  /**
   * Construct a Factorial Sampler.  
   * @param nSamples Number of samples or experiments.  For each
   * experiment, each input factor is assigned a value.  For 
   * example, one sample might be the experiment where I dropped
   * 1 drop of red, 10 drops of green, and 1 drop of blue food
   * coloring into a pitcher.
   * @param nSymbols Number of different values one input factor can have.
   * For example, I could put either 1 drop or 10 drops of red food coloring
   * into a pitcher of water.
   * @param noise If true, then a random amount of noise is mixed in with
   * the input values.  For example, instead of 1.0 drops of red food
   * coloring, you might have 0.95 drops or 1.05 drops.
   * @param dist An array of distributions (e.g. uniform, normal, etc.).  
   * Each input factor must be assigned one distribution.
   */
  DDaceFactorialSampler(int nSamples, int nSymbols, bool noise,
			const std::vector<Distribution>& dist);
  
   /**
   * DDaceFactorialSampler, initailizes the data fields and
   *   checks that all parameters are valid.
   * @param nSamples    the number of samples
   * @param nSymbols    the number of possible values an input can take
   *                    example: if x is a binary input it can take
   *                             values 0 or 1, so the number of symbols
   *                             for x is 2.
   */ 

  DDaceFactorialSampler(int nSamples, int nSymbols);
  
  virtual ~DDaceFactorialSampler(){;}


  /**
   * Generate an array of samples (or experiments or runs).  
   * Each sample contains a value 
   * for each input factor.  For example, if the input factors are
   * the amount of red food coloring, amount of green food coloring,
   * and the amount of blue food coloring, then the array might
   * look like this: <br>
   * {1DropOfRed, 1DropOfGreen, 1DropOfBlue}, <br>
   * {1DropOfRed, 1DropOfGreen, 10DropsOfBlue}, <br>
   * {1DropOfRed, 10DropsOfGreen, 1DropOfBlue}, <br>
   * {1DropOfRed, 10DropsOfGreen, 1DropOfBlue}, <br>
   * {10DropsOfRed, 1DropOfGreen, 1DropOfBlue}, <br>
   * {10DropsOfRed, 1DropOfGreen, 10DropsOfBlue}, <br>
   * {10DropsOfRed, 10DropsOfGreen, 1DropOfBlue}, <br>
   * {10DropsOfRed, 10DropsOfGreen, 1DropOfBlue}, <br>
   * <p>
   * SIDE EFFECT:  The input samplePoints array is re-populated.
   * @param samplePoints An array that will be populated with 
   * experiments or runs.
   */    
  virtual std::vector<DDaceSamplePoint>& getSamples(std::vector<DDaceSamplePoint>& samplePoints) const;
  virtual std::vector<std::vector<int> > getP() const ;

  virtual DDaceSamplerBase* clone() const;
  virtual void print(std::ostream& os) const;
  
  virtual const std::string& typeName() const {return typeName_;}
  virtual int getParameter(const std::string& parameterName) const;

private:
  
  static std::string typeName_;
  
  /**
   *  # of different values that one input factor can have.
   * <br>
   * Also known as the number of levels.
   * <br>
   * EXAMPLE:  How much red food coloring did you put inside the pitcher? <br>
   * If every sample contains either 1 drop or 10 drops,  
   * then you have 2 symbols or 2 levels.
   */
  int nSymbols_;

  mutable std::vector<std::vector<int> > symbolMap_;

};

#endif
