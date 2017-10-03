#ifndef DDACESAMPLER_H
#define DDACESAMPLER_H

#include "SmartPtr.h"
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "DDaceSamplePoint.h"
#include "Distribution.h"


/**
 * DDaceSamplerBase: Abstract Data Type for Samplers. This class
 * contains the most common functions needed by most of the samplers.
 * <br>
 * EXAMPLE: <br>
 * Add 1 or 10 drops of red food coloring to pitcher of water <br> 
 * AND <br>
 * Add 1 or 10 drops of green food coloring to the same pitcher <br>
 * AND <br>
 * Add 1 or 10 drops of blue food coloring to the same pitcher <br>
 * THEN <br>
 * 1) The input factors are:  red, green, blue <br>
 * 2) The number of input factors is 3 <br>
 * 3) If we wanted, we could generate 8 samples <br>
 *    pitcher with  1 drop  of red,  1 drop  of green,  1 drop  of blue <br>
 *    pitcher with  1 drop  of red,  1 drop  of green, 10 drops of blue <br>
 *    pitcher with  1 drop  of red, 10 drops of green,  1 drop  of blue <br>
 *    pitcher with  1 drop  of red, 10 drops of green, 10 drops of blue <br>
 *    pitcher with 10 drops of red,  1 drop  of green,  1 drop  of blue <br>
 *    pitcher with 10 drops of red,  1 drop  of green, 10 drops of blue <br>
 *    pitcher with 10 drops of red, 10 drops of green,  1 drop  of blue <br>
 *    pitcher with 10 drops of red, 10 drops of green, 10 drop  of blue <br>
*/

class DDaceSamplerBase
{
 public:
 
  /**
   * Construct a DDaceSamplerBase.  
   * @param nSamples Number of samples or experiments.  For each
   * experiment, each input factor is assigned a value.  For 
   * example, one sample might be the experiment where I dropped
   * 1 drop of red, 10 drops of green, and 1 drop of blue food
   * coloring into a pitcher.
   * @param nInput Number of different input factors.  For example,
   * if, for each experiment, I change the amount of red food coloring,
   * change the amount of green food coloring, and change the amount
   * of blue food coloring, then I wouuld have 3 input factors.
   * @param noise If true, then a random amount of noise is mixed in with
   * the input values.  For example, instead of 1.0 drops of red food
   * coloring, you might have 0.95 drops or 1.05 drops.
   * @param dist An array of distributions (e.g. uniform, normal, etc.).  
   * Each input factor must be assigned one distribution.
   */ 
  DDaceSamplerBase(int nSamples, int nInputs, bool noise,
		   const std::vector<Distribution>& dist)
    : 
    nSamples_(nSamples), nInputs_(nInputs), noise_(noise),
    dist_(dist){;}


  /**
   * Construct a DDaceSamplerBase.  
   * @param nSamples Number of samples or experiments.  For each
   * experiment, each input factor is assigned a value.  For 
   * example, one sample might be the experiment where I dropped
   * 1 drop of red, 10 drops of green, and 1 drop of blue food
   * coloring into a pitcher.
   * @param nInput Number of different input factors.  For example,
   * if, for each experiment, I change the amount of red food coloring,
   * change the amount of green food coloring, and change the amount
   * of blue food coloring, then I wouuld have 3 input factors.
   * @param noise If true, then a random amount of noise is mixed in with
   * the input values.  For example, instead of 1.0 drops of red food
   * coloring, you might have 0.95 drops or 1.05 drops.
   */ 
  DDaceSamplerBase(int nSamples, int nInputs, bool noise)
    : 
    nSamples_(nSamples), nInputs_(nInputs), noise_(noise),
    dist_(0){;}
    
    
  
  virtual ~DDaceSamplerBase() {;}


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
   * SIDE EFFECT:  The input samplePoints array is re-populated.
   * @param samplePoints An array that will be populated with 
   * experiments or runs.
   */
  virtual std::vector<DDaceSamplePoint>& getSamples(std::vector<DDaceSamplePoint>& samplePoints) const = 0;
  virtual std::vector<std::vector<int> > getP() const = 0;

  virtual DDaceSamplerBase* clone() const = 0;
  
  virtual void print(std::ostream& os) const = 0;
  
  virtual const std::string& typeName() const = 0;
  
  /**
   * Return the number of samples or experiments.  For each
   * experiment, each input factor is assigned a value.  For 
   * example, one sample might be the experiment where I dropped
   * 1 drop of red, 10 drops of green, and 1 drop of blue food
   * coloring into a pitcher.
   */  
  virtual int nSamples() const {return nSamples_;}
  
  /**
   * Return the number of different input factors.  For example,
   * if, for each experiment, I change the amount of red food coloring,
   * the amount of green food coloring, and the amount
   * of blue food coloring, then I wouuld have 3 input factors.
   */  
  virtual int nInputs() const {return nInputs_;}
  
  /**
   * Return an array of distributions (e.g. uniform, normal, etc.).  
   * Each input factor must be assigned one distribution.  
   */
  virtual const std::vector<Distribution>& dist() const {return dist_;}
  
  
  virtual std::vector<double> lowerBounds() const ;
  
  virtual std::vector<double> upperBounds() const ;
  
  /**
   * Is a random amount of noise mixed in with
   * the input values.  For example, instead of 1.0 drops of red food
   * coloring, you might have 0.95 drops or 1.05 drops.
   */  
  virtual bool noise() const {return noise_;}
  
  
  virtual int getParameter(const std::string& parameterName) const ;
  
 protected:
 
  /**
   * Number of samples or experiments.  For each
   * experiment, each input factor is assigned a value.  For 
   * example, one sample might be the experiment where I dropped
   * 1 drop of red, 10 drops of green, and 1 drop of blue food
   * coloring into a pitcher.
   */
  int nSamples_;
  
  /**
   * Number of different input factors.  For example,
   * if, for each experiment, I change the amount of red food coloring,
   * change the amount of green food coloring, and change the amount
   * of blue food coloring, then I wouuld have 3 input factors.
   */
  int nInputs_;
  
  /**
   * If true, then a random amount of noise is mixed in with
   * the input values.  For example, instead of 1.0 drops of red food
   * coloring, you might have 0.95 drops or 1.05 drops.
   */
  bool noise_;
  
  /**
   * The starting seed for a random number generator.
   */
  int seed_;
  
  /**
   * An array of distributions (e.g. uniform, normal, etc.).  
   * Each input factor must be assigned one distribution.
   */
  std::vector<Distribution> dist_;
};


/**
 * DDaceSampler:
 * Base class for Samplers. This class contains the most common functions
 * needed by most of the samplers.
 * 
 */

class DDaceSampler
{
 public:
  DDaceSampler() : ptr_(0) {;}
  DDaceSampler(const DDaceSamplerBase& base);
  
  std::vector<DDaceSamplePoint>& getSamples(std::vector<DDaceSamplePoint>& samplePoints) const ;
  
  std::vector<std::vector<int> > getP() const;
  void print(std::ostream& os) const ;
  const std::string& typeName() const ;
  
  int nSamples() const ;
  int nInputs() const ;
  int getParameter(const std::string& parameterName) const ;
  
  const std::vector<Distribution>& dist() const ;
  std::vector<double> lowerBounds() const ;
  std::vector<double> upperBounds() const ;
  bool noise() const;
 private:
  SmartPtr<DDaceSamplerBase> ptr_;
};

std::ostream& operator<<(std::ostream& os, const DDaceSampler& sampler);

#endif
