#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
  #include "ddace_config.h"
  // lapack's uniform random number generator.
  #define DLARAN_F77 F77_FUNC(dlaran,DLARAN)
#else
  #include "ddace_fortran.h"
  #define DLARAN_F77 DDACE_FORTRAN_GLOBAL(dlaran,DLARAN)
#endif

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <sys/time.h>
//Doesn't appear to be used:
//#include <unistd.h>
#endif

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>

#include "SmartPtr.h"
#include "PseudoRandomTestsOnly.h"

#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif
double DLARAN_F77(int *);

/**\ingroup Random 
 * Base class for distributions for random number generator.
 * @see Distribution
 * @author Leslea Lehoucq
 */
class DistributionBase
{
 public:
  DistributionBase(){;}
  
  virtual ~DistributionBase() {;}
  
  virtual DistributionBase* clone() const = 0;

  /** Return the seed value (set by the user, or generated from the clock) */
  static int seed(); 
  
  /** use seed to create a 48-bit seed to dlaran */
  static int* seed48();

  /** generate a random deviate from this distribution */
  virtual double getDeviate() const = 0;

  /** generate the x-value such that CDF(x) = prob */
  virtual double getDeviate(double prob) const = 0;
  
  /** evaluate the CDF at a point x */
  virtual double getCDF(double x) const = 0;

  /** get the lowest x than can be drawn from this distribution */
  virtual double lowerBound() const = 0;

  /** get the highest x than can be drawn from this distribution */
  virtual double upperBound() const = 0;

  /** get the mean of the distribution */
  virtual double mean() const = 0;

  /** get the standard deviation of the distribution */
  virtual double stdDev() const = 0;

  /** print a description of the distribution */
  virtual void print(std::ostream& os) const = 0;

  virtual void printAttributes(std::ostream& os) const = 0;
  
  virtual const std::string& typeName() const = 0 ;

  /** set the seed to be used in generating uniform deviates */
  static void setSeed(int seed);

  /** call dlaran() to generate a uniform deviate */
  static double uniformUnitDeviate();

  /** get a seed by sampling the milliseconds on the system clock */
  static int timeSeed(); 
  static bool seedSet()  {return seedSet_;}
			
			
  /* When doing unit testings, set this flag to true.    */
  /* It will disable the Random Number Generator.        */
  /* The Number Generator will output the sequence       */
  /* (0.001, 0.002, 0.003, etc.)                         */
  static void usePseudoRandom(bool x);
			
 protected:

  static void initRandom();

  static int seed_;
  static int seed48_[4];
  static bool seedSet_;
			
 private:
		
  /* When doing unit testings, set this flag to true.    */
  /* It will disable the Random Number Generator.        */
  /* The Number Generator will output the sequence       */
  /* (0.001, 0.002, 0.003, etc.)                         */
  static bool usePseudoRandomTestsOnly;
			
  static PseudoRandomTestsOnly 
    pseudoRandomTestsOnly;		
};


/**
 * Handle for distributions for random number generator. 
 * @see DistributionBase
 * @author Leslea Lehoucq
 */
class Distribution 
{
 public:
  Distribution() : ptr_(0) {;}
  Distribution(const DistributionBase& base);
  
  /** draw a random deviate from the distribution */
  double getDeviate() const throw();
	
  /** get the variable x such that CDF(x)=prob */
  double getDeviate(double prob) const throw();

  /** evaluate the cumulative distribution function at a point x */
  double getCDF(double x) const ;

  /** get the lowest value of x that can be drawn from this distribution */
  double lowerBound() const;

  /** get the highest value of x that can be drawn from this distribution */
  double upperBound() const;

  /** get the mean of the distribution */
  double mean() const;

  /** get the standard deviation of the distribution */
  double stdDev() const;

  /** print a description of the distribution */
  void print(std::ostream& os) const;
  void printAttributes(std::ostream& os) const ;

  /** set the seed to be used in the underlying 
   * uniform random number generator */
  static void setSeed(int seed) {DistributionBase::setSeed(seed);}

  /** Create a seed based on the current millisecond count */
  static int timeSeed() {return DistributionBase::timeSeed();}
  static bool seedSet() {return DistributionBase::seedSet();}

  /** return the class name of the distribution */
  const std::string& typeName() const ;
 private:
  SmartPtr<DistributionBase> ptr_;
};
 
inline double Distribution::getDeviate() const throw()
{ return ptr_->getDeviate();}

inline double Distribution::getDeviate(double x) const throw()
{ return ptr_->getDeviate(x);}

inline double Distribution::getCDF(double x) const
{ return ptr_->getCDF(x);}


inline double Distribution::lowerBound() const
{ return ptr_->lowerBound();}

inline double Distribution::upperBound() const
{ return ptr_->upperBound();}

inline double Distribution::mean() const
{ return ptr_->mean();}

inline double Distribution::stdDev() const
{ return ptr_->stdDev();}

inline void Distribution::print(std::ostream& os) const
{ ptr_->print(os);}

inline void Distribution::printAttributes(std::ostream& os) const
{ ptr_->printAttributes(os);}

inline const std::string& Distribution::typeName() const 
{
  return ptr_->typeName();
}

inline std::ostream& operator<<(std::ostream& os, const Distribution& dist)
{
  dist.print(os);
  return os;
}

#endif // DISTRIBUTION_H


