#include "Distribution.h"

// use defines to switch between dlaran and the system rand().
//#define DDACE_USE_DLARAN
#define DDACE_USE_SYSTEMRAND

int DistributionBase::seed_ = 0;
int DistributionBase::seed48_[4] = {0, 0, 0, 0};

/* set up a pseudo random number generator that we can use during unit testing */
PseudoRandomTestsOnly DistributionBase::pseudoRandomTestsOnly;
bool DistributionBase::usePseudoRandomTestsOnly = false;

bool DistributionBase::seedSet_ = false;

void DistributionBase::setSeed(int seed)
{
	
  /* seed the pseudo random number generator */	
  pseudoRandomTestsOnly.setSeed(seed);	
	
  // MSE, 8/13/02: allow seed to be reset.
  //if (seedSet_)
  //  ExceptionBase::raise("DistributionBase::setSeed called twice");

  seed_ = seed;
  initRandom();
  seedSet_ = true;
#ifdef DDACE_USE_SYSTEMRAND
  srand(seed_);
#endif

}

int DistributionBase::seed() 
{
	
  /* seed the pseudo random number generator */
  pseudoRandomTestsOnly.setSeed(0);	
	
  if(!seedSet_)
    {
      setSeed(timeSeed());
    }
  return seed_;
}

int* DistributionBase::seed48() 
{
	
  /* seed the pseudo random number generator */
  pseudoRandomTestsOnly.setSeed(0);	
	
  if(seedSet_)
    return seed48_;
  else 
    {
      throw std::runtime_error("DistributionBase::seed48() : seed has not yet been set.");
      return 0; // -Wall
    }
}
  
int DistributionBase::timeSeed()
{
  struct timeval tp;
  // Modified from http://mywebpage.netscape.com/yongweiwu/timeval.h.txt
#if defined(_WIN32) || defined(_WIN64)
  union {
    __int64 ns100; /*time since 1 Jan 1601 in 100ns units */
    FILETIME ft;
  } now;
  
  GetSystemTimeAsFileTime (&now.ft);
  tp.tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL);
  tp.tv_sec = (long) ((now.ns100 - 116444736000000000LL) / 10000000LL);
#else
  gettimeofday(&tp,NULL);
#endif
  
  // use number of microseconds as seed
  int rtn = tp.tv_usec;
  
  //cerr << "time seed is " << rtn << endl;
  
  return rtn;
}

void DistributionBase::initRandom() 
{
  int r = 4096;
  int r2 = r*r;
  seed48_[0] = 0;
  seed48_[1] = seed_/r2;
  seed48_[2] = (seed_ - seed48_[1]*r2)/r;
  seed48_[3] = (seed_ - seed48_[1]*r2 - seed48_[2]*r);
  if ((seed48_[3] % 2) == 0)
    {
      seed48_[3]++;
    }
}


void DistributionBase::usePseudoRandom(bool x){
	usePseudoRandomTestsOnly = x;
}

double DistributionBase::uniformUnitDeviate()
{
  /* If we are doing a unit test, use a pseudo number generator */
  if (usePseudoRandomTestsOnly) 
    {	
      return(pseudoRandomTestsOnly.getPseudoRandom());
    }	

  if (!seedSet_)
    setSeed(timeSeed());
    
#ifdef DDACE_USE_DLARAN
  //double u = dlaran_((int*) seed48_);
  double u = DLARAN_F77((int *) seed48_);
  return u;
#elif defined DDACE_USE_SYSTEMRAND
  return ((double) rand())/((double) RAND_MAX);
#endif
}

#undef DDACE_USE_DLARAN
#undef DDACE_USE_SYSTEMRAND


Distribution::Distribution(const DistributionBase& base)
  : ptr_(base.clone())
{;}
