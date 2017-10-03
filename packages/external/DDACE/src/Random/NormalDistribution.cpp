#include "NormalDistribution.h"

using namespace std;

extern "C" 
{
  void cdfnor(int *which,double *p,double *q,double *x,double *mean,
	      double *sd,int *status,double *bound);
}


std::string NormalDistribution::typeName_ = "NormalDistribution";

NormalDistribution::NormalDistribution(const Mean& m, 
				       const StdDeviation& sigma, 
				       double numDeviations)
  : DistributionBase(), mean_(m.value()), stdDev_(sigma.value())
{
  if(stdDev_ < 0.0)
    throw std::runtime_error("NormalDistribution : in ctor, sigma must be positive.");
  if(numDeviations < 0.0)
    throw runtime_error("NormalDistribution : in ctor, numDevations must be positive.");
  lower_ = mean_ - numDeviations*stdDev_;
  upper_ = mean_ + numDeviations*stdDev_;
  pLower_ = getUntruncatedCDF(lower_);
  pUpper_ = getUntruncatedCDF(upper_);
}  


NormalDistribution::NormalDistribution(const Mean& m, 
				       const StdDeviation& sigma)
  : DistributionBase(), mean_(m.value()), stdDev_(sigma.value())
{
  if(stdDev_ < 0.0)
    throw runtime_error("NormalDistribution : in ctor, sigma must be positive.");
  lower_ = mean_ - stdDev_*2;
  upper_ = mean_ + stdDev_*2;
  pLower_ = getUntruncatedCDF(lower_);
  pUpper_ = getUntruncatedCDF(upper_);
}  


NormalDistribution::NormalDistribution(double lower, double upper)
  : DistributionBase()
{
  if ( lower > upper)
    throw runtime_error("NormalDistribution : in ctor, lower bound greater \nthan upper bound.");
  
  mean_ = (lower + upper)/2;
  lower_ = lower;
  upper_ = upper;
  
  int dummy = 4; // code number in dcdflib.c for "find standard deviation."
  double p = 0.025;
  double q = 1 - p;
  double tmpStdDev = 0.0;
  int status = 0; // will be reset for success or failure by cdfnor()
  double bound = 0.0; // will be reset for success or failure by cdfnor()
  // for more information on cdfnor() args, see the dcdflib.c 
  // documentation.  Don't really understand this parameter....
  
  cdfnor(&dummy, &p, &q, &lower, &mean_, &tmpStdDev, &status, &bound);
  
  if (!status==0)
    throw runtime_error("NormalDistribution: ctor error in calculating the standard"
	       " deviation.");
  
  stdDev_ = tmpStdDev;
  
  pLower_ = getUntruncatedCDF(lower_);
  pUpper_ = getUntruncatedCDF(upper_);
} 



NormalDistribution::NormalDistribution(double lower, 
			     double upper, 
			     double numDeviations)
  : DistributionBase()
{
  if ( lower > upper)
    throw runtime_error("NormalDistribution : in ctor, lower bound greater \nthan"
	       " upper bound.");
  if (numDeviations < 0)
    throw runtime_error("NormalDistribution : in ctor, numDevations must be positive.");
  
  lower_ = lower;
  upper_ = upper;  
  mean_ = (lower + upper)/2;
  stdDev_ = (upper - lower)/(2*numDeviations);

  pLower_ = getUntruncatedCDF(lower_);
  pUpper_ = getUntruncatedCDF(upper_);
}



// get a random value from the Normal distribution object.
double NormalDistribution::getDeviate() const
{
  return getDeviate(DistributionBase::uniformUnitDeviate());
}

// get a specific value from the Normal distribution object, given
// the cumulative probability of that normal value.
double NormalDistribution::getDeviate(double prob) const
{
  if((prob < 0.0) || (prob > 1.0))
    {
      cerr << "normal distribution " << prob << endl;
      throw runtime_error("NormalDistribution::getDeviate() : probability out of bounds.");
    }
  int dummy = 2; //code number in the dcdflib.c code for "find X."
  double p = pLower_ + prob*(pUpper_ - pLower_);
  double q = 1.0 - p;
  double x = 0.0; // the X we're looking for -- have to give it a 
  // default (what a stupid idea....)
  int status = 0; // will be reset for success or failure by cdfnor()
  double bound = 0.0; // will be reset for success or failure by cdfnor()
  // for more information on cdfnor() args, see the dcdflib.c 
  // documentation.  Don't really understand this parameter....
  
  // I guess the followin is a hack: compiler complains that mean_
  // and stdDev_ are of type const double, while the function cdfnor()
  // expects type double args.
  double tmpMean = mean_;
  double tmpSD = stdDev_;
  cdfnor(&dummy, &p, &q, &x, &tmpMean, 
	 &tmpSD, &status, &bound);
  
  if (!status==0)
    throw runtime_error("NormalDistribution::getDeviate() : error in inverse cdf"
	       " calculation.");
  
  return x;
}

void NormalDistribution::print(ostream& os) const
{
  os << "NORMAL MEAN " << mean() << " DEV " << stdDev() << " CUTOFF " 
		 << (upperBound()-lowerBound())/2.0/stdDev() << endl;
}

void NormalDistribution::printAttributes(ostream& os) const
{
os << "distribution=\"normal\" mean=\"" << mean() 
	 << "\" sigma=\"" << stdDev() << "\" cutoff=\"" << 
	(upperBound()-lowerBound())/2.0/stdDev() << "\"";
}

DistributionBase* NormalDistribution::clone() const
{
  DistributionBase* rtn = new NormalDistribution(*this);
  if (rtn == 0) throw std::bad_alloc();
  return rtn;
}

double NormalDistribution::getCDF(double x) const 
{
  double p = getUntruncatedCDF(x);
  return (p-pLower_)/(pUpper_ - pLower_);
}

double NormalDistribution::getUntruncatedCDF(double x) const
{
  int iwhich = 1; //code number in the dcdflib.c code for "find X."
  double p = 0.0; // these values need to be initialized. Stupid.
  double q = 1.0;

  int status = 0; // will be reset for success or failure by cdfnor()
  double bound = 0.0; // will be reset for success or failure by cdfnor()
  // for more information on cdfnor() args, see the dcdflib.c 
  // documentation.  Don't really understand this parameter....
  
  // I guess the followin is a hack: compiler complains that mean_
  // and stdDev_ are of type const double, while the function cdfnor()
  // expects type double args.
  double tmpMean = mean_;
  double tmpSD = stdDev_;
  cdfnor(&iwhich, &p, &q, &x, &tmpMean, 
	 &tmpSD, &status, &bound);
  
  if (!status==0)
    throw runtime_error("NormalDistribution::getUntruncatedCDF() : error in inverse cdf"
	       " calculation.");
  
  return p;
}

