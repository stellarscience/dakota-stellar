#include "UniformDistribution.h"

std::string UniformDistribution::typeName_ = "UniformDistribution";

UniformDistribution::UniformDistribution()
  : DistributionBase(), lowerBound_(0.0), upperBound_(1.0)
{;}



UniformDistribution::UniformDistribution(double lowerBound, 
			       double upperBound)
  : DistributionBase()
{
  if (lowerBound > upperBound)
    throw std::runtime_error("UniformDistribution : in ctor, lower bound greater \nthan upper bound.");
  else
    {
      lowerBound_ = lowerBound;
      upperBound_ = upperBound;
    }
}



double UniformDistribution::getDeviate() const
{
  return (lowerBound_ + 
	  (upperBound_ - lowerBound_)*DistributionBase::uniformUnitDeviate());
}
double UniformDistribution::getDeviate( double prob) const
{
  return (lowerBound_ + 
	  (upperBound_ - lowerBound_)*prob);
}

double UniformDistribution::getCDF(double x) const 
{
  return (x-lowerBound())/(upperBound()-lowerBound());
}

double UniformDistribution::lowerBound() const
{ return lowerBound_;}

double UniformDistribution::upperBound() const
{ return upperBound_;}

double UniformDistribution::mean() const
{
  return (lowerBound_ + upperBound_)/2;
}

double UniformDistribution::stdDev() const
{
  return sqrt(pow((upperBound_ - lowerBound_), 2)/12);
}

void UniformDistribution::print(std::ostream& os) const
{
  os << "UNIFORM " << lowerBound() << " " << upperBound();
}

void UniformDistribution::printAttributes(std::ostream& os) const
{
  os << "distribution=\"uniform\" lower=\"" << lowerBound() 
		 << "\" upper=\"" << upperBound() << "\"";
}

DistributionBase* UniformDistribution::clone() const
{
  DistributionBase* rtn = new UniformDistribution(*this);
  if (rtn == 0) throw std::bad_alloc();
  return rtn;
}

std::vector<int> randomIVector(const int leng)
{
  std::vector<int> rtn(leng);
  int    k, iran1, iran2, itmp;

  // uses default time seed.
  //  UniformDistribution unif(0.0,1.0);
  for (k=0; k<leng; k++) { rtn[k] = k;}
  for (k=0; k<3*leng; k++) {
    //iran1 = (int) (leng * unif.getDeviate());
    //iran2 = (int) (leng * unif.getDeviate());
    iran1 = (int) (leng * DistributionBase::uniformUnitDeviate());
    iran2 = (int) (leng * DistributionBase::uniformUnitDeviate());
    iran1 = (iran1 == leng) ? 0 : iran1;
    iran2 = (iran2 == leng) ? 0 : iran2;
    itmp  = rtn[iran2];
    rtn[iran2] = rtn[iran1];
    rtn[iran1] = itmp;
  }
  return rtn;
}

