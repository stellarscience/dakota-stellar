#include "DDaceFactorialSampler.h"

using namespace std;

string DDaceFactorialSampler::typeName_ = "DDaceFactorialSampler";

/**
 * DDaceFactorialSampler, initailizes the data fields and
 *   checks that all parameters are valid.
 * @param nSamples    the number of samples
 * @param nSymbols    the number of possible values an input can take
 *                    example: if x is a binary input it can take
 *                             values 0 or 1, so the number of symbols
 *                             for x is 2.
 * @param noise       switches the use the use of noise on or off
 * @param dist        array of distributions, one for each input
 */
DDaceFactorialSampler::DDaceFactorialSampler(int nSamples, 
				     int nSymbols, 
				     bool noise,
				     const std::vector<Distribution>& dist)
  : DDaceSamplerBase(nSamples, dist.size(), noise, dist),
    nSymbols_(nSymbols), symbolMap_(0)
{
  /* compute # of samples or experiments or runs */	
  int testSamples = (int) pow((double)nSymbols_, nInputs_);
  
  /* Compare the computed value with the value the caller gave us */
  /* If the two values do not match, throw an exception       */
  if (testSamples != nSamples_) 
    throw runtime_error("DDaceFactorialSampler ctor : nSymbols ^ nInputs is not equal\n to nSamples.");

  /**
   * Check that nInputs equals dist.size(), this must be true or else
   * there is the possibility of a seg fault in the method getSamples().
   */
  if (nInputs_ != (int) dist_.size())
	throw runtime_error("DDaceFactorialSampler: nInputs not equal to dist.size()");
}


/**
 * DDaceFactorialSampler, initailizes the data fields and
 *   checks that all parameters are valid.
 * @param nSamples    the number of samples
 * @param nSymbols    the number of possible values an input can take
 *                    example: if x is a binary input it can take
 *                             values 0 or 1, so the number of symbols
 *                             for x is 2.
 */
DDaceFactorialSampler::DDaceFactorialSampler(int nSamples, int nSymbols) 
  : DDaceSamplerBase
       (nSamples,
       (int)floor((log10((double)nSamples)/log10((double)nSymbols))+0.5), 
       false ),
    nSymbols_(nSymbols), symbolMap_(0)
{

  /* compute # of samples or experiments or runs */		
  int testSamples = (int) pow( (double) nSymbols_, nInputs_);
  
  /* Compare the computed value with the value the caller gave us */
  /* If the two values do not match, throw an exception       */  
  if (testSamples != nSamples_) 
    throw runtime_error("DDaceFactorialSampler ctor : nSymbols ^ nInputs is not "
	       "equal\n to nSamples.");
	         	       
}

/**
 * If a constructor which doesn't set the distribution is used to create the
 * DDaceFactorialSampler object the getSamples() method will throw an
 * runtime exception unless the setDist() method is used first to
 * pass in an vector of distributions of length nInputs_.
 */
vector<DDaceSamplePoint>& DDaceFactorialSampler::getSamples(std::vector<DDaceSamplePoint>& samplePoints)
  const
{
  samplePoints.resize(nSamples_);
  symbolMap_.resize(nSamples_);
  double denom = (double) nSymbols_;
  std::vector<double> x(nInputs_);

  // symbolInt is a place holder for the distribution domain 
  // value for the s-th symbol of the i-th input.
  int symbolInt;
  double p;  // probability	      

  /**
   * Check that nInputs equals dist.size(), this must be true or else
   * there is the possibility of a seg fault in the method getSamples().
   */
  if (nInputs_ != (int) dist_.size())
    {
      throw runtime_error("DDaceFactorialSampler::getSamples: nInputs not equal to dist.length()");
    }

  int counter = 0;
  // n is the sample number currently being calculated
  for (int n = 0; n < nSamples_; n += nSymbols_)
    {
      // s is the input symbol number being used
      for (int s = 0; s < nSymbols_; s++)
	{
	  // i is the input being used 
	  for (int i = 0; i < nInputs_; i++)
	    {
	      symbolMap_[counter].resize(nInputs_); 
	      // generate all combinations of input variable values.
	      symbolInt = ( (int) ((s + n)/
				   ((int) pow((double)nSymbols_,i)))) % nSymbols_;
	      symbolMap_[counter][i] = symbolInt; 
	      // get the probability of this input.
	      if (noise_)
		{
		  p = (symbolInt + DistributionBase::uniformUnitDeviate())/denom;
		}
	      else
		{
		  p = (symbolInt + 0.5)/denom;
		}
	      // map the probability to its domain value for the distribution of
	      // the i-th input (compute x as CDF^-1(p))
	      /* PROBLEM:  What if the caller never created a distribution array? */
	      if (this->dist_.size()!=0) {
  	          x[i] = dist_[i].getDeviate(p);
	      } else {
	      	  Distribution distribution(UniformDistribution(0,100));
	      	  x[i] = distribution.getDeviate(p);
	      }
	    }//for
	  samplePoints[n + s] = DDaceSamplePoint((n+s), x);
          counter++;
	}
    }
	return samplePoints;
}

DDaceSamplerBase* DDaceFactorialSampler::clone() const
{
  DDaceSamplerBase * rtn = new DDaceFactorialSampler(*this);
  if (rtn==0) throw std::bad_alloc();
  return rtn;
}

void DDaceFactorialSampler::print(ostream& os) const
{
  os << "<Factorial "
		 << "samples=\"" << nSamples_ << "\" ";
	os << "symbols=\"" << nSymbols_ << "\" ";
	os << "perturb=\"";
	if (noise_)
		{
			os << "true\" ";
		}
	else
		{
			os << "false\" ";
		}
	os << "seed=\"" << DistributionBase::seed() << "\"/>";
}

int DDaceFactorialSampler::getParameter(const string& parameterName) const
{
  string tmp(parameterName);
  transform(tmp.begin(),tmp.end(),tmp.begin(),(int(*) (int)) toupper);

  if (tmp == "SYMBOLS")
    {
      return nSymbols_;
    }
  throw runtime_error("DDaceFactorialSampler::getParameter() : invalid parameter "
	     " name.");
  return 0;
}

std::vector<std::vector<int> > DDaceFactorialSampler::getP() const 
{
        return symbolMap_;
}

