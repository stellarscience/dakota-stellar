#include "DDaceRandomSampler.h"

using namespace std;

std::string DDaceRandomSampler::typeName_ = "DDaceRandomSampler";

DDaceRandomSampler::DDaceRandomSampler(int nSamples, 
					 const std::vector<Distribution>& dist)
	: DDaceSamplerBase(nSamples, dist.size(), false, dist)
{
        /**
         * Check that nInputs equals dist_.size(), this must be true or else
         * there is the possibility of a seg fault in the method getSamples().
         */
        if (nInputs_ != (int) dist_.size())
          {
            throw runtime_error("DDaceRandomSampler: nInputs not equal to dist.length()");
          }

}

DDaceRandomSampler::DDaceRandomSampler(int nSamples, int nInputs)
	: DDaceSamplerBase(nSamples, nInputs, false)
{;}

/**
 * If a constructor which doesn't set the distribution is used to create the
 * DDaceRandomSampler object the getSamples() method will throw an
 * ArgumentMisMatchException unless the setDist() method is used first to
 * pass in an array of distributions of length nInputs_.
 */

vector<DDaceSamplePoint>& DDaceRandomSampler::getSamples(std::vector<DDaceSamplePoint>& samplePoints) const
{
	std::vector<double> x(nInputs_);
	samplePoints.resize(nSamples_);

        /**
         * Check that nInputs equals dist.size(), this must be true or else
         * there is the possibility of a seg fault in the method getSamples().
         */
        if (nInputs_ != (int) dist_.size())
          {
            throw runtime_error("DDaceRandomSampler::getSamples: nInputs not equal to dist.length()");
          }

	for (int s=0; s<nSamples_; s++)
		{
			for (int i=0; i<nInputs_; i++)
				{
					x[i] = dist_[i].getDeviate();
				}
			samplePoints[s] = DDaceSamplePoint(s, x);
		}
return samplePoints;
}
	
DDaceSamplerBase* DDaceRandomSampler::clone() const
{
	DDaceSamplerBase* rtn = new DDaceRandomSampler(*this);
	if (rtn==0) throw std::bad_alloc();
	return rtn;
}

void DDaceRandomSampler::print(ostream& os) const
{
	os 
                 << "<Random "
		 << "samples=\"" << nSamples_ << "\" "
		 << "seed=\"" << DistributionBase::seed() << "\"/>";
}
	
std::vector<std::vector<int> > DDaceRandomSampler::getP() const 
{
        throw runtime_error("DDaceSamplerBase::getP not defined for base class");
        std::vector<std::vector<int> > tmp;
        return tmp;
}
	
