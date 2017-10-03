#include "DDaceSampler.h"

using std::bad_alloc;
using std::runtime_error;
using std::ostream;

int DDaceSamplerBase::getParameter(const std::string& parameterName) const
{
  throw std::runtime_error("DDaceSamplerBase::getParameter not defined for base class");
  return 0;
}

std::vector<std::vector<int> > DDaceSamplerBase::getP() const
{
  throw std::runtime_error("DDaceSamplerBase::getP not defined for base class");
        std::vector<std::vector<int> > tmp;
//        for(int row = 0; row < 2; row++){
//           std::vector<int> oneTmp;
//           for(int col = 0; col < 2; col++)
//              oneTmp.push_back(0);
//           tmp.push_back(oneTmp);
//        }
	return tmp;
           
}

std::vector<double> DDaceSamplerBase::upperBounds() const
{
	std::vector<double> ub(dist_.size());
	for (int i=0; i< (int) dist_.size(); i++)
		{
			ub[i] = dist_[i].upperBound();
		}
	return ub;
}

std::vector<double> DDaceSamplerBase::lowerBounds() const
{
	std::vector<double> lb(dist_.size());
	for (int i=0; i< (int) dist_.size(); i++)
		{
			lb[i] = dist_[i].lowerBound();
		}
	return lb;
}

DDaceSampler::DDaceSampler(const DDaceSamplerBase& base)
	: ptr_(base.clone())
{;}


std::vector<DDaceSamplePoint>& DDaceSampler::getSamples(std::vector<DDaceSamplePoint>& samplePoints) const
{
	if(ptr_)
	{
		return ptr_->getSamples(samplePoints);	
	}
	else
	{
	  throw std::bad_alloc();
	}
}

void DDaceSampler::print(std::ostream& os) const
{
	if(ptr_)
	{
	ptr_->print(os);
	}
	else
	{
	  throw std::bad_alloc();
	}
}

const std::string& DDaceSampler::typeName() const 
{
	return ptr_->typeName();
}

std::vector<std::vector<int> > DDaceSampler::getP() const 
{
	return ptr_->getP();
}


int DDaceSampler::nSamples() const 
{
	return ptr_->nSamples();
}

int DDaceSampler::nInputs() const 
{
	return ptr_->nInputs();
}

const std::vector<Distribution>& DDaceSampler::dist() const 
{
	return ptr_->dist();
}

std::vector<double> DDaceSampler::upperBounds() const 
{
	return ptr_->upperBounds();
}

std::vector<double> DDaceSampler::lowerBounds() const 
{
	return ptr_->lowerBounds();
}

bool DDaceSampler::noise() const
{
  return ptr_->noise();
}


int DDaceSampler::getParameter(const std::string& parameterName) const 
{
	return ptr_->getParameter(parameterName);
}



std::ostream& operator<<(std::ostream& os, const DDaceSampler& sampler)
{
	sampler.print(os);
	return os;
}
