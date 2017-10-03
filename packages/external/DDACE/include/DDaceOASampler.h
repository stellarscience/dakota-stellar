#ifndef DDACEOASAMPLER_H
#define DDACEOASAMPLER_H

/**

Orthogonal Array sampling

*/

#include "DDaceSampler.h"
#include "UniformDistribution.h"


class DDaceOASampler : public DDaceSamplerBase
{
 public:
	DDaceOASampler(int nSamples,bool noise,const std::vector<Distribution>& dist);
	DDaceOASampler(int nSamples, int nInputs, bool noise); 

	virtual ~DDaceOASampler(){;}
	
	virtual std::vector<DDaceSamplePoint>& getSamples(std::vector<DDaceSamplePoint>& samplePoints) const ;
        virtual std::vector<std::vector<int> > getP() const {return symbolMap_;}
	
	virtual DDaceSamplerBase* clone() const ;
	virtual void print(std::ostream& os) const ;
	virtual const std::string& typeName() const {return typeName_;}
 private:
	void initPattern();
	
	std::vector<std::vector<int> > symbolMap_;
	int nSymbols_;

	static std::string typeName_;
};

#endif

