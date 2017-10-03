#ifndef DDACELHSAMPLER_H
#define DDACELHSAMPLER_H

/**

Latin hypercube sampling

*/

#include "DDaceSampler.h"
#include <algorithm>
#include <vector>
#include "UniformDistribution.h"

class DDaceLHSampler : public DDaceSamplerBase
{
 public:
	DDaceLHSampler(int nSamples, int nReplications, 
		bool noise, const std::vector<Distribution>& dist);

	DDaceLHSampler(int nSamples, int nInputs,
		int nReplications, bool noise);

	virtual ~DDaceLHSampler(){;}
	
	virtual std::vector<DDaceSamplePoint>& getSamples(std::vector<DDaceSamplePoint>& samplePoints) const ;
        virtual std::vector<std::vector<int> > getP() const 
            {return permutationMatrix_;}


	virtual DDaceSamplerBase* clone() const ;
	virtual void print(std::ostream& os) const ;

	virtual const std::string& typeName() const {return typeName_;}
	virtual int getParameter(const std::string& parameterName) const ;
 private:
	void initPattern();
	
	static std::string typeName_;
	std::vector<std::vector<int> > permutationMatrix_;
	int nSymbols_;
	int nReplications_;
};

#endif
