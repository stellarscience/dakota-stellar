#ifndef DDACEBOXBEHNKENSAMPLER_H
#define DDACEBOXBEHNKENSAMPLER_H

/**

Box-Behnken sampling - a deterministic design of
experiments technique for parameter space sampling

*/

/*

Tony Giunta, 19 March 2001
             01 May   2002, modified constructor arguments and
                            getSamples() function

*/

#include "DDaceSampler.h"


class DDaceBoxBehnkenSampler : public DDaceSamplerBase
{
 public:
	DDaceBoxBehnkenSampler(int nSamples, int nInputs,
			       const std::vector<Distribution>& dist);
	virtual ~DDaceBoxBehnkenSampler(){;}
	
	virtual std::vector<DDaceSamplePoint>& getSamples(std::vector<DDaceSamplePoint>& samplePoints) const ;
	virtual std::vector<std::vector<int> > getP() const ;
	
	virtual DDaceSamplerBase* clone() const ;
	virtual void print(std::ostream& os) const ;
	virtual const std::string& typeName() const {return typeName_;}

 private:
	static std::string typeName_;
};

#endif
