#ifndef DDACECENTRALCOMPOSITESAMPLER_H
#define DDACECENTRALCOMPOSITESAMPLER_H

/**

Central-Composite-Design sampling - a deterministic design
of experiments technique for parameter space sampling

*/

/*

Tony Giunta, 30 April 2001
             01 May   2002, modified constructor arguments and
                            getSamples() function

*/

#include "DDaceSampler.h"
#include "UniformDistribution.h"
#include <cmath>

class DDaceCentralCompositeSampler : public DDaceSamplerBase
{
 public:
	DDaceCentralCompositeSampler(int nSamples, int nInputs, 
				     const std::vector<Distribution>& dist);
	virtual ~DDaceCentralCompositeSampler(){;}
	
	virtual std::vector<DDaceSamplePoint>& getSamples(std::vector<DDaceSamplePoint>& samplePoints) const ;
        virtual std::vector<std::vector<int> > getP() const ;

	virtual DDaceSamplerBase* clone() const ;
	virtual void print(std::ostream& os) const ;
	virtual const std::string& typeName() const {return typeName_;}
 private:
	static std::string typeName_;
	void fullFactorialPoints(int, 
				 int, 
				 const std::vector<double>&, 
				 const std::vector<double>&, 
				 const std::vector<double>&, 
				 std::vector<double>&, 
				 std::vector<DDaceSamplePoint>&,
				 int&) const ;
};

#endif
