#ifndef DDACERANDOMSAMPLER_H
#define DDACERANDOMSAMPLER_H

/**

 Sampler class for generating pure (quasi)random samples

*/

#include "DDaceSampler.h"
#include "UniformDistribution.h"

class DDaceRandomSampler : public DDaceSamplerBase
{
 public:
  DDaceRandomSampler(int nSamples, 
		     const std::vector<Distribution>& dist);
  
  DDaceRandomSampler(int nSamples, int nInputs); 

  virtual ~DDaceRandomSampler(){;}
  
  virtual std::vector<DDaceSamplePoint>& getSamples(std::vector<DDaceSamplePoint>& samplePoints) const ;
  virtual std::vector<std::vector<int> > getP() const ;

  virtual DDaceSamplerBase* clone() const ;
  virtual void print(std::ostream& os) const ;
  virtual const std::string& typeName() const {return typeName_;}
 private:
  static std::string typeName_;
};

#endif
