#ifndef DDACEUSERINPUTSAMPLER_H
#define DDACEUSERINPUTSAMPLER_H

/**

 Sampler class for sampling at a set of points specified by the user.

*/

#include "DDaceSampler.h"
#include "DDaceSamplePoint.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

class DDaceUserInputSampler : public DDaceSamplerBase
{
 public:
	DDaceUserInputSampler(const std::string& ptFilename);
	virtual ~DDaceUserInputSampler(){;}
	
	virtual std::vector<DDaceSamplePoint>& getSamples(std::vector<DDaceSamplePoint>& samplePoints) const 
		{samplePoints = pts_;
		 return samplePoints; }
        virtual std::vector<std::vector<int> > getP() const ;

	virtual const std::vector<Distribution>& dist() const ;
	virtual std::vector<double> lowerBounds() const ;
	virtual std::vector<double> upperBounds() const ;

	virtual DDaceSamplerBase* clone() const ;
	virtual void print(std::ostream& os) const ;

	virtual const std::string& typeName() const {return typeName_;}
	virtual int getParameter(const std::string& parameterName) const ;

	static std::vector<std::vector<std::string> > tokenizeFile(std::istream& is, char comment);
                                                                                
        static std::vector<std::string> stringTokenizer(const std::string& str);                                                                                
        static int findNextWhitespace(const std::string& str, int offset);
                                                                                
        static int findNextNonWhitespace(const std::string& str, int offset);
                                                                                
 private:
	
	static std::string typeName_;
	std::string ptFilename_;
	std::vector<DDaceSamplePoint> pts_;
	std::vector<double> lowerBounds_;
	std::vector<double> upperBounds_;
};

#endif
