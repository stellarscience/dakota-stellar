#ifndef DDACEARRAYSAMPLER_H
#define DDACEARRAYSAMPLER_H

/**

 Sampler class for sampling at a set of points specified by the user.

*/

#include <iostream>
#include <fstream>
#include "DDaceSampler.h"
#include "DDaceSamplePoint.h"

class DDaceArraySampler : public DDaceSamplerBase
{
 public:

        /**
         * Construct a sampler.
         * @param array The input data.   
         * array[i][] holds one value for each
         * input variable (x1, x2, x3, ...) for
         * the ith run or the ith experiment.  
         * array[i][j] holds the value of the jth input variable
         * during the ith run.  
         */
	DDaceArraySampler(std::vector< std::vector < double > > array);


        /**
         * Set the input data for the sampler.
         * @param array The input data.   
         * array[i][] holds one value for each
         * input variable (x1, x2, x3, ...) for
         * the ith run or the ith experiment.  
         * array[i][j] holds the value of the jth input variable
         * during the ith run.  
         */
        virtual void setInputData(std::vector < std::vector < double > >& array);


	virtual ~DDaceArraySampler(){;}
	
	virtual std::vector<DDaceSamplePoint>& getSamples(std::vector<DDaceSamplePoint>& samplePoints) const 
		{samplePoints = pts_;
		return samplePoints;}

	virtual const std::vector<Distribution>& dist() const ;
	virtual std::vector<double> lowerBounds() const ;
	virtual std::vector<double> upperBounds() const ;

	virtual DDaceSamplerBase* clone() const ;
	virtual void print(std::ostream& os) const ;

	virtual const std::string& typeName() const {return typeName_;}
	virtual int getParameter(const std::string& parameterName) const ;
	virtual std::vector<std::vector<int> > getP() const ;
 private:
	
	static std::string typeName_;
	std::string ptFilename_;
	std::vector<DDaceSamplePoint> pts_;
	std::vector<double> lowerBounds_;
	std::vector<double> upperBounds_;
};

#endif 
