#ifndef DDACEUNIPROC_H
#define DDACEUNIPROC_H

#include "DDaceSampler.h"
#include "DDaceServer.h"

class DDaceUniproc : public DDaceServer
{
 public:
	DDaceUniproc(const DDaceSampler& sampler,
		 const std::vector<double>& lower,
		const vector<double>& upper, 
		const vector<std::string>& varNames,
		const vector<string>& outputNames, 
		const string& archiveFilename);

	DDaceUniproc(const DDaceSampler& sampler, 
			const vector<double>& lower,
			const vector<double>& upper);

	virtual ~DDaceUniproc(){;}

	// sampling functions
	virtual bool getNextSample(DDaceSamplePoint& pt) ;
	virtual void storeFunctionValue(const DDaceSamplePoint& pt, 
					const vector<double>& values) ;
	virtual void recordRunStatus(const DDaceSamplePoint& pt,							 DDaceRunStatus status);
 protected:
};

#endif


