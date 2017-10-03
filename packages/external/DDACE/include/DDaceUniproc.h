#ifndef DDACEUNIPROC_H
#define DDACEUNIPROC_H

#include "DDaceBool.h"
#include "DDaceSampler.h"
#include "DDaceServer.h"
#include "XMLObject.h"

/**

 DDaceUniproc defines the interface that will be used by both the
 client and server in a uniprocessor context.

 */
class DDaceUniproc : public DDaceServer
{
 public:
	DDaceUniproc(const DDaceSampler& sampler, const Array<String>& varNames,
							 const Array<String>& outputNames, 
							 const String& archiveFilename);
	DDaceUniproc(const DDaceSampler& sampler);
	DDaceUniproc(const XMLObject& xmlObj);

	virtual ~DDaceUniproc(){;}

	// sampling functions
	virtual bool getNextSample(DDaceSamplePoint& pt) ;
	virtual void storeFunctionValue(const DDaceSamplePoint& pt, 
																	const Array<double>& values) ;
	virtual void recordRunStatus(const DDaceSamplePoint& pt, 
															 DDaceRunStatus status);
 protected:
};

#endif


