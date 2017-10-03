#ifndef DDACESERVER_H
#define DDACESERVER_H

#include "DDaceBool.h"
#include "DDaceMachineBase.h"
#include "DDaceSampler.h"
#include "DDaceRunStatus.h"
#include "XMLObject.h"

/**
   DDaceServer class is very simple:
   its only overloaded methods are a trivial virtual dtor and 
   the methods to get sample values and set results.
  
   The getNextSample() method sends the client a new sample point.  If
   all the points are used up, getNextSample() notifies the client by
   returning false.
   
   storeFunctionValue() stores a function evaluation sent by the client
 
  */

class DDaceServer : public DDaceMachineBase
{
 public:
	DDaceServer(const DDaceSampler& sampler, const Array<String>& varNames,
							const Array<String>& outputNames, const String& archiveFilename);
	DDaceServer(const DDaceSampler& sampler);
	DDaceServer(const XMLObject& xmlObj);

	virtual ~DDaceServer(){;}

	// sampling functions
	virtual bool getNextSample(DDaceSamplePoint& pt) ;
	virtual void getSampler(DDaceSampler& sampler) const;
	virtual void getRunStatus(Array<DDaceRunStatus>& status) const;
	virtual void getArchiveFilename(String& archiveFilename) const;
	virtual void getVarNames(Array<String>& names) const;
	virtual void getOutputNames(Array<String>& names) const;
	virtual int  getNumOutputs() const;
	virtual void storeFunctionValue(const DDaceSamplePoint& pt, 
					const Array<double>& values) ;
	virtual void recordRunStatus(const DDaceSamplePoint& pt, 
															 DDaceRunStatus status);
	virtual void writeToArchive();

	virtual void setArchivePath(const String& archivePath);
	virtual void setArchiveName(const String& archiveFilename);
	virtual void setVariableNames(const Array<String>& name);
	virtual void setOutputNames(const Array<String>& name);

	virtual void getResults(Array<DDaceSamplePoint>& pts, 
													Array<Array<double> >& funcValues) const ;
 protected:

	int stackPtr_;
	Array<Array<double> > results_;
	Array<DDaceRunStatus> status_;
	DDaceSampler sampler_;
	String archiveFilename_;
	String archivePath_;
	Array<String> varNames_;
	Array<String> outputNames_;
};

#endif
