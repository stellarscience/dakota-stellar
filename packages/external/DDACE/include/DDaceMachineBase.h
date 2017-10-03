#ifndef DDACEMACHINEBASE_H
#define DDACEMACHINEBASE_H

#include "Array.h"
#include "Array.inl"
#include "String.h"
#include "DDaceSamplePoint.h"
#include "DDaceRunStatus.h"
#include "DDaceSampler.h"

/// Some human-readable message tags and status flags

enum MsgID {RequestTag, IndexTag, ResultIndexTag, ResultValueTag};
enum ProcStatus {Available, Busy, AskedForData, 
								 AskedForRequest, Faulty, Done};


/**

DDaceMachineBase defines the interface that will be used by both client
and server in a multiprocessor context, and by the single machine in
uniprocessor context. 

From this base class, we will derive DDaceClient, DDaceServer, 
and DDaceUniproc. In multiprocessor mode, client machines will construct 
DDaceClient objects while the server constructs a DDaceServer. 

Server-specific methods, such as I/O methods, do nothing when 
called on the base class.

 */
class DDaceMachineBase
{
 public:
	DDaceMachineBase() : pts_(0) {;}
	virtual ~DDaceMachineBase(){;}

	// sampling and result storage methods, to be called in main sampling loop.
	virtual bool getNextSample(DDaceSamplePoint& pt) = 0;
	virtual void getSampler(DDaceSampler& sampler) const;
	virtual void getRunStatus(Array<DDaceRunStatus>& status) const;
	virtual void getArchiveFilename(String& archiveFilename) const;
	virtual void getVarNames(Array<String>& names) const;
	virtual void getOutputNames(Array<String>& names) const;
	virtual int  getNumOutputs() const;

	virtual void storeFunctionValue(const DDaceSamplePoint& pt, 
					const Array<double>& values) = 0;
	virtual void recordRunStatus(const DDaceSamplePoint& pt, 
				     DDaceRunStatus status);
	// Output.
	virtual void writeToArchive();
	virtual void writeSamples(ostream& os) const ;
	virtual void getResults(Array<DDaceSamplePoint>& pts, 
													Array<Array<double> >& funcValues) const ;
	// Some mutators
	virtual void setArchivePath(const String& archivePath);
	virtual void setArchiveName(const String& archiveName);
	virtual void setVariableNames(const Array<String>& name);
	virtual void setOutputNames(const Array<String>& name);
 protected:
	
	// All processors will get the list of sample points.
	Array<DDaceSamplePoint> pts_;
};

#endif




