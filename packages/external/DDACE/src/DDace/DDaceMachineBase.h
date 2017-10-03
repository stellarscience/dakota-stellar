#ifndef DDACEMACHINEBASE_H
#define DDACEMACHINEBASE_H

#include "DDaceSamplePoint.h"
#include "DDaceRunStatus.h"


/*

DDaceMachineBase defines the interface that will be used by both client
and server in a multiprocessor context, and by the single machine in
uniprocessor context. 

From this base class, we will derive DDaceClient, DDaceServer, 
and DDaceUniproc. In multiprocessor mode, client machines will construct 
DDaceClient objects while the server constructs a DDaceServer. 

Server-specific methods, such as I/O methods, do nothing when 
called on the base class.
 */



// Some human-readable message tags and status flags
enum MsgID {RequestTag, IndexTag, ResultIndexTag, ResultValueTag};
enum ProcStatus {Available, Busy, AskedForData, 
		 AskedForRequest, Faulty, Done};


class DDaceMachineBase
{
 public:
	DDaceMachineBase() : pts_(0) {;}
	virtual ~DDaceMachineBase(){;}

	// sampling and result storage methods, to be called in main sampling loop.
	virtual bool getNextSample(DDaceSamplePoint& pt) = 0;
	virtual void storeFunctionValue(const DDaceSamplePoint& pt, 
					const std::vector<double>& values) = 0;
	virtual void recordRunStatus(const DDaceSamplePoint& pt, 
					 DDaceRunStatus status);
	// Output.
	virtual void writeToArchive();
	virtual void writeSamples(ostream& os) const ;
	virtual void getResults(std::vector<DDaceSamplePoint>& pts, 
			std::vector<std::vector<double> >& funcValues) const ;
	// Some mutators
	virtual void setArchivePath(const std::string& archivePath);
	virtual void setArchiveName(const std::string& archiveName);
	virtual void setVariableNames(const std::vector<std::string>& name);
	virtual void setOutputNames(const std::vector<std::string>& name);
 protected:
	
	// All processors will get the list of sample points.
	std::vector<DDaceSamplePoint> pts_;
};

#endif




