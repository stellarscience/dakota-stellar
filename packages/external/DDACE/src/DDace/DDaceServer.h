#ifndef DDACESERVER_H
#define DDACESERVER_H

#include "DDaceBool.h"
#include "DDaceMachineBase.h"
#include "DDaceSampler.h"
#include "DDaceRunStatus.h"

class DDaceServer : public DDaceMachineBase
{
 public:
	DDaceServer(const DDaceSampler& sampler,
		const std::vector<double>& lower,
		const std::vector<double>& upper,
		const std::vector<std::string>& varNames,
		const std::vector<std::string>& outputNames,
		const std::string& archiveFilename);
	DDaceServer(const DDaceSampler& sampler,
			const std::vector<double>& lower,
			const std::vector<double>& upper);

	virtual ~DDaceServer(){;}

	// sampling functions
	virtual bool getNextSample(DDaceSamplePoint& pt) ;
	virtual void storeFunctionValue(const DDaceSamplePoint& pt, 
					const std::vector<double>& values) ;
	virtual void recordRunStatus(const DDaceSamplePoint& pt, 
					 DDaceRunStatus status);
	virtual void writeToArchive();

	virtual void setArchivePath(const std::string& archivePath);
	virtual void setArchiveName(const std::string& archiveFilename);
	virtual void setVariableNames(const std::vector<std::string>& name);
	virtual void setOutputNames(const std::vector<std::string>& name);

	virtual void getResults(std::vector<DDaceSamplePoint>& pts, 
			std::vector<std::vector<double> >& funcValues) const ;
 protected:

	int stackPtr_;
	std::vector<std::vector<double> > results_;
	std::vector<DDaceRunStatus> status_;
	std::vector<double> lower_;
	std::vector<double> upper_;
	DDaceSampler sampler_;
	std::string archiveFilename_;
	std::string archivePath_;
	std::vector<std::string> varNames_;
	std::vector<std::string> outputNames_;
};

#endif
