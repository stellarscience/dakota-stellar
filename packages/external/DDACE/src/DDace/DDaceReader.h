#ifndef DDACEREADER_H
#define DDACEREADER_H

#include "DDaceReaderBase.h"
#include "SmartPtr.h"

class DDaceReader
{
 public:
	DDaceReader(DDaceReaderBase* ptr) : ptr_(ptr) {;}

	void getBounds(Array<double>& lower, Array<double>& upper) const 
		{ptr_->getBounds(lower, upper);}

	void getSampler(DDaceSampler& sampler) const 
		{ptr_->getSampler(sampler);}

	void getVarNames(Array<String>& varNames) const 
		{ptr_->getVarNames(varNames);}
			
	void getOutputNames(Array<String>& outputNames) const 
		{ptr_->getOutputNames(outputNames);}

	void getArchiveFilename(String& archiveFilename) const 
		{ptr_->getArchiveFilename(archiveFilename);}
	
	void getArchivedData(Array<DDaceSamplePoint>& pts,
											 Array<Array<double> >& results,
											 Array<DDaceRunStatus>& status) const 
		{ptr_->getArchivedData(pts, results, status);}

 protected:
	SmartPtr<DDaceReaderBase> ptr_;
};
	


#endif
