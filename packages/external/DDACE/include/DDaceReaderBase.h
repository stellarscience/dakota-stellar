#ifndef DDACEREADERBASE_H
#define DDACEREADERBASE_H

#include "Array.h"
#include "String.h"
#include "DDaceSampler.h"
#include "DDaceSamplePoint.h"
#include "DDaceRunStatus.h"

class DDaceReaderBase
{
 public:
	DDaceReaderBase(const Array<String>& cmds);
	DDaceReaderBase(const String& filename);
	virtual ~DDaceReaderBase(){;}

	virtual void getBounds(Array<double>& lower, Array<double>& upper) const ;
	virtual void getSampler(DDaceSampler& sampler) const ;
	virtual void getVarNames(Array<String>& varNames) const ;
	virtual void getOutputNames(Array<String>& outputNames) const = 0 ;
	virtual void getArchiveFilename(String& archiveFilename) const ;
	
	virtual void getArchivedData(Array<DDaceSamplePoint>& pts,
															 Array<Array<double> >& results,
															 Array<DDaceRunStatus>& status) const ;
 protected:
	bool lookup(const String& keyword, String& cmd) const ;
	int countInputs() const ;

	int nInputs_;
	Array<String> cmds_;
	Array<Array<String> > tokens_;
};
	


#endif
