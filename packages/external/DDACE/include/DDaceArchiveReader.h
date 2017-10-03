#ifndef DDACEARCHIVEREADER_H
#define DDACEARCHIVEREADER_H

#include "DDaceReaderBase.h"

class DDaceArchiveReader : public DDaceReaderBase
{
 public:
	DDaceArchiveReader(const Array<String>& cmds);
	DDaceArchiveReader(const String& filename);

	virtual void getOutputNames(Array<String>& outputNames) const ;

	virtual void getArchivedData(Array<DDaceSamplePoint>& pts,
															 Array<Array<double> >& results,
															 Array<DDaceRunStatus>& status) const ;
 private:
	void readResults(Array<DDaceSamplePoint>& pts,
									 Array<Array<double> >& results,
									 Array<DDaceRunStatus>& status,
									 int startLine) const ;
	int nSamples_;
	int nOutputs_;
};
	


#endif
