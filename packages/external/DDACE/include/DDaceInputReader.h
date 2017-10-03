#ifndef DDACEINPUTREADER_H
#define DDACEINPUTREADER_H

#include "DDaceReaderBase.h"

class DDaceInputReader : public DDaceReaderBase
{
 public:
	DDaceInputReader(const Array<String>& cmds);
	DDaceInputReader(const String& filename);
	virtual ~DDaceInputReader(){;}

	virtual void getOutputNames(Array<String>& outputNames) const ;
	virtual void getArchiveFilename(String& archiveFilename) const ;
 private:
	//	void splitList(const String& big, Array<String>& list) const ;
};
	


#endif
