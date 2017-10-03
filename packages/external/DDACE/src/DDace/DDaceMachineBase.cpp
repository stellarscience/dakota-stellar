#include "DDaceMachineBase.h"


void DDaceMachineBase::writeSamples(ostream& os) const
{
	for (int i=0; i<pts_.length(); i++)
		{
			os << i << " " << pts_[i] << endl;
		}
}

void DDaceMachineBase::recordRunStatus(const DDaceSamplePoint& /* pt */, 
					 DDaceRunStatus /* status */)
{
	return; // client does nothing
}
	
void DDaceMachineBase::writeToArchive() 
{
	return; // client does nothing
}

void DDaceMachineBase::setArchivePath(const std::string& /* archiveFilename */)
{
	return; // client does nothing
}

void DDaceMachineBase::setArchiveName(const std::string& /* archiveFilename */)
{
	return; // client does nothing
}

void DDaceMachineBase::setVariableNames(const std::vector<std::string>&)
{
	return; // client does nothing
}

void DDaceMachineBase::setOutputNames(const std::vector<std::string>&)
{
	return; // client does nothing
}

void DDaceMachineBase::getResults(std::vector<DDaceSamplePoint>&, 
				std::vector<std::vector<double> >& ) const
{
	return; // does nothing
}
		
void DDaceMachineBase::getSampler(DDaceSampler& sampler) const
{
  throw runtime_error("DDaceMachineBase::getSampler() : cannot call this "
		       "method directly.\nOnly call through DDaceServer or"
		       "DDaceUniproc.");
  return; //does nothing
}

void DDaceMachineBase::getRunStatus(std::vector<DDaceRunStatus>& status) const
{
  throw runtime_error("DDaceMachineBase::getRunStatus() : cannot call "
		       "this method directly.\nOnly call through "
		       "DDaceServer or DDaceUniproc.");
  return; //does nothing
}

void DDaceMachineBase::getArchiveFilename(std::string& archiveFilename) const
{
  throw runtime_error("DDaceMachineBase::getArchiveFilename() : cannot "
		       "call this method directly.\nOnly call through "
		       "DDaceServer or DDaceUniproc.");
  return; //does nothing
}
  
void DDaceMachineBase::getVarNames(std::vector<std::string>& names) const
{
  throw runtime_error("DDaceMachineBase::getVarNames() : cannot call "
		       "this method\ndirectly. Only call through DDaceServer "
		       "or DDaceUniproc.");
}


void DDaceMachineBase::getOutputNames(std::vector<std::string>& names) const
{
  throw runtime_error("DDaceMachineBase::getOutputNames() : cannot call "
		       "this method\ndirectly. Only call through DDaceServer "
		       "or DDaceUniproc.");
}

int DDaceMachineBase::getNumOutputs() const
{
  throw runtime_error("DDaceMachineBase::getNumOutputs() : cannot call "
		       "this method\ndirectly. Only call through DDaceServer "
		       "or DDaceUniproc.");
  return 0;
}
