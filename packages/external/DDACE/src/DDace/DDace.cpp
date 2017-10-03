#include "DDace.h"
#include "DDaceClient.h"
#include "DDaceServer.h"
#include "DDaceUniproc.h"
#include "XMLObject.h"

#include <iostream>
#include <fstream>
using namespace std;

std::vector<std::string> DDace::keywords_;

DDace::DDace()
  : ptr_(0)
{
  //ExceptionBase::raise("DDace::DDace() : call to empty constructor for\n"
  //		       "DDace suggests error in creating DDace object.\n"
  //		       "This ctor should never be called, except as an\n"
  //		       "error flag.");
}

DDace::DDace(const DDaceSampler& sampler)
	: ptr_(0)
{
	ptr_ = new DDaceUniproc(sampler);
}

	
DDace::DDace(const DDaceSampler& sampler, 
		const std::vector<std::string>& varNames,
		const std::vector<std::string>& outputNames,
		const std::string& archiveFilename)
  : ptr_(0)
{

	  ptr_ = new DDaceUniproc(sampler, varNames, 
				  outputNames, archiveFilename);
}
	
DDace::DDace(const XMLObject& xmlObj)
  : ptr_(0)
{
  ptr_ = new DDaceUniproc(xmlObj);
}

	
	


bool DDace::getNextSample(DDaceSamplePoint& pt)
{
  checkPtr();
  return ptr_->getNextSample(pt);
}

void DDace::storeFunctionValue(const DDaceSamplePoint& pt,
			       const std::vector<double>& values)
{
  checkPtr();
  ptr_->storeFunctionValue(pt, values);
}

void DDace::recordRunStatus(const DDaceSamplePoint& pt,
			    DDaceRunStatus status)
{
  checkPtr();
  ptr_->recordRunStatus(pt, status);
}

void DDace::getSampler(DDaceSampler& sampler) const
{
  checkPtr();
  ptr_->getSampler(sampler);
}

void DDace::getArchiveFilename(std::string& archiveFilename) const
{
  checkPtr();
  ptr_->getArchiveFilename(archiveFilename);
}

void DDace::getVarNames(std::vector<std::string>& names) const
{
  checkPtr();
  ptr_->getVarNames(names);
}

void DDace::getOutputNames(std::vector<std::string>& names) const
{
  checkPtr();
  ptr_->getOutputNames(names);
}

int DDace::getNumOutputs() const
{
  checkPtr();
  return ptr_->getNumOutputs();
}

void DDace::getRunStatus(std::vector<DDaceRunStatus>& status) const
{
  checkPtr();
  ptr_->getRunStatus(status);
}

void DDace::writeToArchive() 
{
  checkPtr();
  ptr_->writeToArchive();
}

void DDace::writeSamples(ostream& os) const 
{
  checkPtr();
  ptr_->writeSamples(os);
}

void DDace::getResults(std::vector<DDaceSamplePoint>& pts, 
		       std::vector<std::vector<double> >& funcValues) const
{
  ptr_->getResults(pts, funcValues);
}

void DDace::setArchivePath(const std::string& archivePath)
{
  checkPtr();
  ptr_->setArchivePath(archivePath);
}

void DDace::setArchiveName(const std::string& archiveFilename)
{
  checkPtr();
  ptr_->setArchiveName(archiveFilename);
}

void DDace::setVariableNames(const std::vector<std::string>& varNames)
{
  checkPtr();
  ptr_->setVariableNames(varNames);
}

void DDace::setOutputNames(const std::vector<std::string>& outputNames)
{
  checkPtr();
  ptr_->setOutputNames(outputNames);
}

void DDace::checkPtr() const 
{
  if (ptr_==0) ExceptionBase::raise("DDace::checkPtr() failed");
}


const std::vector<std::string>& DDace::keywords()
{
  if (keywords_.length() == 0)
    {
      keywords_.resize(4);
      keywords_[0] = "DDACE";
      keywords_[1] = "BOUNDS";
      keywords_[2] = "VARIABLE";
      keywords_[3] = "RETURN";
    }
  return keywords_;
}
