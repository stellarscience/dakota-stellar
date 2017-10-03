#include "DDaceAnalyzerBase.h"



DDaceAnalyzerBase::DDaceAnalyzerBase(const String& archiveFile, const String& outputName)
  : archiveFile_(archiveFile), 
  outputNames_(1), 
  outputData_(0), 
  ddaceObj_(), 
  outputIndex_(-1), 
  nSamples_(0)
{
  try {
    DDaceXMLReader reader(archiveFile);
    
    //create the DDace object from the information in the archive file
    ddaceObj_ = reader.createObject();

    Array<String> outputNameList(0);
    ddaceObj_.getOutputNames(outputNameList);
    
    // go through the list of output names extracted from the DDace
    // object, and find out if the user selected output name matches
    // with one for the list.  If so, assign the outputName to the
    // data member outputName_, and keep track of the index (relative
    // to the list of output names in the DDace object) in outputIndex_.
    int i = 0;
    for(i=0; i < outputNameList.length(); i++)
      {
	if(outputNameList[i] == outputName)
	  {
	    outputIndex_ = i;
	    outputNames_[0] = outputName;
	    break;
	  }
      }
    
    // if we didn't match the user specified output with something in
    // the list, then return an error.
    if(outputIndex_ < 0) ExceptionBase::raise("DDaceAnalyzerBase() ctor : "
					      "outputName provided does "
					      "not\nmatch any output names"
					      " in the ddaceObj.");
    
    //get the number of samples in the DDace object.
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    nSamples_ = sampler.nSamples();

    // get the data corresponding to outputName_ from ddaceObj, and put 
    // it in the outputData_ Array.
    Array<DDaceSamplePoint> tmpPts;
    Array<Array<double> > tmpOutput;
    ddaceObj_.getResults(tmpPts, tmpOutput);
    outputData_.resize(nSamples_);
    int j = 0;
    for(j = 0; j < nSamples_; j++)
      {
	outputData_[j] = tmpOutput[j][outputIndex_];
      }
  }
  catch(ExceptionBase& e)
    { e.trace("in DDaceAnalyzerBase::ctor()"); }
}


DDaceAnalyzerBase::DDaceAnalyzerBase(const String& archiveFile)
  : archiveFile_(archiveFile), 
  outputNames_(0), 
  outputData_(0), 
  ddaceObj_(), 
  outputIndex_(-1), 
  nSamples_(0)
{
  try{
    DDaceXMLReader reader(archiveFile);
    ddaceObj_ = reader.createObject();
  
    ddaceObj_.getOutputNames(outputNames_);
  
    //get the number of samples in the DDace object.
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    nSamples_ = sampler.nSamples();
 
    //Note: the outputData_ will be gathered straight form ddaceObj_ when
    //it is needed.
  }
  catch(ExceptionBase &e)
    { e.trace("in DDaceAnalyzerBase::ctor()");}
  
} 



DDaceAnalyzerBase::DDaceAnalyzerBase(const DDace& ddaceObj, const String& outputName)
  : archiveFile_(""), 
  outputNames_(1), 
  outputData_(0), 
  ddaceObj_(ddaceObj),
  outputIndex_(-1), 
  nSamples_(0)
{

  //Note: we created the ddaceObj_ member of DDaceAnalyzerBase directly
  // from a DDace object, instead of from an archive file.  So there
  // will be an empty string stored for the archiveFile_ name.

  try{
    Array<String> outputNameList(0);
    ddaceObj_.getOutputNames(outputNameList);
  
    // go through the list of output names extracted from the DDace
    // object, and find out if the user selected output name matches
    // with one for the list.  If so, assign the outputName to the
    // data member outputName_, and keep track of the index (relative
    // to the list of output names in the DDace object) in outputIndex_.
    for(int i=0; i < outputNameList.length(); i++)
      {
	if(outputNameList[i] == outputName)
	  {
	    outputIndex_ = i;
	    outputNames_[0] = outputName;
	    break;
	  }
      }
    if (outputIndex_ < 0) ExceptionBase::raise("DDaceAnalyzerBase() ctor : "
					       "outputName provided does "
					       "not\nmatch any output names"
					       " in the ddaceObj.");
    
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    int nSamples_ = sampler.nSamples();
    
    Array<DDaceSamplePoint> tmpPts;
    Array<Array<double> > tmpOutput;
    ddaceObj_.getResults(tmpPts, tmpOutput);
    outputData_.resize(nSamples_);
    int j = 0;    
    // store data corresponding to outputName_ in outputData_.
    for(j = 0; j < nSamples_; j++)
      {
	outputData_[j] = tmpOutput[j][outputIndex_];
      }
  }
  catch(ExceptionBase& e)
    { e.trace("in DDaceAnalyzerBase() ctor."); }
}

DDaceAnalyzerBase::DDaceAnalyzerBase(const DDace& ddaceObj)
    : archiveFile_(""), outputNames_(0),  outputData_(0), ddaceObj_(ddaceObj),
      outputIndex_(-1), nSamples_(0)
    { 
      try
	{
	  ddaceObj_.getOutputNames(outputNames_);
	  DDaceSampler sampler;
	  ddaceObj_.getSampler(sampler);
	  nSamples_ = sampler.nSamples();
	}
      catch(ExceptionBase& e)
	{ e.trace("in DDaceAnalyzerBase() ctor()."); }
    }


