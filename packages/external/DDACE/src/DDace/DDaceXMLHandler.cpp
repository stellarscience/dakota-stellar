#include "DDaceXMLHandler.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <hash_set>
#include "UniformDistribution.h"
#include "NormalDistribution.h"
#include "DDaceLHSampler.h"
#include "DDaceRandomSampler.h"
#include "DDaceOASampler.h"
#include "DDaceFactorialSampler.h"
#include "DDaceUserInputSampler.h"
#include "SmartPtr.h"
#include "DDaceRunStatus.h"

using namespace std;

struct eqstr
{
  bool operator()(const std::string s1, const string s2) const
  {
    return s1.compare(s2) == 0;
  }
};

std::hash_set<std::string, hash<string>, eqstr> DDaceXMLHandler::samplerTypes_;


DDaceXMLHandler::DDaceXMLHandler()
  : samplerName_("null"), samplerParams_(), varNames_(), varDist_(),
    outputNames_(), archiveFilename_("ddaceArchive"), currentTag_(-1), 
    charStr_(), pts_(), results_(), status_(),isSamplePt_(false), 
    isSampleResult_(false), isArchiveFile_(false), getPtsComplete_(false)
{
  try
    {
      if (samplerTypes_.size()==0)
	{
	  samplerTypes_.insert("Random", "blah");
	  samplerTypes_.insert("LatinHypercube", "blah");
	  samplerTypes_.insert("OrthogonalArray", "blah");
	  samplerTypes_.insert("Factorial", "blah");
	  samplerTypes_.insert("UserInputSampler", "blah");

	}
    }
  catch(std::exception& e)
    {
	std::cerr << e.what();
    }
}



DDace DDaceXMLHandler::createObject()  
{
  try
    {
      DDaceSampler sampler = createSampler();
      DDace ddaceObj(sampler, varNames_, outputNames_, archiveFilename_);

      // Through the process of creating ddaceObj, the pts_ have been
      // created:

      getPtsComplete_ = true;
      // if we are reading from an archive file, then we need to 
      // store the pts_, results_, and status_ values in the DDace
      // object.  Otherwise, these will be given default values
      // with the above creation of the ddaceObj, and will be
      // filled in as the runs are completed.
      if(isArchiveFile_)
	{
	  if(getPtsComplete_)
	    {
	      for (int i=0; i<pts_.length(); i++)
		{
		  ddaceObj.storeFunctionValue(pts_[i], results_[i]);
		  ddaceObj.recordRunStatus(pts_[i], status_[i]);
		}
	    }
	}
      return ddaceObj;
    }
  catch(std::exception& e)
    {
	std::cerr << "An error occurred: " << e.what();
    }
  return DDace(); // -Wall
}

DDaceSampler DDaceXMLHandler::createSampler() const
{
  try
    {
      DDaceSampler sampler;

      if (samplerParams_.containsKey("seed"))
	{
	  int seed = atoi(samplerParams_.get("seed"));
	  DistributionBase::setSeed(seed);
	}
      
      if (samplerName_ == "Random")
	{
	  int nSamples = atoi(samplerParams_.get("samples"));
	  sampler = DDaceRandomSampler(nSamples, varDist_);
	}
      else if (samplerName_ == "LatinHypercube")
	{
	  int nSamples = atoi(samplerParams_.get("samples"));
	  int nReps = atoi(samplerParams_.get("replications"));
	  bool noise = false;
	  if (samplerParams_.containsKey("perturb"))
	    {
	      noise = atob(samplerParams_.get("perturb"));
	    }
	  sampler = DDaceLHSampler(nSamples, nReps, noise, varDist_);
	}
      else if (samplerName_ == "OrthogonalArray")
	{
	  int nSamples = atoi(samplerParams_.get("samples"));
	  bool noise = false;
	  if (samplerParams_.containsKey("perturb"))
	    {
	      noise = atob(samplerParams_.get("perturb"));
	    }
	  sampler = DDaceOASampler(nSamples, noise, varDist_);
	}
      else if (samplerName_ == "Factorial")
	{
	  int nSamples = atoi(samplerParams_.get("samples"));
	  int nReps = atoi(samplerParams_.get("symbols"));
	  bool noise = false;
	  if (samplerParams_.containsKey("perturb"))
	    {
	      noise = atob(samplerParams_.get("perturb"));
	    }
	  sampler = DDaceFactorialSampler(nSamples, nReps, noise, varDist_);
	}
      else if (samplerName_ == "UserInputSampler")
	{
	  std::string filename = samplerParams_.get("filename");
	  sampler = DDaceUserInputSampler(filename);
	}
      else
	{
	  throw runtime_error("unrecognized sampler type in DDaceXMLHandler::createSampler()");
	}
      return sampler;
    }
  catch(std::exception& e)
    {
	cerr << "In DDaceXMLHandler::createSampler(): " << e.what() << endl;
    }

  return DDaceSampler();
}

Array<std::string> DDaceXMLHandler::getVarNames() const
{
  return varNames_;
}

Array<std::string> DDaceXMLHandler::getOutputNames() const
{
  return outputNames_;
}

void DDaceXMLHandler::startElement(const std::string& name, 
				   const Hashtable& attributes)
{
  try
    {
      if (samplerTypes_.containsKey(name))
	{
	  samplerParams_ = attributes;
	  samplerName_ = name;
	  
	  int samples = atoi(attributes.get("samples"));
	  pts_.resize(samples);
	  results_.resize(samples);
	  status_.resize(samples);
	  // give default values to elements of status_
	  for(int j=0; j < samples; j++)
	    status_[j] = DDaceRunNotStarted;
	}
      else if (name == "Variable")
	{
	  varNames_.append(attributes.get("name"));
	  
	  if (attributes.containsKey("distribution") 
	      && attributes.get("distribution") != "uniform")
	    {
	      std::string distType = attributes.get("distribution");
	      if (distType == "normal")
		{
		  if (attributes.containsKey("mean") 
		      && attributes.containsKey("sigma"))
		    {
		      double mean = atof(attributes.get("mean"));
		      double sigma = atof(attributes.get("sigma"));
		      double nDev = atof(attributes.get("cutoff"));
		      varDist_.append(NormalDistribution(Mean(mean), 
							 StdDeviation(sigma),
							 nDev));
		    }
		  else if (attributes.containsKey("upper") && attributes.containsKey("lower"))
		    {
		      double lower = atof(attributes.get("lower"));
		      double upper = atof(attributes.get("upper"));
		      if (attributes.containsKey("nDev"))
			{
			  double nDev = atof(attributes.get("cutoff"));
			  varDist_.append(NormalDistribution(lower, upper, nDev));
			}
		      else
			{
			  varDist_.append(NormalDistribution(lower, upper));
			}
		    }
		  else
		    {
		      ExceptionBase::raise("invalid xml specification of normally distributed variable");
		    }
		}
	      else
		{
		  ExceptionBase::raise("unrecognized xml specification of distribution type");
		}
	    }			
	  else if ( samplerName_ == "UserInputSampler" )
	    {
	      // no distribution information required for UserInputSampler
	    }	      
	  else /* uniform distribution */
	    {
	      double lower = atof(attributes.get("lower"));
	      double upper = atof(attributes.get("upper"));
	      varDist_.append(UniformDistribution(lower, upper));
	      
	    }
	}
      else if (name == "Archive")
	{
	  archiveFilename_ = attributes.get("name");
	}
      else if (name == "Output")
	{
	  outputNames_.append(attributes.get("name"));
	}
      else if (name == "Sample")
	{
	  currentTag_ = atoi(attributes.get("tag"));
	}
      else if (name == "SamplePoint")
	{
	  isSamplePt_ = true;
	  charStr_ = "";
	}
      else if (name == "SampleResult")
	{
	  isArchiveFile_ = true;
	  if (currentTag_ == -1)
	    ExceptionBase::raise("Index tag corresponding to current "
				 "sample is -1.");
	  String statusString = attributes.get("status");
	  if (statusString == "Run OK")
	    status_[currentTag_] = DDaceRunOK;
	  else if (statusString == "Run failed")
	    status_[currentTag_] = DDaceRunFailed;
	  else if (statusString == "Post processing failed")
	    status_[currentTag_] = DDacePostProcFailed;
	  else if (statusString == "Run not started")
	    status_[currentTag_] = DDaceRunNotStarted;
	  else if (statusString == "Run pending")
	    status_[currentTag_] = DDaceRunPending;
	  else 
	    ExceptionBase::raise("unrecognized sample status inside "
				 "SampleResult for "
				 "DDaceXMLHandler::startElement()");
	  isSampleResult_ = true;
	  charStr_= "";
	}
    }
  catch(ExceptionBase& e)
    {
      e.trace("in DDaceXMLHandler::startElement( tag = " + name + ")");
    }
}

void DDaceXMLHandler::characters(const String& chars, 
				 const unsigned int length)
{

  if(isSamplePt_ || isSampleResult_)
    charStr_ += chars;
}

void DDaceXMLHandler::endElement(const String& name )
{
  if (name == "SamplePoint")
    {
      if (currentTag_ == -1)
	ExceptionBase::raise("Index tag corresponding to current "
			     "sample is -1.");
      Array<double> tempPt;
      Array<String> charStrLines = StrUtils::splitIntoLines(charStr_);
      tempPt.resize(charStrLines.length());
      for(int i=0; i < charStrLines.length(); i++)
	tempPt[i] = atof(charStrLines[i]);
      DDaceSamplePoint point(currentTag_, tempPt);
      pts_[currentTag_] = point;

      // reset these variables as we exit the SamplePoint block...
      isSamplePt_ = false;
      charStr_="";
    }
  if (name == "SampleResult")
    {
      if (currentTag_ == -1)
				ExceptionBase::raise("Index tag corresponding to current "
														 "sample is -1.");
      Array<String> charStrLines = StrUtils::splitIntoLines(charStr_);
      results_[currentTag_].resize(charStrLines.length());
      for(int i=0; i < charStrLines.length(); i++)
	results_[currentTag_][i] = atof(charStrLines[i]);

      // reset these variables as we exit the SamplePoint block...
      isSampleResult_ = false;
      charStr_="";
    }
  if (name == "Results")
    {
      // if we've reached the end of the Results section, and
      // the currentTag_ is one less than the number of pts_,
      // we've read in all (and the appropriate number) of the 
      // pts_ values from the archive file.
      if(currentTag_ == (pts_.length() - 1))
	getPtsComplete_ = true;
    }
}



