#include "DDaceServer.h"
#include "ArrayComm.h"
#include <iostream>
#include <fstream>
using namespace std;
#include "System.h"
#include "PComm.h"
#include "DDace.h"
#include "SmartPtr.h"
#include "XMLObject.h"
#include "DDaceSampler.h"
#include "UniformDistribution.h"
#include "NormalDistribution.h"
#include "DDaceLHSampler.h"
#include "DDaceRandomSampler.h"
#include "DDaceOASampler.h"
#include "DDaceFactorialSampler.h"
#include "DDaceUserInputSampler.h"

// ------------------------------------------------------------------------
// Ctors for DDaceServer. DDaceServer holds all problem data, names, etc.
// 
// The DDaceServer object is constructed from existing objects rather
// than by parsing an input file; we assume that the input file (if any)
// has been parsed at a higher level.
// ------------------------------------------------------------------------

DDaceServer::DDaceServer(const DDaceSampler& sampler, 
			 const Array<String>& varNames,
			 const Array<String>& outputNames, 
			 const String& archiveFilename)
  : DDaceMachineBase(),
    stackPtr_(0),
    results_(0), 
    status_(0),
    sampler_(sampler),
    archiveFilename_(archiveFilename),
    varNames_(varNames),
    outputNames_(outputNames)
{
  sampler_.getSamples(pts_);
  results_.resize(pts_.length());
  status_.resize(pts_.length());
  for (int i=0; i<status_.length(); i++) status_[i] = DDaceRunNotStarted;
  pts_.bcast(0, DDace::getComm());
}



DDaceServer::DDaceServer(const DDaceSampler& sampler) 
  : DDaceMachineBase(),
    stackPtr_(0),
    results_(0), 
    status_(0),
    sampler_(sampler),
    archiveFilename_("ddaceArchive"),
    varNames_(sampler.nInputs()),
    outputNames_(0)
{
  int i;
  sampler_.getSamples(pts_);
  results_.resize(pts_.length());
  status_.resize(pts_.length());
  for (i=0; i<status_.length(); i++) status_[i] = DDaceRunNotStarted;
  // sould we initialize the values in results_?  

  pts_.bcast(0, DDace::getComm());
  
  // make up some dummy variable names if the user hasn't given any.
  for (i=0; i<varNames_.length(); i++)
    {
      char tmp[20];
      sprintf(tmp, "var%d", i);
      varNames_[i] = tmp;
    }
}

DDaceServer::DDaceServer(const XMLObject& xmlObj)
  : DDaceMachineBase(),
    stackPtr_(0),
    results_(0), 
    status_(0),
    sampler_(),
    archiveFilename_(""),
    varNames_(0),
    outputNames_(0)
{

try
{
  int nc = 0;

  int nsamples = 0, nreplications = 0, nsymbols = 0;

  bool noise;
  std::string userInputFile, archiveFile;
  std::vector<std::string> varNames(0);
  std::vector<std::string> outNames(0);
  std::vector<Distribution> distributions(0);
  std::string samplerName;

  if(xmlObj.getTag() == "DDace")
  {
 	nc = xmlObj.numChildren(); // should be 4 -- sampler, archive name,
  	 // input variable names, output names
 	for(int i=0; i < nc; i++)
 	{
 	  XMLObject xKid = xmlObj.getChild(i);

 	// Settings associated with the samplers:
	  if((xKid.getTag() == "Random") ||
   	(xKid.getTag() == "Factorial") ||
   	(xKid.getTag() == "OrthogonalArray") ||
   	(xKid.getTag() == "LatinHypercube") || 
   	(xKid.getTag() == "UserInputSampler"))
  	  {
  	    //keep track of the sampler name....
  	    samplerName = xKid.getTag();
  	
  	    if((xKid.getTag() == "Random") ||
  		 (xKid.getTag() == "Factorial") ||
  		 (xKid.getTag() == "OrthogonalArray") ||
  		 (xKid.getTag() == "LatinHypercube"))
  		{
  			nsamples = xKid.getAttribute("samples").atoi();
  		}
  	
  	if((xKid.getTag() == "Factorial") ||
  		 (xKid.getTag() == "LatinHypercube") ||
  		 (xKid.getTag() == "OrthogonalArray"))
  		{
  			noise = xKid.getAttribute("noise").atob();
  		}
  	
  	if(xKid.getTag() == "Factorial")
  		nsymbols = xKid.getAttribute("symbols").atoi();
  	else if(xKid.getTag() == "LatinHypercube")
  		nreplications = xKid.getAttribute("replications").atoi();
  	else if(xKid.getTag() == "UserInputSampler")
  		{
  			XMLObject tmpKid = xKid.getChild(0);
  			userInputFile = tmpKid.getAttribute("filename");
  		}
  }
 			// Settings associated with the archive file:
 			else if(xKid.getTag() == "Archive")
  archiveFile = xKid.getAttribute("name");

 			// Array: could be input variable names (varNames) or 
 			// output variable names (outVarNames)....
 			else if(xKid.getTag() == "Array")
  {
  	// set up an Array of input variable names....
  	//if(xKid.getAttribute("name") == "inputVariables")
  	if(xKid.getAttribute("name") == "inputVariables")
  		{
  			int nc1 = xKid.numChildren();
  			if(nc1 == 0)
   ExceptionBase::raise("DDaceServer XML ctor : number "
     		 "of input variables must be\n"
     		 "greater than zero.");

  			varNames.resize(nc1);
  			distributions.resize(nc1);
  			for(int j = 0; j < nc1; j++)
   {
   	XMLObject varKid = xKid.getChild(j);
   	if(varKid.getTag() == "DDaceVariable")
   		{
   			// get variable name
   			varNames[j] = varKid.getAttribute("variable");

   			// get distribution information
   			XMLObject distKid = varKid.getChild(0);
   			if(distKid.getTag() == "UniformDistribution")
    {
    	double lwr = distKid.getAttribute("lower").atof();
    	double upr = distKid.getAttribute("upper").atof();
    		
    	distributions[j] = UniformDistribution(lwr, upr);
    	//distributions.append(UniformDistribution(lwr, upr));
    }
   			else if(distKid.getTag() == "NormalDistribution")
    {
    	double mean = distKid.getAttribute("mean").atof();
    	double sigma = distKid.getAttribute("sigma").atof();
    	
    	//distributions.append(NormalDistribution(mean, sigma));
    	distributions[j] = NormalDistribution(mean, sigma);
    }
   			else
    ExceptionBase::raise("DDaceServer XML ctor : "
      		 "only support for Normal\n"
      		 "and Uniform distributions"
      		 " in current version of "
      		 "DDace.");
   		}
   	else
   		ExceptionBase::raise("DDaceServer XML ctor : "
       "Array with name "
       "inputVariables\ndoes not"
       "contain DDaceVariables...");
   }
  		}
  	else if(xKid.getAttribute("name") == "outVarNames")
  		{
  			int nc2 = xKid.numChildren();
  			if(nc2 == 0)
   ExceptionBase::raise("DDaceServer XML ctor : number "
     		 "of output variables must be\n"
     		 "greater than zero.");
  			outNames.resize(nc2);
  			for(int k=0; k < nc2; k++)
   {
   	XMLObject outKid = xKid.getChild(k);
   	outNames[k] = outKid.getAttribute("value");
   }
  		}
  }
 		}
 }
			else
 {
 	ExceptionBase::raise("DDaceServer XML ctor : XML object not of"
   			 "type DDace.");
 }


			// we should have all the pieces to build a DDace object, so build it:

			// First, build the sampler:
			if(samplerName == "Random")
 sampler_ = DDaceRandomSampler(nsamples, distributions);
			else if(samplerName == "Factorial")
 sampler_ = DDaceFactorialSampler(nsamples,nsymbols, 
      noise, distributions);
			else if(samplerName == "LatinHypercube")
 sampler_ = DDaceLHSampler(nsamples, nreplications, noise, distributions);
			else if(samplerName == "OrthogonalArray")
 sampler_ = DDaceOASampler(nsamples, noise, distributions);
			else if(samplerName == "UserInputSampler")
 sampler_ = DDaceUserInputSampler(userInputFile);
			else
 ExceptionBase::raise("DDaceServer XML ctor : sampler name not "
   		 "recognized.");
			
			sampler_.getSamples(pts_);
			results_.resize(pts_.length());
			status_.resize(pts_.length());
			for (int i=0; i<status_.length(); i++) status_[i] = DDaceRunNotStarted;

			pts_.bcast(0, DDace::getComm());
  
			varNames_.resize(sampler_.nInputs());
			if(varNames.length() != varNames_.length())
 ExceptionBase::raise("DDaceServer XML ctor : number of variables read "
   		 "from XMLObject \nis not the same as number of "
   		 "variables expected from XMLObject.");
			else
 varNames_ = varNames;
			
			outputNames_.resize(outNames.length());
			outputNames_ = outNames;
			
			if(archiveFile != "")
 archiveFilename_ = archiveFile;
			else
 archiveFilename_ = "ddaceArchive";
		}
	catch(ExceptionBase& e)
		{
			e.trace("in DDaceServer XML ctor");
		}
}
	
	
// ------------------------------------------------------------------------
// On the server, getNextSample doles out sample points to the processors
// and then waits for data to be returned. 
// It loops until everything is done, and then returns false to terminate
// the calling loop.
// ------------------------------------------------------------------------	
	
	
bool DDaceServer::getNextSample(DDaceSamplePoint& /* pt */)
{
  if (pts_.length()==0) return false;
  
  int nProc = DDace::getComm().getNProc();
  
  Array<int> dataIndex(nProc);
  int request;
  
  // status variables: 
  // bool done indicates whether all processors are finished
  bool done = false;
  
  // masterDone indicates whether all sample points have been handed out.
  // After the master is done, we need to keep running to collect data
  // from those processors that are still working.
  bool masterDone = false;
  
  // processorNotified indicates whether processor #i has been told that
  // all samples are done.
  Array<bool> processorNotified(nProc, false);
  
  // procStatus indicated the current state of processor #i. 
  // "Available" means that it is ready to send a request for a new sample pt.
  // "AskedForRequest" means that the master is waiting for a request from it.
  // "AskedForData" means that the master is waiting for data from it.
  // "Busy" means that it is running a function evaluation.
  // "Faulty" means that we've detected an error in a transmission from it.
  Array<ProcStatus> procStatus(nProc, Available);
  
  while (!done)
    {
      // look for requests from all available processors.
      // look for results from all busy processors.
      for (int i=1; i<nProc; i++)
	{
	  cerr << "Server: Check for Proc Status" << endl;
	  if (procStatus[i] == Available)
	    {
	      cerr << "Server: Proc Available; receive request" << endl;
	      PMachine::irecv(&request, 1, PMachine::INT, (int) RequestTag, i, 
			      DDace::getComm());
	      procStatus[i] = AskedForRequest;
	      cerr << "Server: Proc Asked for Request" << endl;
	    }
	  else if (procStatus[i] == Busy)
	    {
	      dataIndex[i] = -1;
	      PMachine::irecv(&(dataIndex[i]), 1, PMachine::INT, 
			      (int) ResultIndexTag, i, DDace::getComm());
	      procStatus[i] = AskedForData;
	    }
	}
      
      // pick up next incoming transmission
      
      cerr << "Server: waitAny pick up next incoming transmission" << endl;
      int srcID = PMachine::waitAny(DDace::getComm());
      cerr << "Server: after waitAny srcID = " << srcID << endl;
      
      // check for bad processor
      if (srcID < 0)
	{
	  procStatus[-srcID] = Faulty;
	  break;
	}
      
      
      // if we've been waiting for a request for a sample pt, send it
      // the index for the next sample.
      if (procStatus[srcID] == AskedForRequest)
	{
	  int sendIndex = stackPtr_;
	  procStatus[srcID] = Busy;
	  if (stackPtr_ >= pts_.length())
	    {
	      masterDone = true;
	      sendIndex = -1;
	      processorNotified[srcID] = true;
	      procStatus[srcID] = Done;
	    }
	  cerr << "Server: send sample point to p" << srcID << endl;

	  PMachine::send(&sendIndex, 1, PMachine::INT, 
			 (int) IndexTag, srcID, DDace::getComm());
	  stackPtr_++;
	} 
      // if we've been waiting for data, receive the next data point.
      else if (procStatus[srcID] == AskedForData)
	{
	  int status;
	  cerr << "Server: waiting for data" << endl;
	  PMachine::recv(&status, 1, PMachine::INT, (int) ResultIndexTag, 
			 srcID, DDace::getComm());
	  cerr << "Server: received data from p " << srcID << endl;
	  status_[dataIndex[srcID]] = (DDaceRunStatus) status;
	  if (status_[dataIndex[srcID]] == DDaceRunOK)
	    {
	      Array<double> values;
	      ArrayComm::recv(values, (int) ResultValueTag, srcID, 
			      DDace::getComm());
	      results_[dataIndex[srcID]] = values;
	    }
	  procStatus[srcID] = Available;
	  writeToArchive();
	}
      else 
	{
	  done = true;
	  break;
	}
      
      // if all sample points have been sent, keep going until all processors
      // have finished.
      if (masterDone)
	{
	  done = true;
	  for (int i=1; i<nProc; i++)
	    {
	      if (!processorNotified[i]) 
		{
		  done = false;
		  break;
		}
	    }
	}
    }
  return false;
	
}


void DDaceServer::getSampler(DDaceSampler& sampler) const
{
  sampler = sampler_;
}

void DDaceServer::getRunStatus(Array<DDaceRunStatus>& status) const
{
  status = status_;
}

void DDaceServer::getArchiveFilename(String& archiveFilename) const
{
  archiveFilename = archiveFilename_;
}

void DDaceServer::getVarNames(Array<String>& names) const
{
  names = varNames_;
}

void DDaceServer::getOutputNames(Array<String>& names) const
{
  names = outputNames_;
}

int DDaceServer::getNumOutputs() const
{
  return outputNames_.length();
}

// ------------------------------------------------------------------------
// The server's copy of storeFunctionValue is a dummy. The actual
// storing of function values on the server side is done in the getSample
// function call. Wierd but true.
// ------------------------------------------------------------------------

void DDaceServer::storeFunctionValue(const DDaceSamplePoint& /* pt */,
      const Array<double>& /* values */)
{
  ExceptionBase::raise("internal error: DDaceServer::storeFunctionValue should never be called");
}

void DDaceServer::recordRunStatus(const DDaceSamplePoint& pt, 
   DDaceRunStatus status)
{
  ExceptionBase::raise("internal error: DDaceServer::recordRunStatus should never be called");
}


// ------------------------------------------------------------------------
// Write all information to the archive file. This is intended to be a 
// dump that can be used for restart.
// ------------------------------------------------------------------------

void DDaceServer::writeToArchive()
{
  int i;
  int j;
  String file;
  
  if (archivePath_.length() > 0)
    {
      file = archivePath_ + "/" + archiveFilename_;
    }
  else
    {
      file = archiveFilename_;
    }
  ofstream of(file.cString());
  
  of << "<DDaceArchive date=\"" << System::date() << "\">" << endl;
  of << sampler_ << endl;
  
  if(sampler_.typeName() == "DDaceUserInputSampler")
    {
      for (i=0; i<sampler_.nInputs(); i++)
	{ 
	  of << "<Variable name=\"" << varNames_[i] 
	     << "\" lower=\"" << sampler_.lowerBounds()[i]
	     << "\" upper=\"" << sampler_.upperBounds()[i]
	     << "\"/>" << endl;
	}
    }
  else
    {
      for (i=0; i<sampler_.nInputs(); i++)
	{
	  of << "<Variable name=\"" << varNames_[i] << "\" ";
	  sampler_.dist()[i].printAttributes(of);
	  of << "/>" << endl;
	}
    }
  for (i=0; i<outputNames_.length(); i++)
    {
      of << "<Output name=\"" << outputNames_[i] << "\"/>" << endl;
    }
  of << "<Results>" << endl;
  for (i=0; i<pts_.length(); i++)
    {
      of << "<Sample tag=\"" << i << "\">" << endl;
      of << "<SamplePoint> " << endl;
      for (j=0; j<pts_[i].length(); j++)
	{
	  of << '\t' << pts_[i][j] << endl;
	}
      of << "</SamplePoint>" << endl;
      of << "<SampleResult status=\"";
      
      if (status_[i]==DDaceRunOK)
	{
	  of << "Run OK\">" << endl;  
	  for (j=0; j<results_[i].length(); j++)
	    {
	      of << '\t' << results_[i][j] << endl;
	    }
	  of << "</SampleResult>" << endl;
	}
      else
	{
	  for (j=0; j<outputNames_.length(); j++)
	    {
	      switch(status_[i])
		{
		case DDaceRunFailed:
		  of << "Run failed\"/>" << endl;
		  break;
		case DDacePostProcFailed:
		  of << "Post processing failed\"/>" << endl;
		  break;
		case DDaceRunPending:
		  of << "Run pending\"/>" << endl;
		  break;
		case DDaceRunNotStarted:
		  of << "Run not started\"/>" << endl;
		  break;
		default:
		  ExceptionBase::raise("DDaceServer::writeToArchive unrecognized run status");
		}
	    }
	}
      of << "</Sample>" << endl;
      
    }
  of << "</Results>" << endl;
  of << "</DDaceArchive>" << endl;
}



// ------------------------------------------------------------------------
// mutators for variable names and such
// ------------------------------------------------------------------------

void DDaceServer::setArchivePath(const String& archivePath)
{
	archivePath_ = archivePath;
}

void DDaceServer::setArchiveName(const String& archiveFilename)
{
	archiveFilename_ = archiveFilename;
}

void DDaceServer::setVariableNames(const Array<String>& varNames)
{
	if (varNames.length() != sampler_.nInputs())
		{
			ExceptionBase::raise("DDaceServer::setVariableNames: mismatched array dimension");
		}
	varNames_ = varNames;
}

void DDaceServer::setOutputNames(const Array<String>& outputNames)
{
	outputNames_ = outputNames;
}


//----------------------------------------------------------------------
// getting results
//----------------------------------------------------------------------

void DDaceServer::getResults(Array<DDaceSamplePoint>& pts, 
   		 Array<Array<double> >& funcValues) const
{
	pts = pts_;
	funcValues = results_;
}





