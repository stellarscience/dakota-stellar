#include "DDaceReaderBase.h"
#include "StrUtils.h"
#include "DDaceLHSampler.h"
#include "DDaceOASampler.h"
#include "DDaceRandomSampler.h"
#include "Array.cpp"
#include <iostream>
#include <fstream>
using namespace std;

DDaceReaderBase::DDaceReaderBase(const String& filename)
	: nInputs_(0), cmds_(0), tokens_(0)
{
	ifstream in(filename.cString());
	if (!in) fatalError("DDaceReaderBase::DDaceReaderBase failed to open file");
	
	cmds_ = readFile(in, '#');
	tokens_.resize(cmds_.length());
	
	for (int i=0; i<cmds_.length(); i++)
		{
			tokens_[i] = stringTokenizer(cmds_[i]);
		}
	nInputs_ = countInputs();
}
	
DDaceReaderBase::DDaceReaderBase(const Array<String>& cmds)
	: nInputs_(0), cmds_(cmds), tokens_(cmds.length())
{
	for (int i=0; i<cmds.length(); i++)
		{
			tokens_[i] = stringTokenizer(cmds[i]);
		}
	nInputs_ = countInputs();
}

bool DDaceReaderBase::lookup(const String& keyword, String& cmd) const 
{
  // dumb linear search for match for keyword
  
  for (int i=0; i<tokens_.length(); i++)
    {
      if (tokens_[i].length() == 0) continue;
      if (tokens_[i][0].allCaps() == keyword.allCaps())
	{
	  if (tokens_[i].length()==1) cmd="";
	  else
	    {
	      Array<String> remainder(tokens_[i].length()-1);
	      for (int j=1; j<tokens_[i].length(); j++)
		{
		  remainder[j-1] = tokens_[i][j];
		}
	      cmd = reassembleFromTokens(remainder);
	    }
	  return true;
	}
    }
  cmd = "";
  return false;
}

void DDaceReaderBase::getBounds(Array<double>& lower, Array<double>& upper) const 
{
  lower.reserve(10);
  upper.reserve(10);
  
  for (int i=0; i<tokens_.length(); i++)
    {
      if (tokens_[i].length() < 2) continue;
      if (tokens_[i][0].allCaps()=="BOUNDS")
	{
	  if (tokens_[i].length() != 4) 
	    fatalError("DDaceReaderBase::getBounds: invalid bounds specification");
	  lower.append(atof(tokens_[i][2].cString()));
	  upper.append(atof(tokens_[i][3].cString()));
	}
    }
  cerr << "lower bounds: " << lower << endl;
  cerr << "upper bounds: " << upper << endl;
}

void DDaceReaderBase::getVarNames(Array<String>& varNames) const 
{
  varNames.reserve(10);
  
  for (int i=0; i<tokens_.length(); i++)
    {
      if (tokens_[i].length() < 2) continue;
      if (tokens_[i][0].allCaps()=="BOUNDS")
	{
	  varNames.append(tokens_[i][1]);
	}
    }
}


void DDaceReaderBase::getArchiveFilename(String& archiveFilename) const 
{
  fatalError("DDaceReaderBase::getArchiveFilename shouldn't be called");
}

void DDaceReaderBase::getArchivedData(Array<DDaceSamplePoint>& /* pts */,
																			Array<Array<double> >& /*results*/,
																			Array<DDaceRunStatus>& /*status*/) const
{
  fatalError("DDaceReaderBase::getArchivedData shouldn't be called");
}

int DDaceReaderBase::countInputs() const 
{
  int count = 0;
  for (int i=0; i<tokens_.length(); i++)
    {
      if (tokens_[i].length() < 2) continue;
      if (tokens_[i][0].allCaps()=="BOUNDS") count++;
    }
  return count;
}



void DDaceReaderBase::getSampler(DDaceSampler& sampler) const 
{
  String method;
  String nSampleString;
  String perturbString;
  String repString;
  
  if (!lookup("METHOD", method)) 
    fatalError("DDaceReaderBase::getSampler: method spec not found");
  
  if (!lookup("SAMPLES", nSampleString))
    fatalError("DDaceReaderBase::getSampler: num sample spec not found");
  int nSamples = atoi(nSampleString.cString());
  
  bool perturb = lookup("PERTURB", perturbString);
  
  if (method.allCaps()=="LH")
    {
      if (perturb && perturbString.allCaps()=="TRUE")
	{
	  perturb = true;
	}
      else
	{
	  perturb = false;
	}
      if (!lookup("REPLICATIONS", repString))
	fatalError("DDaceReaderBase::getSampler: LH replications spec not found");
      int nReps = atoi(repString.cString());
      sampler = DDaceLHSampler(nSamples, nInputs_, nReps, perturb);
      return;
    }
  
	if (method.allCaps()=="RANDOM")
	  {
	    sampler = DDaceRandomSampler(nSamples, nInputs_);
	    return;
	  }
	
	if (method.allCaps()=="OA")
		{
			if (perturb && perturbString.allCaps()=="TRUE")
				{
					perturb = true;
				}
			else
				{
					perturb = false;
				}
			sampler = DDaceOASampler(nSamples, nInputs_, perturb);
			return;
		}

	fatalError("DDaceReaderBase::getSampler: no match to method found");
}








