#include "DDaceUniproc.h"


//-----------------------------------------------------------------------
// Ctors for DDaceUniproc. These do nothing beyond invoking the ctors 
// for DDaceServer.
//-----------------------------------------------------------------------

DDaceUniproc::DDaceUniproc(const DDaceSampler& sampler, 
			 const std::vector<std::string>& varNames,
			 const std::vector<string>& outputNames, 
			 const string& archiveFilename)
	: DDaceServer(sampler, varNames, outputNames, archiveFilename)
{
	;
}

DDaceUniproc::DDaceUniproc(const DDaceSampler& sampler)
	: DDaceServer(sampler)
{
	;
}

DDaceUniproc::DDaceUniproc(const XMLObject& xmlObj)
	: DDaceServer(xmlObj)
{
	;
}

//-----------------------------------------------------------------------
// In uniprocessor mode, getNextSample() just grabs the next point out
// of the list of sample points and returns true. 
// When all points have been used, return
// false to terminate the sampling loop.
//-----------------------------------------------------------------------

bool DDaceUniproc::getNextSample(DDaceSamplePoint& pt)
{
	if (stackPtr_ < pts_.length())
		{
			pt = pts_[stackPtr_];
			stackPtr_++;
			writeToArchive();
			return true;
		}
	else
		{
			return false;
		}
}

//-----------------------------------------------------------------------
// In uniprocessor mode, storeFunctionValue() just puts the returned function
// value into the results array. The index() method of DDaceSamplePoint
// is used to determine its place in the results array.
//-----------------------------------------------------------------------

void DDaceUniproc::storeFunctionValue(const DDaceSamplePoint& pt,
				 const vector<double>& values)
{
	results_[pt.index()] = values;
	status_[pt.index()] = DDaceRunOK;
}

void DDaceUniproc::recordRunStatus(const DDaceSamplePoint& pt,
					 DDaceRunStatus status)
{
	status_[pt.index()] = status;
}
