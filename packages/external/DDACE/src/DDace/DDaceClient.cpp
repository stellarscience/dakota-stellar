#include "DDaceClient.h"
#include "ArrayComm.h"
#include "DDace.h"

DDaceClient::DDaceClient()
	: DDaceMachineBase()
{
	pts_.bcast(0, DDace::getComm());
}

bool DDaceClient::getNextSample(DDaceSamplePoint& pt)
{

	if (pts_.length()==0) return false;

	// send a request for the index to the next point
	int request = 1;
	cerr << "Client sending request for sample point" << endl;
	PMachine::send(&request, 1, PMachine::INT, (int) RequestTag, 0,
								 DDace::getComm());

	// wait for an answer; use blocking recv
	int index;
	cerr << "Client waiting for sample point" << endl;
	PMachine::recv(&index, 1, PMachine::INT, (int) IndexTag, 0,
								 DDace::getComm());

	// if index < 0, we've used up all point. Return false to terminate
	// sampling loop.
	cerr << "Client received sample point" << endl;
	if (index < 0)
		{
			return false;
		}
	
	// else look up the sample point at the given index, and return it.
	pt = pts_[index];

	return true;
}

void DDaceClient::storeFunctionValue(const DDaceSamplePoint& pt,
																		 const Array<double>& values)
{
	// send results back to host for storage
	int index = pt.index();
	//	cerr << "sending " << index << " " << values << endl;
	PMachine::send(&index, 1, PMachine::INT, (int) ResultIndexTag, 0,
								 DDace::getComm());

	// send status=OK message

	int status = (int) DDaceRunOK;
	PMachine::send(&status, 1, PMachine::INT, (int) ResultIndexTag, 0,
								 DDace::getComm());

	// send result array
	ArrayComm::send(values, (int) ResultValueTag, 0, DDace::getComm());
}

void DDaceClient::recordRunStatus(const DDaceSamplePoint& pt, 
																	DDaceRunStatus status)
{
	// send results back to host for storage
	int index = pt.index();
	//	cerr << "sending " << index << " " << values << endl;
	PMachine::send(&index, 1, PMachine::INT, (int) ResultIndexTag, 0,
								 DDace::getComm());

	// send status=OK message

	int tmp = (int) status; 
	PMachine::send(&tmp, 1, PMachine::INT, (int) ResultIndexTag, 0,
								 DDace::getComm());
}






