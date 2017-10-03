#ifndef DDACECLIENT_H
#define DDACECLIENT_H

#include "Array.h"
#include "String.h"
#include "DDaceSamplePoint.h"
#include "DDaceMachineBase.h"

/**

DDaceClient derives from DDaceMachineBase. This class is very simple:
its only overloaded methods are a trivial virtual dtor and 
the methods to get sample values and set results.

The getNextSample() method sends the server a request for a new sample point.
When it gets a response, it returns true and returns the point by reference.
If the server has used up all points, getNextSample() is notified of that
and returns false.

storeFunctionValue() sends function evaluation to the server for storage.

 */

class DDaceClient : public DDaceMachineBase
{
 public:
	DDaceClient();
	virtual ~DDaceClient(){;}

	// overload the methods for getting sample and setting value.
	virtual bool getNextSample(DDaceSamplePoint& pt) ;
	virtual void storeFunctionValue(const DDaceSamplePoint& pt, 
																	const Array<double>& values) ;
	virtual void recordRunStatus(const DDaceSamplePoint& pt, 
															 DDaceRunStatus status);
 protected:
};

#endif
