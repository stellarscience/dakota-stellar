#ifndef DDACE_H
#define DDACE_H

#include "Array.h"
#include "String.h"
#include "DDaceSampler.h"
#include "SmartPtr.h"
#include "DDaceMachineBase.h"
#include "PComm.h"
#include <iostream>
#include "XMLObject.h"

/**

class DDace is the main user-level object in the DDace system. A DDace
object can be constructed from a DDace input file, list of DDace command
strings, an XMLObject, or directly from components such as DDaceSampler 
objects and bounds arrays. Optional methods for setting variable names,
output names, and the archive filename are available; for problems driven 
by input files these names will often be provided in the input file.

Internally, DDace just contains a smart pointer to a DDaceMachineBase object
which actually does the work of generating samples and storing results. 
In uniprocessor mode, the pointer points to a DDaceUniproc object. In
multiprocessor mode, it points to either a DDaceClient (if myPid != 0) or
DDaceServer (if myPid==0) object.

 */


class DDace
{
 public:

  // Empty ctor.  This constructor is only called to throw an execption
  // if a "real" DDace object cannot be constructed.

  // L 7-6-2000 : don't throw exception for empty ctor.  Need the empty
  // ctor to be able to create a DDace object as a member of another
  // class.
  DDace();
  
	// bare-bones ctor, given sampling method and bounds.
	// Defaults are used for variable and output names and the archive filename.
	DDace(const DDaceSampler& sampler);

	// complete ctor, incorporating all relevant information.
	DDace(const DDaceSampler& sampler, const Array<String>& varNames,
				const Array<String>& outputNames, const String& archiveFilename);

	// XMLObject ctor -- translate info in an XMLObject into a DDace object.
	DDace(const XMLObject& xmlObj);

	// getNextSample(pt) returns true if a new sample has been obtained,
	// false if no samples are left.
	bool getNextSample(DDaceSamplePoint& pt);

	// put a function evaluation into the table of results. 
	void storeFunctionValue(const DDaceSamplePoint& pt, 
													const Array<double>& values);
	void recordRunStatus(const DDaceSamplePoint& pt, DDaceRunStatus status);
	

	// optional variable, output, and archive file name mutators. 
	// It's OK to call these at any time. If the names have been set,
	// these will overwrite the existing values.
	void setArchivePath(const String& archivePath);
	void setArchiveName(const String& archiveName);
	void setVariableNames(const Array<String>& name);
	void setOutputNames(const Array<String>& name);
	
	// Output.
	void getSampler(DDaceSampler& sampler) const;
	void getRunStatus(Array<DDaceRunStatus>& status) const ;
	void getArchiveFilename(String& archiveFilename) const ;
	void getVarNames(Array<String>& names) const;
	void getOutputNames(Array<String>& names) const;
	int  getNumOutputs() const;

	void writeToArchive();	
	void writeSamples(ostream& os) const ;
	void getResults(Array<DDaceSamplePoint>& pts, 
			Array<Array<double> >& funcValues) const ;
	
	// support for IDEA: publish a list of keywords so that DDace-related
	// commands can be grabbed from an IDEA input file. 
	static const Array<String>& keywords() ;
	
	// set parallel communicator
	static void setComm(const PComm& comm) {comm_ = comm; commReady_=true;}
	static const PComm& getComm() {return myComm();}


 private:
	// Make sure pointer is not zero and throw an exception if it is. 
	void checkPtr() const ;

	// list of DDace keywords for publication to IDEA applications.
	static Array<String> keywords_;

	// The DDaceMachineBase object will actually do the work.
	SmartPtr<DDaceMachineBase> ptr_;

	static const PComm& myComm();
	static PComm comm_;
	static bool commReady_;
};

#endif




