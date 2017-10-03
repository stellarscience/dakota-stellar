#ifndef DDACE_H
#define DDACE_H

#include "DDaceSampler.h"
#include "SmartPtr.h"
#include "DDaceMachineBase.h"
#include "DDaceReader.h"
#include <iostream>

/*

class DDace is the main user-level object in the DDace system. A DDace
object can be constructed from a DDace input file, list of DDace command
strings, or directly from components such as DDaceSampler objects
and bounds arrays. Optional methods for setting variable names, output names,
and the archive filename are available; for problems driven by input
files these names will often be provided in the input file.

Internally, DDace just contains a smart pointer to a DDaceMachineBase object
which actually does the work of generating samples and storing results. 
In uniprocessor mode, the pointer points to a DDaceUniproc object. In
multiprocessor mode, it points to either a DDaceClient (if myPid != 0) or
DDaceServer (if myPid==0) object.

 */




class DDace
{
 public:
	// ctor from input file
	DDace(const std::string& filename);

	// ctor from list of input file command strings
	DDace(const std::vector<std::string>& cmds);

	// bare-bones ctor, given sampling method and bounds.
	// Defaults are used for variable and output names and the archive filename.
	DDace(const DDaceSampler& sampler, const std::vector<double>& lower,
				const std::vector<double>& upper);

	// complete ctor, incorporating all relevant information.
	DDace(const DDaceSampler& sampler, 
		const std::vector<double>& lower,
		const std::vector<double>& upper, 
		const std::vector<std::string>& varNames,
		const std::vector<std::string>& outputNames, 
		const std::string& archiveFilename);

	
	// functions to be called in the main sampling loop:

	// getNextSample(pt) returns true if a new sample has been obtained,
	// false if no samples are left.
	bool getNextSample(DDaceSamplePoint& pt);
	// put a function evaluation into the table of results. 
	void storeFunctionValue(const DDaceSamplePoint& pt, 
				const std::vector<double>& values);
	void recordRunStatus(const DDaceSamplePoint& pt, DDaceRunStatus status);
	

	// optional variable, output, and archive file name mutators. 
	// It's OK to call these at any time. If the names have been set,
	// these will overwrite the existing values.
	void setArchivePath(const std::string& archivePath);
	void setArchiveName(const std::string& archiveName);
	void setVariableNames(const std::vector<std::string>& name);
	void setOutputNames(const std::vector<std::string>& name);

	// Output.
	void writeToArchive();	
	void writeSamples(ostream& os) const ;
	void getResults(std::vector<DDaceSamplePoint>& pts, 
			std::vector<std::vector<double> >& funcValues) const ;

	// support for IDEA: publish a list of keywords so that DDace-related
	// commands can be grabbed from an IDEA input file. 
	static const std::vector<std::string>& keywords() ;


 private:
	// Make sure pointer is not zero and throw an exception if it is. 
	void checkPtr() const ;

	// Auxiliary method to save code duplication in construction. Construction
	// from an input file generates a list of command strings, so we can
	// funnel both infile and command-list ctors through this init() method.
	void init(const std::vector<std::string>& cmds);

	// list of DDace keywords for publication to IDEA applications.
	static std::vector<std::string> keywords_;

	// The DDaceMachineBase object will actually do the work.
	SmartPtr<DDaceMachineBase> ptr_;
};

#endif




