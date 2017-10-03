#include "AsciiDataWriter.h"
#include <iomanip>
#include <fstream>
using namespace std;

#include "System.h"

DataWriterBase* AsciiDataWriter::clone() const
{
  DataWriterBase* rtn = new AsciiDataWriter(*this);
  if ( rtn == 0 ) MemoryException::raise("AsciiDataWriter::clone()");
  return rtn;
}

void AsciiDataWriter::writeToFile(const String& filename) const
{
  try{
    
    ofstream of(filename.cString());
    
    write(of);
  }
  catch(ExceptionBase& e)
    { e.trace("in AsciiDataWriter::writeToFile(String&)");}
}

void AsciiDataWriter::write(ostream& os) const
{
  try {
    os << "date = " << System::date() << endl;
    if(archiveName_.length() > 0)
      os << "Associated archive file name = " << archiveName_ << endl;
    os << "Number of input variables : " << inputNames_.length() << endl;
    os << "Number of output variables : " << outputNames_.length() << endl;

    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    int nsamps = sampler.nSamples();
    
    os << "Number of samples : " << sampler.nSamples() << endl;
    os << endl << endl;
    
    if(inputNames_.length() == 0)
      ExceptionBase::raise("AsciiDataWriter::write(ostream&) : no input variables"
			  "found.");
    Array<int> maxInputLen(inputNames_.length());
    Array<int> maxOutputLen(outputNames_.length());
    
    int i;
    for(i = 0; i < inputNames_.length(); i++)
      {
	inputNames_[i].length() > 15 ?
	  maxInputLen[i] = inputNames_[i].length() :
	  maxInputLen[i] = 15;
      }

    for(i = 0; i < outputNames_.length(); i++)
      {
	outputNames_[i].length() > 15 ?
	  maxOutputLen[i] = outputNames_[i].length() :
	  maxOutputLen[i] = 15;
      }
    
    // right justify each of the columns 
    os << setiosflags(ios::right);
    for(i = 0; i < inputNames_.length(); i++)
      os << setw(maxInputLen[i] + 1) << inputNames_[i];
    for(i = 0; i < outputNames_.length(); i++)
      os << setw(maxOutputLen[i] + 1) << outputNames_[i];
    os << endl << endl;
    
    Array<DDaceSamplePoint> pts(nsamps);
    Array<Array<double> > outputValues(nsamps);
    ddaceObj_.getResults(pts, outputValues);
    
    for(int j = 0; j < nsamps; j++)
      {
	for(int i = 0; i < inputNames_.length(); i++)
	  os << setw(maxInputLen[i] + 1) << pts[j][inputIndices_[i]];
	for(int k = 0; k < outputNames_.length(); k++)
	  os << setw(maxOutputLen[k] + 1) << outputValues[j][outputIndices_[k]];
	os << endl;
      }
  }
  catch(ExceptionBase& e)
    { e.trace("in AsciiDataWriter::write(ostream&)");}

}


