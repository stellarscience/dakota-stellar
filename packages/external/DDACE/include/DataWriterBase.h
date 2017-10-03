#ifndef DATAWRITERBASE_H
#define DATAWRITERBASE_H

#include "String.h"
#include "DDace.h"
#include "DDaceXMLReader.h"

/**
 * DataWriterBase subclasses translate the internal DDace archive format into
 * something that can be understood by a third-party visualization
 * package (VTK, Grace, Matlab), or by a third-party analysis tool (S-Plus, 
 * Matlab, etc.).  
 * This base class contains all the accessor methods:
 * 1) single input/single output variable (2-D plotting)
 * 2) two inputs/one output (3-D plotting)
 * 3) an arbitrary set of inputs and outputs (manipulation of data 
 *    by another package)
 * 4) selection of all inputs and outputs recorded in the archive file
 *    or DDace object.
 */

class DataWriterBase 
{
 public:
  /** construct from the DDace object containing the archive */
  DataWriterBase(const DDace& ddaceObject)
    : ddaceObj_(ddaceObject), inputNames_(0), inputIndices_(0), 
    outputNames_(0), outputIndices_(0), archiveName_("") {;}

  /** construct from a DDace archive file */
  DataWriterBase(const String& archiveFilename);
  
  virtual ~DataWriterBase(){;}
  
  /** choose the input variable that will be plotted on the x axis */
  virtual void chooseInputVariable(const String& xVarName);
  /** choose the input variables for plotting on the x and y axes */
  virtual void chooseInputVariables(const String& xVarName, 
				    const String& yVarName);
  /** choose all input variables (usually for writing to an ASCII 
      formatted file that will be passed on to another application). */
  virtual void chooseAllInputVariables();

  /** choose the output to be plotted on the vertical axis */
  virtual void chooseOutputVariable(const String& outputName);
  
  /** choose any set of the input variables and output variables contained
      in the archive file or the DDace object. */
  virtual void chooseVariables(Array<String>& inputNames, 
			       Array<String>& outputNames);

  /** In some cases (most likely formatting data for further manipulation 
      instead of just plotting), want all input and output varibles in
      the DDace object written to a file in a specified format. */
  virtual void chooseAllVariables();

  virtual DataWriterBase* clone() const = 0;

  /** write the plot data to a file */
  virtual void writeToFile(const String& filename) const = 0 ;
  
  /** write the plot data to a stream */
  virtual void write(ostream& os) const = 0 ;
  
 protected:
  DDace ddaceObj_;
  Array<String> inputNames_;
  Array<int> inputIndices_;
  Array<String> outputNames_;
  Array<int> outputIndices_;
  String archiveName_;
};

#endif
