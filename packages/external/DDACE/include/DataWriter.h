// * Written by Leslea Lehoucq, 9/5/2000

#ifndef DATAWRITER_H
#define DATAWRITER_H

#include "SmartPtr.h"
#include "DataWriterBase.h"

/**
 * DataWriter is a reference-counted handle for DataWriterBase objects.
 */

class DataWriter 
{
 public:
  /** construct with a ptr to a DataWriterBase derived class */
  DataWriter(DataWriterBase& ptr);
  
  virtual ~DataWriter(){;};
  /** choose the input variable that will be plotted on the x axis */
  void chooseInputVariable(const String& xVarName);
  /** choose the input variables for plotting on the x and y axes */
  void chooseInputVariables(const String& xVarName, 
			    const String& yVarName);
  
  /** choose the output to be plotted on the vertical axis */
  void chooseOutputVariable(const String& outputName);
  
  /** choose any set of the input variables and output variables contained
      in the archive file or the DDace object. */
  virtual void chooseVariables(Array<String>& inputNames, 
			       Array<String>& outputNames);

  /** In some cases (most likely formatting data for further manipulation 
      instead of just plotting), want all input and output varibles in
      the DDace object written to a file in a specified format. */
  virtual void chooseAllVariables();

  /** write the plot data to a file */
  void writeToFile(const String& filename) const;
  
  /** write the plot data to a stream */
  void write(ostream& os) const ;
  
 protected:
  SmartPtr<DataWriterBase> ptr_;
};

#endif
