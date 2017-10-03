// * Written by Leslea Lehoucq, 9/5/2000

#include "DataWriter.h"

/**
 * DataWriter is a reference-counted handle for DataWriterBase objects.
 */


  /** construct with a ptr to a DataWriterBase derived class */
  DataWriter::DataWriter(DataWriterBase& ptr)
    : ptr_(ptr.clone())
  {;}
  
  /** choose the input variable that will be plotted on the x axis */
  void DataWriter::chooseInputVariable(const String& xVarName)
  {
    ptr_->chooseInputVariable(xVarName);
  }

  /** choose the input variables for plotting on the x and y axes */
  void DataWriter::chooseInputVariables(const String& xVarName, 
			    const String& yVarName)
  {
    ptr_->chooseInputVariables(xVarName, yVarName);
  }
  
  /** choose the output to be plotted on the vertical axis */
  void DataWriter::chooseOutputVariable(const String& outputName)
  {
    ptr_->chooseOutputVariable(outputName);
  }
  
  /** choose any set of the input variables and output variables contained
      in the archive file or the DDace object. */
  void DataWriter::chooseVariables(Array<String>& inputNames, 
			       Array<String>& outputNames)
  {
    ptr_->chooseVariables(inputNames, outputNames);
  }

  /** In some cases (most likely formatting data for further manipulation 
      instead of just plotting), want all input and output varibles in
      the DDace object written to a file in a specified format. */
  void DataWriter::chooseAllVariables()
  {
    ptr_->chooseAllVariables();
  }

  /** write the plot data to a file */
  void DataWriter::writeToFile(const String& filename) const
  {
    ptr_->writeToFile(filename);
  }
  
  
  /** write the plot data to a stream */
  void DataWriter::write(ostream& os) const
  {
    ptr_->write(os);
  }
