#ifndef DDACEANALYZERBASE_H
#define DDACEANALYZERBASE_H

#include "DDace.h"
#include "String.h"
#include "StrUtils.h"
#include "DDaceXMLReader.h"
#include "DDaceSamplePoint.h"
#include "DDaceSampler.h"
#include <cstdlib>
#include <cstring>

/**
 * DDaceAnalyzerBase: 
 * Base class for analyzers.
 *
 * A "restricted analysis" in the following documentation means that for
 * the analysis, we're using only one set of input/output data.  "One
 * set" for a main effects analysis is one input variable and one 
 * output variable.  "One set" for a response surface analysis is two
 * input variables and one output variable.
 *
 * An "exhaustive analysis" in the following documentation means that for
 * the analysis, we're using all the input/output data contained in
 * the DDace object ddaceObj_.  This means that we carry out an analysis
 * for each and every analysis-appropriate pairing of the input/output
 * variables (all one input, one output combinations for main effects,
 * all two input, one output combinations for response surfaces.
 */
class DDaceAnalyzerBase
{
public:
  
  DDaceAnalyzerBase() : archiveFile_(), outputNames_(0), outputData_(0),
    ddaceObj_(), outputIndex_(-1), nSamples_(0) {;}
    

  /**
   * Construct a base object that requires the name of the archive
   * file from which to gather the data, and the name of the output
   * variable which will be used in the analysis.
   */
  DDaceAnalyzerBase(const String& archiveFile, const String& outputName);
  
  /**
   * Create a base object with only the archive file name.  This means
   * that analysis data will be generated for all feasible combinations
   * of input/output required by the current analyzer tool (i.e. main
   * effects, response surface, etc.).
   */
  DDaceAnalyzerBase(const String& archiveFile);

  /**
   * Construct a base object that requires a DDace object
   * from which to gather the data, and the name of the output
   * variable which will be used in the analysis.
   */
  DDaceAnalyzerBase(const DDace& ddaceObj, const String& outputName);
   
  /**
   * Create a base object with only a DDace object.  This means
   * that analysis data will be generated for all feasible combinations
   * of input/output required by the current analyzer tool (i.e. main
   * effects, response surface, etc.).
   */
  DDaceAnalyzerBase(const DDace& ddaceObj);

  virtual ~DDaceAnalyzerBase(){;}
  
  virtual DDaceAnalyzerBase* clone() const = 0;
  
  
  /**
   * For an analysis focused on a single output, return the name of
   * the output data vector.
   */
  virtual void getOutputName(String& outputName) const 
    { outputName = outputNames_[0];}
  
  /**
   * For an exhaustive analysis, (all possible combinations of the 
   * input/output, as appropriate to the analysis) return the names of
   * all the output variables.
   */
  virtual void getOutputNames(Array<String>& outputNames) const
    { outputNames = outputNames_;}
  
  /**
   * For an analysis with only one output variable, get the data
   * corresponding to that output varaible
   */
  virtual void getOutputData(Array<double>& outputData) const
    { outputData = outputData_;}

 protected:

  /**
   * archiveFile_ : if we generated our data from an archive file, save
   * the archive file name; otherwise, archiveFile_ will be an empty
   * string ("").
   */
  String archiveFile_;

  /**
   * outputNames_ : an array of strings of output variable names; if 
   * this is a restricted analysis, the array will be of length 1;
   * if this is an exhaustive analysis, the array will be of length
   * ddaceObj.getOutputNames(names).length()
   */

  Array<String> outputNames_;

  /**
   * outputData_ : for a restricted analysis, the outputData for the 
   * selected output variable is contained in outputData_.
   */

  Array<double> outputData_;

  /** 
   * ddaceObj_ : the DDace object from which we can get the input/output 
   * data
   */

  DDace ddaceObj_;

  /**
   * outputIndex_ : contains the index of the single outputNames_ 
   * entry if a restricted analysis; if an exhaustive analysis,
   * outputIndex_ will remain -1.
   */

  int outputIndex_;

  /// nSamples_ is the number of samples generated for each input variable
  int nSamples_;
}; 

#endif //DDACEANALYZERBASE_H


