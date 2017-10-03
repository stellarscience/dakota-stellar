#ifndef MAINEFFECTSANALYZER_H
#define MAINEFFECTSANALYZER_H

#include <cmath>
#include "DDace.h"
#include "Array.inl"
#include "String.h"
#include "StrUtils.h"
#include "Exceptions.h"
#include "DDaceXMLReader.h"
#include "DDaceAnalyzerBase.h"
#include "DDaceSampler.h"
#include "DDaceLHSampler.h"
#include "TreeBuildingXMLHandler.h"
#include "XMLParser.h"
#include "FileInputSource.h"
#include "XMLUtils.h"
#include "XMLObject.h"

/**
 * MainEffectsAnalyzer : an analysis tool for studying the variance of
 * an output relative to the unique sample values of a single input 
 * variable.
 *
 */

class MainEffectsAnalyzer : public DDaceAnalyzerBase
{
 public:

  MainEffectsAnalyzer();

  /**
   * MainEffectsAnalyzer ctor for creating a MainEffectsAnalzyer object
   * from the data associated with one input variable (specified by
   * inputName), and one output variable (specified by outputName).  The
   * data is retrieved from the archive file named archiveName.
   */
  MainEffectsAnalyzer(const String& archiveName, const String& inputName, 
		      const String& outputName);

  /**
   * MainEffectsAnalyzer ctor for creating a MainEffectsAnalzyer object
   * from the data associated with a selected output (named outputName),
   * and all varied inputs in the archive file named archiveName.  The
   * data for creating the MainEffectsAnalyzer comes from the archive 
   * file named archiveName.
   */
  MainEffectsAnalyzer(const String& archiveName,const String& outputName);

  /**
   * MainEffectsAnalyzer ctor for creating a MainEffectsAnalzyer object
   * from the data associated with one input variable (specified by
   * inputName), and one output variable (specified by outputName).  The
   * data is retrieved from the ddaceObj.
   */
  MainEffectsAnalyzer(const DDace& ddaceObj, const String& inputName, 
		      const String& outputName);

  /**
   * MainEffectsAnalyzer ctor for creating a MainEffectsAnalzyer object
   * from the data associated with a selected output (named outputName),
   * and all varied inputs in the ddaceObj.  The data for creating the 
   * MainEffectsAnalyzer comes from the DDace object named ddaceObj.
   */
  MainEffectsAnalyzer(const DDace& ddaceObj, const String& outputName);

  virtual ~MainEffectsAnalyzer(){;}

  DDaceAnalyzerBase* clone() const;
  
  /**
   * If we constructed the the MainEffectsAnalyzer object with each
   * input paired with the selected output, we get back the 2-D Array 
   * of ordered, unique inputs (by column), and a 2-D Array of the 
   * corresponding output means. Each row corresponds to a unique input 
   * variable value.
   */
  void getGroupMeans(Array<Array<double> >& xOrdered,
		     Array<Array<double> >& groupMeans) const;

  /**
   * Get the Array of ordered, unique inputs, and the Array of corresponding 
   * output means for an input variable specified by inputName.  This 
   * version of getGroupMeans is to get grouped mean data for a single
   * input variable, when data for all inputs is available.
   */
  void getGroupMeans(String& inputName, Array<double>& xOrdered,
		     Array<double>& groupMeans) const;

   /**
   * If we constructed the the MainEffectsAnalyzer object with each of the
   * inputs paired with the selected output, we get back the 2-D Array 
   * of ordered, unique inputs, and the corresponding 2-D Array of 
   * standard deviations.  Each column of each 2-D Array corresponds to 
   * an input variable.
   */
 void getGroupStdDevs(Array<Array<double> >& xOrdered,
		       Array<Array<double> >& groupStdDevs) const;

  /**
   * Get the Array of ordered, unique inputs, and the Array of corresponding 
   * output standard deviations for an input variable specified by 
   * inputName.  This version of getGroupStdDevs is to get grouped standard
   * deviation data (grouped by unique input) for a single input variable,
   * when data for all inputs is available.
   */
  void getGroupStdDevs(String& inputName, 
		       Array<double>& xOrdered,
		       Array<double>& groupStdDevs) const;
  /**
   * For a main effects analysis restricted to a single specified
   * input and a single specified output, get the name of the input
   * variable and copy the name to inputName.
   */
  void getInputName(String& inputName) const;

  /**
   * For a main effects analysis for all inputs, record the names of
   * all the input variables in inputNames.
   */
  void getInputNames(Array<String>& inputNames) const 
    { inputNames = inputNames_;}

  /**
   * Get the input data associated with this main effects analysis.
   */
  void getInputData(Array<Array<double> >& inputData) const 
    {inputData = inputData_;}
  
  /**
   * Get the number of intervals into which each input variable's domain
   * was divided for the Latin Hypercube sampling.
   */
  int nIntervals() const {return nIntervals_;}

  /**
   * Reformat the main effects data to an XML format.
   */
  void mainEffectsToXML(XMLObject& meXML) const;
  
 private:
  
  //
  // Check that this is a Latin Hypercube sample; if LH, calculate 
  // the number of unique intervals into which each input is divided.
  void getSampleInformation();
  
  // If the user has provided an input name to create the Analyzer,
  // make sure it's a valid name (that it matches a variable name
  // contained in the ddaceObj_), and record it's index.  If
  // we should get stats for all inputs, then name = "ALL_INPUTS"
  // and we just copy the input variable array from the ddaceObj_
  // to inputNames_.
  void getMainEffectsInputNames(const String& name);

  // Get the input data associated with this analysis.  This might
  // be input data for only one input variable, or it may be the
  // input data for all the input variables in ddaceObj_.
  void getData();

  // For a restricted main effects analysis, calculate the mean and
  // standard deviation of the output (named outputName) for the
  // each value of the input (named inputName).
  void calculateStats();
  
  // For an exhaustive main effects analysis, calculate the mean and
  // standard deviation of the output (specified by outputName_), for 
  // the each input.
  void calculateStats(const DDace& ddaceObj);

  //  Array<Array<Array<double> > > calculateGroupMeans();

  // For an exhaustive main effects analysis, calculate the mean of the
  // output, grouped by each value of an input variable.  Do this for
  // all input/output variable combinations (thus the 2-D Array).
  Array<Array<double> > calculateGroupMeans();

  // For an exhaustive main effects analysis, calculate and return the
  // grouped output standard deviations from the output groupMeans and
  // the output data (grouped by unique input value).
  Array<Array<double> > calculateGroupStdDev(const Array<Array<double> >& 
					      groupMeans);
  
  // sortDouble takes Arrays array1, array2, and array3 and rearranges
  // the elements of all three arrays relative to the increasing
  // order of the elements in array1
  void sortDouble(Array<Array<double> >& array1, 
		  Array<Array<double> >& array2,
		  Array<Array<double> >& array3) const; 
 protected:

  // an array containing the names of all the input variables for
  // this object;  if this is an exhaustive analysis, this array
  // will contain the names of all the inputs in the ddaceObj_ (see
  // DDaceAnalyzerBase.h).  if this is a restricted analysis, this
  // array will be of length one, containing only the single input
  // name provided by the user.
  Array<String> inputNames_;

  // If this is a restricted analysis, inputData_ will contain a 1-D
  // Array of input data associated with the single name in the
  // inputNames_ Array.  If this is an exhaustive analysis, this Array
  // will contain k columns of data where k=(# if input variables in
  // the ddaceObj_). 
  Array<Array<double> > inputData_;

  // If this is a restricted analysis, then this will be a 2-D Array
  // with only one column of the unique values in inputData_, in 
  // ascending order.  If this is an exhaustive analysis, this will 
  // be an (i x j)  Array, where i indexes the unique input values in 
  // ascending order, and j indexes the input variables, gathered
  // from the ddaceObj_ (see DDaceAnalyzerBase.h).
  Array<Array<double> > orderedUniqueInput_;

  // If this is a restricted analysis, then this will be a 2-D Array
  // with only one column of the unique values in orderedGroupMeans_, 
  // ordered to correspond with the orderedUniqueInput_.  If this 
  // is an exhaustive analysis, this will be an (i x j)  Array, where i 
  // indexes the output means corresponding to values in
  // orderedUniqueInput_, and j indexes the input variables.
  Array<Array<double> > orderedGroupMeans_;

  // If this is a restricted analysis, then this will be a 2-D Array
  // with only one column of the unique values in orderedGroupStdDevs_, 
  // ordered to correspond with the orderedUniqueInput_.  If this 
  // is an exhaustive analysis, this will be an (i x j)  Array, where i 
  // indexes the output standard deviations corresponding to values in
  // orderedUniqueInput_, and j indexes the input variables.
  Array<Array<double> > orderedGroupStdDevs_;

  // nIntervals_ is the number of intervals into which the input data
  // was divided (this should be the same number for each input variable 
  // in a Latin Hypercube sampling scheme).
  int nIntervals_;

  // If this is a restricted analysis, this will be the index of the 
  // input variable used in this analysis.  Otherwise inputIndex_ will
  // be set to -1 to signify that this is an exhaustive analysis.
  int inputIndex_;
};

#endif //MAINEFFECTSANALYZER_H


