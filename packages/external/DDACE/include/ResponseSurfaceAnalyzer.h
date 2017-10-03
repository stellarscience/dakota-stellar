#ifndef RESPONSE_SURFACE_ANALYZER_H
#define RESPONSE_SURFACE_ANALYZER_H

#include "DDace.h"
#include "Mars.h"
#include "Array.inl"
#include "Exceptions.h"
#include "DDaceXMLReader.h"
#include "TreeBuildingXMLHandler.h"
#include "XMLParser.h"
#include "FileInputSource.h"
#include "XMLUtils.h"
#include "XMLObject.h"
#include "DDaceAnalyzerBase.h"


/**
 * ResponseSurfaceAnalyzer : an analysis tool for computing the response
 * surface of a function.  The input is usually taken from the output
 * generated from a previous DDace sampling run.
 *
 */

class ResponseSurfaceAnalyzer : public DDaceAnalyzerBase
{
 public:
  
  ResponseSurfaceAnalyzer();
  
  ResponseSurfaceAnalyzer(const String& archiveName, const String& inputName1, 
			  const String& inputName2, const String& outputName);

  ResponseSurfaceAnalyzer(const DDace& ddaceObj, const String& inputName1, 
			  const String& inputName2, const String& outputName);

  ResponseSurfaceAnalyzer(const String& archiveName, const String& inputName1, 
			  const String& inputName2, const String& outputName,
			  const int nPtsPerDimension);

  ResponseSurfaceAnalyzer(const DDace& ddaceObj, const String& inputName1, 
			  const String& inputName2, const String& outputName,
			  const int nPtsPerDimension);

  virtual ~ResponseSurfaceAnalyzer(){;}

  DDaceAnalyzerBase* clone() const;

  int getNGridPts();

  void getXYGridData(Array<double>& inputData1, 
		     Array<double>& inputData2);
  
  //Array<double> getInputGridData(String& inputName);

  void getOutputGridData(Array<double>& outputGridData);

  void responseSurfaceToXML(XMLObject& rsXML) const;

  void responseSurfacePlotXML(String& plotTitle, XMLObject& rsXML) const;

 private:

  void setRSAnalyzerData(const String& inputName1, const String& inputName2);
  
 protected:
  
  Mars                  marsObj_;
  Array<String>         inputNames_;
  Array<int>            inputIndices_;

  Array<double>         fixedInputValues_;
  Array<double>         inputLowerBounds_;
  Array<double>         inputUpperBounds_;

  Array<double>         marsInput1_;
  Array<double>         marsInput2_;
  Array<double>         marsOutput_;

  int                   nGridPts_;
};

#endif //RESPONSE_SURFACE_ANALYZER_H
