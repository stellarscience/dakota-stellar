#include "ResponseSurfaceAnalyzer.h"

ResponseSurfaceAnalyzer::ResponseSurfaceAnalyzer()
  : DDaceAnalyzerBase(), marsObj_(), inputNames_(2), 
    inputIndices_(2,-1), fixedInputValues_(0), inputLowerBounds_(0), 
    inputUpperBounds_(0), marsInput1_(0), marsInput2_(0), 
    marsOutput_(0)
{;}  

ResponseSurfaceAnalyzer::ResponseSurfaceAnalyzer(const String& archiveName, 
						 const String& inputName1, 
						 const String& inputName2, 
						 const String& outputName)
  : DDaceAnalyzerBase(archiveName, outputName), inputNames_(2),
    inputIndices_(2, -1), fixedInputValues_(0), inputLowerBounds_(0), 
    inputUpperBounds_(0), marsInput1_(0), marsInput2_(0), 
    marsOutput_(0)
{
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    marsObj_ = Mars(sampler.nInputs(), sampler.nSamples());

  // By default, put 20 grid points per dimension (i.e. 400 grid
  // points for a 2-D grid).  If you want control over number
  // of grid points, use the ResponseSurfaceAnalyzer ctor() with the
  // nPtsPerDim input.
  marsObj_.setNPtsPerDim(20);
  nGridPts_ = 20*20;
  setRSAnalyzerData(inputName1, inputName2);

  // check to see how all the variables are set (this is error checking...):
  
  for(int k = 0; k < sampler.nInputs() - 1; k++)
    cerr << fixedInputValues_[k] << ", ";

  cerr << fixedInputValues_[sampler.nInputs() -1] << endl;
  cerr << "Input lower bounds: " << inputLowerBounds_ << endl;
  cerr << "Input upper bounds: " << inputUpperBounds_ << endl;
}

ResponseSurfaceAnalyzer::ResponseSurfaceAnalyzer(const DDace& ddaceObj,
						 const String& inputName1,
						 const String& inputName2,
						 const String& outputName)
  : DDaceAnalyzerBase(ddaceObj, outputName), inputNames_(2), 
  inputIndices_(2,-1), fixedInputValues_(0), inputLowerBounds_(0), inputUpperBounds_(0),
  marsInput1_(0), marsInput2_(0), marsOutput_(0)
{
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    marsObj_ = Mars(sampler.nInputs(), sampler.nSamples());
    
  // By default, put 20 grid points per dimension (i.e. 400 grid
  // points for a 2-D grid).  If you want control over number
  // of grid points, use the ResponseSurfaceAnalyzer with the
  // nPtsPerDim input.
  marsObj_.setNPtsPerDim(20);
  nGridPts_ = 20*20;
  setRSAnalyzerData(inputName1, inputName2);
}

ResponseSurfaceAnalyzer::ResponseSurfaceAnalyzer(const String& archiveName, 
						 const String& inputName1, 
						 const String& inputName2, 
						 const String& outputName,
						 const int nPtsPerDimension)
  : DDaceAnalyzerBase(archiveName, outputName), marsObj_(), inputNames_(2), 
  inputIndices_(2,-1), fixedInputValues_(0),inputLowerBounds_(0), inputUpperBounds_(0),
  marsInput1_(0), marsInput2_(0), marsOutput_(0)
{
  // Assign nPtsPerDimension grid points per dimension (i.e. 
  // (nPtsPerDimension * nPtsPerDimension) grid points for a 
  // 2-D grid).  If you want control over number of grid 
  // points, use the ResponseSurfaceAnalyzer with the
  // nPtsPerDim input.
  
  DDaceSampler sampler;
  ddaceObj_.getSampler(sampler);
  marsObj_ = Mars(sampler.nInputs(), sampler.nSamples());
  
  int tmpPtsPerDim = nPtsPerDimension;
  if(nPtsPerDimension > 30)
    {
      cerr << "You have selected a  number of points per dimension \nthat "
	   << "exceeds the current maximum allowable, 30.  Your number \n"
	   << "of points per dimension is reset to 30.\n";
      tmpPtsPerDim = 30;
    }
  marsObj_.setNPtsPerDim(tmpPtsPerDim);
  nGridPts_ = tmpPtsPerDim * tmpPtsPerDim;
  setRSAnalyzerData(inputName1, inputName2);
}

ResponseSurfaceAnalyzer::ResponseSurfaceAnalyzer(const DDace& ddaceObj,
						 const String& inputName1,
						 const String& inputName2,
						 const String& outputName,
						 const int nPtsPerDimension)
  : DDaceAnalyzerBase(ddaceObj, outputName), marsObj_(), inputNames_(2), 
  inputIndices_(2,-1), fixedInputValues_(0), inputLowerBounds_(0), inputUpperBounds_(0),
  marsInput1_(0), marsInput2_(0),  marsOutput_(0)
{
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    marsObj_ = Mars(sampler.nInputs(), sampler.nSamples());

    int tmpPtsPerDim = nPtsPerDimension;
  // Assign nPtsPerDimension grid points per dimension (i.e. 
  // (nPtsPerDimension * nPtsPerDimension) grid points for a 
  // 2-D grid).  If you want control over number of grid 
  // points, use the ResponseSurfaceAnalyzer with the
  // nPtsPerDim input.
  if(nPtsPerDimension > 30)
    {
      cerr << "You have selected a  number of points per dimension \nthat "
	   << "exceeds the current maximum allowable, 30.  Your number \n"
	   << "of points per dimension is reset to 30.\n";
      tmpPtsPerDim = 30;
    }
  marsObj_.setNPtsPerDim(tmpPtsPerDim);
  nGridPts_ = tmpPtsPerDim * tmpPtsPerDim;
  setRSAnalyzerData(inputName1, inputName2);
}

void ResponseSurfaceAnalyzer::setRSAnalyzerData(const String& inputName1, 
					   const String& inputName2)
{
  try{
    //create the Mars objects marsObj_.
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);

    // If the user has provided legal input names (that is, the names
    // provided correspond to varNames in the ddaceObj_), assign
    // these names to inputNames_, and assign their indices to 
    // inputIndices_.


    if(inputNames_.length() == 0) inputNames_.resize(2);
    Array<String> varNames(sampler.nInputs());
    ddaceObj_.getVarNames(varNames);

    for(int i = 0; i < varNames.length(); i++)
      {
	if((inputIndices_[0] == -1) || (inputIndices_[1] == -1))
	  {
	    if(varNames[i] == inputName1)
	      {
		inputNames_[0] = inputName1;
		inputIndices_[0] = i;
		continue;
	      }
	    else if(varNames[i] == inputName2)
	      {
		inputNames_[1] = inputName2;
		inputIndices_[1] = i;
		continue;
	      }
	    else continue;
	  }
	else break;
      }
    if((inputIndices_[0] == -1) || (inputIndices_[1] == -1))
      {
	ExceptionBase::raise("ResponseSurfaceAnalyzer ctor() : at least one\n"
			    "of your inputVariables is erroneous.  Please\n"
			    "check that you have seleced valid input names,\n"
			    "and try again.");
	//exit;
      }

    // get the bounds on the input variables
    inputLowerBounds_.resize(sampler.nInputs());
    inputUpperBounds_.resize(sampler.nInputs());
    
    inputLowerBounds_ = sampler.lowerBounds();
    inputUpperBounds_ = sampler.upperBounds();
    
    marsObj_.setBounds(inputLowerBounds_, inputUpperBounds_);

    fixedInputValues_.resize(sampler.nInputs());

    // set the fixedInputValues_ to the mean of the distribution for
    // each input variable:
    for(int j=0; j < sampler.nInputs(); j++)
      {
	fixedInputValues_[j] = sampler.dist()[j].mean();
      }
    
    Array<DDaceSamplePoint> inputs;
    Array<Array<double> > results;
    ddaceObj_.getResults(inputs, results);

    Array<Array<double> > marsInputGrid;
    cerr << "nGridPts (for setting rows of marsInputGrid): " << nGridPts_ 
	 << endl;

    marsInputGrid.resize(nGridPts_);

    for(int k=0; k < nGridPts_; k++)
      marsInputGrid[k].resize(sampler.nInputs());
    marsOutput_.resize(nGridPts_);
    marsInput1_.resize(nGridPts_);
    marsInput2_.resize(nGridPts_);

    marsObj_.generate2DGridData(inputs, outputData_, 
			    inputIndices_[0], inputIndices_[1], 
			    fixedInputValues_, marsInputGrid,
			    marsOutput_);
    cerr << "nGridPts: " << nGridPts_ << endl;
    cerr << "marsInput length: " << marsInput1_.length() << endl;
    cerr << "marsInputGrid rows: " << marsInputGrid.length() << endl;

    for(int l=0; l < nGridPts_; l++)
      {
	marsInput1_[l] = marsInputGrid[l][inputIndices_[0]];
	marsInput2_[l] = marsInputGrid[l][inputIndices_[1]];
      }
    cerr << "Mars data:" << marsInputGrid << endl;
    cerr << "Mars output:" << marsOutput_ << endl;
  }
  catch(ExceptionBase& e)
    { e.trace("in ResponseSurfaceAnalyzer ctor()"); }
}

int ResponseSurfaceAnalyzer::getNGridPts()
{
  return nGridPts_;
}

void ResponseSurfaceAnalyzer::getXYGridData(Array<double>& inputData1, 
					   Array<double>& inputData2)
{
  inputData1 = marsInput1_;
  inputData2 = marsInput2_;
}

void ResponseSurfaceAnalyzer::getOutputGridData(Array<double>& outputGridData)
{
  outputGridData = marsOutput_;
}


DDaceAnalyzerBase* ResponseSurfaceAnalyzer::clone() const
{
  DDaceAnalyzerBase* rtn = new ResponseSurfaceAnalyzer(*this);
  if (rtn == 0) MemoryException::raise("ResponseSurfaceAnalyzer::clone()");
  return rtn;
}


void ResponseSurfaceAnalyzer::responseSurfaceToXML(XMLObject& rsXML) const
{
  try
    {
      XMLUtils::init();
      
      XMLObject rs("RespnseSurfaceAnalyzer");
      
      XMLObject archive("Archive");
      archive.addAttribute("name", archiveFile_);
      rs.addChild(archive);
      
      XMLObject inputName1("Input1");
      inputName1.addAttribute("name", inputNames_[0]);
      rs.addChild(inputName1);
      
      XMLObject inputName2("Input2");
      inputName2.addAttribute("name", inputNames_[1]);
      rs.addChild(inputName2);
      
      XMLObject outputName("Output");
      outputName.addAttribute("name", outputNames_[0]);
      rs.addChild(outputName);
      
      XMLObject rsData("ResponseSurfaceData");
      rs.addChild(rsData);
      
      for(int i=0; i < nGridPts_; i++)
	{
	  XMLObject dataPoint("Data");      
	  dataPoint.addAttribute("input1", Double(marsInput1_[i]).toString());
	  dataPoint.addAttribute("input2", Double(marsInput2_[i]).toString());
	  dataPoint.addAttribute("output", Double(marsOutput_[i]).toString());
	  
	  rsData.addChild(dataPoint);
	}
      rsXML = rs;
    }
  catch(ExceptionBase& e)
    {
      e.trace("in ResponseSurfaceAnalyzer::responseSurfaceToXML()");
    }
}

void ResponseSurfaceAnalyzer::responseSurfacePlotXML(String& plotTitle,
						     XMLObject& rsXML) const
{
  try
    {
      int i, j, k;
      Array<double> xGrid;
      Array<double> yGrid;
      Array<Array<double> > zGrid;

      const int nPtsPerDim = marsObj_.getNPtsPerDim();

      xGrid.resize(nPtsPerDim);
      yGrid.resize(nPtsPerDim);
      zGrid.resize(nPtsPerDim);

      for(i=0; i < nPtsPerDim; i++)
	{
	  xGrid[i] = marsInput1_[i*nPtsPerDim];
	  cerr << "xGrid[" << i << "] = " << xGrid[i] << endl;
	  yGrid[i] = marsInput2_[i];
	  cerr << "yGrid[" << i << "] = " << yGrid[i] << endl;
	}

      //int l=0;
      cerr << "zGrid: " << endl;
      for(j=0; j < nPtsPerDim; j++){
	zGrid[j].resize(nPtsPerDim);
	cerr << "Row [" << j << "] : ";
	for(k=0; k < nPtsPerDim; k++){
	  //zGrid[j][k] = marsOutput_[l];
	  zGrid[j][k] = marsOutput_[j*nPtsPerDim + k];
	  cerr << zGrid[j][k] << " ";
	  //l++;
	}
	cerr << endl;
      }


      XMLUtils::init();
      
      XMLObject pointSet("PointSet");
      pointSet.addAttribute("xLabel", inputNames_[0]);
      pointSet.addAttribute("yLabel", inputNames_[1]);
      pointSet.addAttribute("zLabel", outputNames_[0]);
      pointSet.addAttribute("title", plotTitle);

      XMLObject xVector("XGridVector");
      XMLObject yVector("YGridVector");
      for(i=0; i < nPtsPerDim; i++){
	XMLObject x("X");
	x.addAttribute("value", Double(xGrid[i]).toString());
	xVector.addChild(x);
	
	XMLObject y("Y");
	y.addAttribute("value", Double(yGrid[i]).toString());
	yVector.addChild(y);
      }
      pointSet.addChild(xVector);
      pointSet.addChild(yVector);
      
      XMLObject zGridXML("ZGrid");
      for(i=0; i < nPtsPerDim; i++){
	XMLObject zRow("Row");      
	for(j=0; j < nPtsPerDim; j++)
	{
	  XMLObject z("Z");
	  z.addAttribute("value", Double(zGrid[i][j]).toString());
	  zRow.addChild(z);
	}
	zGridXML.addChild(zRow);
      }
      pointSet.addChild(zGridXML);
      
      rsXML = pointSet;
    }
  catch(ExceptionBase& e)
    {
      e.trace("in ResponseSurfaceAnalyzer::responseSurfacePlotXML()");
    }
}








