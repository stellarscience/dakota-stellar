// This class is based on the MainEffectAnalyzer created by 
// Charles Tong (December 1998).  Current version by Leslea Lehoucq
// and Kevin Long.

#include "MainEffectsAnalyzer.h"

#define abs(x) (((x) > 0.0) ? (x) : -(x))


MainEffectsAnalyzer::MainEffectsAnalyzer()
  :DDaceAnalyzerBase(), inputNames_(0), inputData_(0), 
  orderedUniqueInput_(0), orderedGroupMeans_(0), 
  orderedGroupStdDevs_(0), nIntervals_(0), inputIndex_(-1)
{
  ExceptionBase::raise("MainEffectsAnalyzer() empty ctor : call to the\n"
		       "empty construtor suggests error.");
}

MainEffectsAnalyzer::MainEffectsAnalyzer(const String& archive,
					 const String& inputName,
					 const String& outputName)
  : DDaceAnalyzerBase(archive, outputName), inputNames_(1), 
  inputData_(0), orderedUniqueInput_(0), orderedGroupMeans_(0), 
  orderedGroupStdDevs_(0), nIntervals_(0), inputIndex_(-1)
{

  // get number of samples, number of replications information.
  getSampleInformation();

  // verify that inputName is a feasible input name in the ddaceObj_
  // (generated from the archive), and if so, get the index of this 
  // inputName in the varibleNames Array in ddaceObj_.
  getMainEffectsInputNames(inputName);

  // Get the input data corresponding to the user selected input variable,
  // and store it in  (nSamples_) x 1 Array called inputData_.
  getData();
  
  calculateStats();
}
  
MainEffectsAnalyzer::MainEffectsAnalyzer(const String& archive,
					 const String& outputName)
  : DDaceAnalyzerBase(archive, outputName), inputNames_(0), 
  inputData_(0), orderedUniqueInput_(0), orderedGroupMeans_(0), 
  orderedGroupStdDevs_(0), nIntervals_(0), inputIndex_(-1)
{
  getSampleInformation();
  getMainEffectsInputNames("ALL_INPUTS");
  getData();
  calculateStats();
}

MainEffectsAnalyzer::MainEffectsAnalyzer(const DDace& ddaceObj,
					 const String& inputName,
					 const String& outputName)
  : DDaceAnalyzerBase(ddaceObj, outputName), inputNames_(1), 
  inputData_(0), orderedUniqueInput_(0), orderedGroupMeans_(0), 
  orderedGroupStdDevs_(0), nIntervals_(0), inputIndex_(-1)
{

  // get number of samples, number of replications information.
  getSampleInformation();

  // verify that inputName is a feasible input name in the ddaceObj_
  // (generated from the archive), and if so, get the index of this 
  // inputName in the varibleNames Array in ddaceObj_.
  getMainEffectsInputNames(inputName);

  // Get the input data corresponding to the user selected input variable,
  // and store it in  (nSamples_) x 1 Array called inputData_.
  getData();
  
  calculateStats();
}
  
MainEffectsAnalyzer::MainEffectsAnalyzer(const DDace& ddaceObj,
					 const String& outputName)
  : DDaceAnalyzerBase(ddaceObj, outputName), inputNames_(0), 
  inputData_(0), orderedUniqueInput_(0), orderedGroupMeans_(0), 
  orderedGroupStdDevs_(0), nIntervals_(0),inputIndex_(-1)
{
  getSampleInformation();
  getMainEffectsInputNames("ALL_INPUTS");
  getData();
  calculateStats();
}

DDaceAnalyzerBase* MainEffectsAnalyzer::clone() const
{
  DDaceAnalyzerBase* rtn = new MainEffectsAnalyzer(*this);
  if (rtn == 0) MemoryException::raise("MainEffectsAnalyzer::clone()");
  return rtn;
}

// *******************************************************************
// Verify Latin Hypercube sampling; get the number of unique input 
// values per input variable, nIntervals_.
// -------------------------------------------------------------------
void MainEffectsAnalyzer::getSampleInformation() 
{
  try{
    // First, make sure that the sampler from the ddaceObj is Latin
    // Hypercube -- this analysis method requires samples to be
    // from a Latin Hypercube.
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    if (sampler.typeName() != "DDaceLHSampler")
      {
	ExceptionBase::raise("Incorrect sampler associated with samples in "
			     "MainEffectsAnalyzer ctor.\n"
			     "Samples must be generated from a Latin "
			     "Hypercube scheme.");
      }
    
    // Next, make sure that then number of samples is a multiple of
    // the number of replications (this should have been checked when
    // the sampler member of the ddaceObj was created, but we'll
    // double check to be safe).
    
    int nreplications = sampler.getParameter("replications");
    int nintervals = nSamples_/nreplications;
    
    if ( nreplications*nintervals != nSamples_ )
      {
	ExceptionBase::raise("MainEffectsAnalyzer ctor : "
			     "nSamples_ should \nbe a multiple of "
			     "nreplications.\n");
      } 
    nIntervals_ = nintervals;
  }
  catch(ExceptionBase& e)
    { e.trace("in MainEffectsAnalyzer ctor"); }
}

// *******************************************************************
// If the user requested analysis for a specific input variable, verify
// that the variable is one available in the ddaceObj_.  If so, record
// the index of the input variable in inputIndex_, and save the name in
// inputNames_[0] (inputNames_, for a single input variable analysis,
// is an Array of length 1).  
//
// If all input variables are to be used in the analysis, the inputIndex_
// retains the value -1, and the inputNames_ Array contains all the 
// input names.
// -------------------------------------------------------------------
void MainEffectsAnalyzer::getMainEffectsInputNames(const String& inputName)
{
  try{
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    Array<String> varNames(sampler.nInputs());
    ddaceObj_.getVarNames(varNames);
    if(inputName != "ALL_INPUTS")
      {
	// Make sure that the input variable name (in inputName) 
	// matches with one of the input variable names in the ddaceObj_.
	// If so, assign inputName to inputName_, and record the index
	// in inputIndex_ (for locating the right column of input data
	// in the ddaceObj_).
	for(int i=0; i < varNames.length(); i++)
	  {
	    if(varNames[i] == inputName)
	      {
		inputIndex_ = i;
		inputNames_[0] = inputName;
		break;
	      }
	  }
	if(inputIndex_ < 0) ExceptionBase::raise("MainEffectsAnalyzer::"
						 "getMainEffectsInputNames() :"
						 " inputName provided does "
						 "not\nmatch any input names "
						 "in the ddaceObj.");
      }
    else
      {
	inputNames_.resize(varNames.length());
	inputNames_=varNames;
      }
  }
  catch(ExceptionBase& e)
    { e.trace("in MainEffectsAnalyzer::getMainEffectsInputNames()"); }
}

// ********************************************************************
// Get the input data corresponding to the one or more requested input
// variables for this analysis.
// --------------------------------------------------------------------
void MainEffectsAnalyzer::getData()
{
  try{
    inputData_.resize(nSamples_);

    // get the input and output (pts and ouputValues, respectively) 
    // arrays corresponding to the DDace object ddaceObj:
    Array<DDaceSamplePoint> pts(nSamples_);
    Array<Array<double> > outputValues(nSamples_);
    ddaceObj_.getResults(pts, outputValues);
    
    if(outputIndex_==-1)
      {
	ExceptionBase::raise("MainEffectsAnalyzer::getData() : user must "
			     "select\na single output variable for this "
			     "analysis.");
	exit(1);
      }
    
    // inputIndex_==1 : do analysis for all inputs, so collect input data
    // for all input varibles.
    if(inputIndex_==-1)
      {
	for (int j = 0; j < nSamples_; j++) 
	  {
	    inputData_[j].resize(pts[j].length());
	    for (int k = 0; k < pts[j].length(); k++)
	      inputData_[j][k] = pts[j][k];
	  }
      }
    // inputIndex_!=1 : do analysis for input in column inputIndex_
    // of the pts array.
    else
      {
	for (int j = 0; j < nSamples_; j++) 
	  {
	    // extract the column of input samples corresponding to
	    // the selected input variable, as well as the column of
	    // the output samples corresponding to the selected output
	    // variable.  Put the input in inputData_, and the output
	    // in outputData_.
	    
	    inputData_[j].resize(1);
	    inputData_[j][0] = pts[j][inputIndex_];
	  }
      }
  }
  catch(ExceptionBase& e)
    { e.trace("in MainEffectsAnalyzer::getData()"); }
}


// *******************************************************************
// generate group means & standard deviations of outputs, segregated
// by unique input (for each input).
// -------------------------------------------------------------------

void MainEffectsAnalyzer::calculateStats() 
{
  try
    {
      
      // ---------------------------------------------------------------
      // get the group means, standard deviations, unique x's.  Order
      // all according to the sorted unique x's.
      // ---------------------------------------------------------------
      
      

      // groupMean will hold the averaged value of the output (y) at a 
      // fixed value of the input (x)
      Array<Array<double> > groupMean(nIntervals_);
      Array<Array<double> > groupStdDev(nIntervals_);
      
      for(int i = 0; i < nIntervals_; i++)
	{
	  groupMean[i].resize(inputNames_.length());
	  groupStdDev[i].resize(inputNames_.length());
	}
      groupMean = calculateGroupMeans();
      groupStdDev = calculateGroupStdDev(groupMean);
      
      Array<Array<double> > xOrdered(nIntervals_);
      
      // first get all the unique X (input) values, and then we will
      // rearrange according to these values to get xOrdered.
      Array<Array<double> > uniqueX(nIntervals_);

      orderedUniqueInput_.resize(uniqueX.length());
      orderedGroupMeans_.resize(groupMean.length());
      orderedGroupStdDevs_.resize(groupStdDev.length());

      // get the nIntervals_ unique input values (the first nIntervals_
      // entries in the inputData_ array)
      for (int j = 0; j < nIntervals_; j++ ) 
	{
	  uniqueX[j].resize(inputData_[j].length());
	  uniqueX[j] = inputData_[j];
	}
      
      // sort the rows of the groupMean and groupStdDev by the increasing
      // ordering of the uniqueXs
      sortDouble(uniqueX, groupMean, groupStdDev);

      orderedUniqueInput_ = uniqueX;
      orderedGroupMeans_ = groupMean;
      orderedGroupStdDevs_ = groupStdDev;
      
    }
  catch(ExceptionBase& e)
    {
      e.trace("in MainEffectsAnalyzer::calculateStats()");
    }
}

// *******************************************************************
// generate group means 
// -------------------------------------------------------------------

Array<Array<double> > MainEffectsAnalyzer::calculateGroupMeans()
{
  try
    {
      int nreplications = nSamples_ / nIntervals_;
      Array<double> ySumByFixedX(nIntervals_);
      Array<double> fixedInputX(nIntervals_);
      Array<Array<double> > rtn(nIntervals_);
      int ncount;

      // carry out the the calculation for each of the i columns of the
      // inputData_ matrix.
      for (int i = 0; i < inputNames_.length(); i++)
	{
	  // carry out the calculation for each of the nIntervals_ unique
	  // input value in each column.
	  for (int j = 0; j < nIntervals_; j++ )
	    {
	      rtn[j].resize(inputNames_.length());
      
	      // start collecting yvalues for each x[j] (get first one)
	      ySumByFixedX[j] = outputData_[j];
	      fixedInputX[j] = inputData_[j][i];
	      
	      ncount = 1;
	      // NOTE: each replication for an input value, x, is nIntervals_
	      // long.  So after we've gone through the first nIntervals_ of
	      // the x-values, we start repeating these values in the next
	      // nIntervals_ segment of the total nSamples_.  However, each
	      // set of nIntervals_ values is a different permutation of the
	      // inputs.  Here, after we read in the first nIntervals_ of
	      // unique inputs, as we go through the remaining 
	      // (nSamples_ - nIntervals_), we add the y-value corresponding
	      // to a unique x-value to the appropriate element of 
	      // ySumByFixedX.
	      for (int k = nIntervals_; k < nSamples_; k++ )
		{

		  // if the current input is the same as the fixed input,
		  // add the corresponding output to the array item
		  // ySumByFixedX[j] that (currently) represents the sum 
		  // of the outputs when the input is held fixed at input 
		  // value j....

		  if ( abs(inputData_[k][i] - fixedInputX[j] ) < 1.0E-8 )
		    {
		      ySumByFixedX[j] += outputData_[k];
		      ncount++;
		    }
		}
	      if ( ncount != nreplications )
		{
		  ExceptionBase::raise("MainEffectsAnalyzer::calculateGroupMeans() :"
				       " samples are not \nreplicated LHS.");
		}
	      // Each indexed entry in ySumByFixedX[] is the sum of the 
	      // outputs for a fixed value of the input.  To get the mean 
	      // of the output for each fixed input, divide each element 
	      // of ySumByFixedX[] by the number of replications.
	      rtn[j][i] = ySumByFixedX[j] / (double) nreplications;
	    }
	}
      
      return rtn;
    }
  catch(ExceptionBase& e)
    {
      e.trace("in MainEffectsAnalyzer::calculateGroupMeans()");
    }
  return 0;
}

// *******************************************************************
// generate group standard deviations
// -------------------------------------------------------------------

Array<Array<double> > MainEffectsAnalyzer::calculateGroupStdDev(const Array<Array<double> >& groupMeans)
{
  try
    {
      int ncount;
      int nreplications = nSamples_ / nIntervals_;
      Array<double> ySqrdErrorSum(nIntervals_);
      Array<double> fixedInputX(nIntervals_);
      Array<Array<double> > rtn(nIntervals_);
      
      // carry out the the calculation for each of the i columns of the
      // inputData_ matrix.
      for (int i = 0; i < inputNames_.length(); i++)
	{
	  // carry out the calculation for each of the nIntervals_ unique
	  // input value in each column.
	  for (int j = 0; j < nIntervals_; j++)
	    {
	      rtn[j].resize(inputNames_.length());
	      ySqrdErrorSum[j] = pow((outputData_[j] - 
				      groupMeans[j][i]), 2);
	      fixedInputX[j] = inputData_[j][i];
	      
	      ncount = 1;
	      for(int k = nIntervals_; k < nSamples_; k++)
		{
		  // see calculateGroupMeans() (above) for further 
		  // explanation of the procedures in this function.
		  if ( abs(inputData_[k][i] - fixedInputX[j] ) < 1.0E-8 )
		    {
		      ySqrdErrorSum[j] += pow((outputData_[k] - 
					       groupMeans[j][i]), 2);
		      ncount++;
		    }
		}
	      if ( ncount != nreplications )
		{
		  ExceptionBase::raise("MainEffectsAnalyzer::calculateGroupStdDev() : "
				       "samples are not \nreplicated LHS.");
		}
	      // get the standard deviation of the output corresponding to
	      // the j-th unique input of the i-th input variable.
	      rtn[j][i] = sqrt(ySqrdErrorSum[j]/ (double) (nreplications - 1));
	    }
	}
      return rtn;
    }
  catch(ExceptionBase& e)
    {
      e.trace("in MainEffectsAnalyzer::calculateGroupStdDev()");
    }
  return 0;
}

// ********************************************************************
// For a one input variable analysis, return the name of the input 
// variable used in the analysis.
// --------------------------------------------------------------------
void MainEffectsAnalyzer::getInputName(String& inputName) const 
{ 
  try
    {
      if(inputIndex_ == -1)
	ExceptionBase::raise("MainEffectsAnalyzer::getInputName() : "
			     "there is more than one\ninput for this "
			     "analysis.  Use getInputNames().");
      else
	inputName = inputNames_[0];
    }
  catch(ExceptionBase& e)
    { e.trace("in MainEffectsAnalyzer::getInputName()"); }
}
  
// ********************************************************************
// Return the ordered, unique input values and the averaged outputs
// corresponding to the ordered inputs.
// --------------------------------------------------------------------

void MainEffectsAnalyzer::getGroupMeans(Array<Array<double> >& xOrdered,
		   Array<Array<double> >& groupMeans) const
{
  xOrdered = orderedUniqueInput_;
  groupMeans = orderedGroupMeans_;
}

// --------------------------------------------------------------------

void MainEffectsAnalyzer::getGroupMeans(String& inputName, 
				       Array<double>& xOrdered,
				       Array<double>& groupMeans) const
{
  int tmpIndex = -1;

  try
    {
      if (inputIndex_ != -1)
	{
	  if(inputName != inputNames_[0])
	    ExceptionBase::raise("MainEffectsAnalyzer::getGroupMeans() : "
				 "input name provided is not valid.");
	  else
	    for(int i=0; i < orderedUniqueInput_.length(); i++)
	      {
		xOrdered[i] = orderedUniqueInput_[i][0];
		groupMeans[i] = orderedGroupMeans_[i][0];
	      }
	}
      else
	{
	  for(int j=0; j < inputNames_.length(); j++)
	    {
	      if(inputName == inputNames_[j])
		{
		  tmpIndex = j;
		  break;
		}
	    }
	  if ((tmpIndex = -1))
	    ExceptionBase::raise("MainEffectsAnalyzer::getGroupMeans() : "
				   "input name provided is not valid.");
	  else 
	    for(int k=0; k < orderedUniqueInput_.length(); k++)
	      {
		xOrdered[k] = orderedUniqueInput_[k][tmpIndex];
		groupMeans[k] = orderedGroupMeans_[k][tmpIndex];
	      }
	}
    }
  catch(ExceptionBase& e)
    {
      e.trace("in MainEffectsAnalyzer::getGroupMeans()");
    }
}

// ********************************************************************
// Return the ordered, unique input values and the averaged outputs
// corresponding to the ordered inputs.
// --------------------------------------------------------------------

void MainEffectsAnalyzer::getGroupStdDevs(Array<Array<double> >& xOrdered,
		      Array<Array<double> >& groupStdDevs) const
{
  xOrdered = orderedUniqueInput_;
  groupStdDevs = orderedGroupStdDevs_;
}

// --------------------------------------------------------------------

void MainEffectsAnalyzer::getGroupStdDevs(String& inputName, 
					  Array<double>& xOrdered,
					  Array<double>& groupStdDevs) const
{
  int tmpIndex = -1;

  try
    {
      if (inputIndex_ != -1)
	{
	  if(inputName != inputNames_[0])
	    ExceptionBase::raise("MainEffectsAnalyzer::getGroupStdDevs() : "
				 "input name provided is not valid.");
	  else
	    for(int i=0; i < orderedUniqueInput_.length(); i++)
	      {
		xOrdered[i] = orderedUniqueInput_[i][0];
		groupStdDevs[i] = orderedGroupStdDevs_[i][0];
	      }
	}
      else
	{
	  for(int j=0; j < inputNames_.length(); j++)
	    {
	      if(inputName == inputNames_[j])
		{
		  tmpIndex = j;
		  break;
		}
	    }
	  if ((tmpIndex = -1))
	    ExceptionBase::raise("MainEffectsAnalyzer::getGroupStdDevs() : "
				   "input name provided is not valid.");
	  else 
	    for(int k=0; k < orderedUniqueInput_.length(); k++)
	      {
		xOrdered[k] = orderedUniqueInput_[k][tmpIndex];
		groupStdDevs[k] = orderedGroupStdDevs_[k][tmpIndex];
	      }
	}
    }
  catch(ExceptionBase& e)
    {
      e.trace("in MainEffectsAnalyzer::getGroupStdDevs()");
    }
  
}

// *******************************************************************
// sortDouble : sort the columns of array1 so that the elements of 
// each column are in ascending order.  The elements in array2 and 
// array3 are ordered so that if array element array1[i1][j] was 
// moved to array1[i3][j] in the reordering, element array2[i1][j]
// was moved to array2[i3][j] in the reordering of array2.  The 
// reordering of array3 will follow the same procedure as the
// reordering of array2.
// -------------------------------------------------------------------
void MainEffectsAnalyzer::sortDouble(Array<Array<double> >& array1, 
				     Array<Array<double> >& array2,
				     Array<Array<double> >& array3) const
{
  int    i, j, k, j1, leng;
  double dtmp;
  
  leng = array1.length();
  
  if (leng == 1) return;
  
  
  for (k=0; k < array1[0].length(); k++)
    {
      for (i=0; i<leng; i++) 
	{
	  for (j=1; j<leng-i; j++) 
	    {
	      j1 = j - 1;
	      // sort all three arrays relative to the sorting of 
	      // the first array (corresponding indices from each
	      // array reordered according to first array)
	      if (array1[j1][k] > array1[j][k]) 
		{
		  dtmp           = array1[j1][k];
		  array1[j1][k]  = array1[j][k];
		  array1[j][k]   = dtmp;
		  
		  dtmp           = array2[j1][k];
		  array2[j1][k]  = array2[j][k];
		  array2[j][k]   = dtmp;
		  
		  dtmp           = array3[j1][k];
		  array3[j1][k]  = array3[j][k];
		  array3[j][k]   = dtmp;
		}
	    }
	}
    }
}

// ********************************************************************
// mainEffectsToXML : create an XML MainEffectsAnalyzer object from
// the member data of the MainEffectsAnalyzer.
// --------------------------------------------------------------------
void MainEffectsAnalyzer::mainEffectsToXML(XMLObject& meXML) const
{
  try
    {
	XMLUtils::init();
	
	// The encompassing object: MainEffectAnalyzer
	XMLObject effect("MainEffectAnalyzer");
	
	// The imbedded archive object: Archive
	XMLObject archive("Archive");
	archive.addAttribute("name", archiveFile_);
	effect.addChild(archive);
	
	// The imbedded output name object: OutputName
	XMLObject outputName("Output");
	outputName.addAttribute("name", outputNames_[0]);
	effect.addChild(outputName);
	
	// NIntervals -- the number of intervals corresponding with the
	// number of unique input values in inputData_
	XMLObject nintervals("NIntervals");
	nintervals.addAttribute("value", Int(nIntervals_).toString());
	effect.addChild(nintervals);
	
	// Statistics, per unique, ordered input values:
  	XMLObject statistics("Statistics");
  	effect.addChild(statistics);

	// Embed a InputVariable object for each element of inputNames_
	for(int j=0; j < inputNames_.length(); j++)
	  {
	    XMLObject inVarStats("InputVariable");
	    inVarStats.addAttribute("name", inputNames_[j]);
	    statistics.addChild(inVarStats);

	    // Each InputVariable will have the following:
	    // (1) a tagged UniqueInput object for each of the ordered,
	    //     unique input values; 
	    // (2) an Input object which contains the unique input value; 
	    // (3) an Output object which holds the mean and standard 
	    //     deviation of the output, corresponding to the unique
	    //     input value.
	    for(int k=0; k < nIntervals_; k++)
	      {
		XMLObject uniqInp("UniqueInput");
		uniqInp.addAttribute("tag", Int(k).toString());
		inVarStats.addChild(uniqInp);

		XMLObject inpVal("Input");
		inpVal.addAttribute("value", Double(orderedUniqueInput_[k][j]).toString());
		XMLObject outVal("Output");
		outVal.addAttribute("mean", Double(orderedGroupMeans_[k][j]).toString());
		outVal.addAttribute("stdDev", Double(orderedGroupStdDevs_[k][j]).toString());
		uniqInp.addChild(inpVal);
		uniqInp.addChild(outVal);
	      }
	  }
  	meXML = effect;
    }
  catch(ExceptionBase& e)
    {
      e.trace("in MainEffectsAnalyzer::mainEffectsToXML()");
    }
  
}

