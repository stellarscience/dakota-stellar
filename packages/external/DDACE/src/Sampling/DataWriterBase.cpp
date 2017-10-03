#include "DataWriterBase.h"

DataWriterBase::DataWriterBase(const String& archiveFile)
  : ddaceObj_(), inputNames_(0), inputIndices_(0), 
  outputNames_(0), outputIndices_(0), archiveName_("")
{
  try {
    DDaceXMLReader reader(archiveFile);

    ddaceObj_ = reader.createObject();
    archiveName_ = archiveFile;
  }
  catch(ExceptionBase& e)
    { e.trace("in DataWriterBase::DataWriterBase(String&)"); }
}


void DataWriterBase::chooseInputVariable(const String& xVarName)
{
  try{
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);

    if( sampler.nInputs() == 0)
      ExceptionBase::raise("DataWriterBase:chooseInputVariable() : there are "
			   "\nno input variables for the specified "
			   "DDace object.");
    else
      {
	Array<String> varNames(sampler.nInputs());
	ddaceObj_.getVarNames(varNames);
    
	for(int i = 0; i < sampler.nInputs(); i++)
	  {
	    if(varNames[i] == xVarName)
	      {
		inputIndices_.resize(1);
		inputNames_.resize(1);
		
		inputIndices_[0] = i;
		inputNames_[0] = varNames[i];
		break;
	      }
	  }
	if (inputIndices_.length() == 0) 
	  ExceptionBase::raise("DataWriterBase::chooseInputVariable"
			       "(String&): input name provided does not"
			       "\nmatch any input names in the "
			       "DDace Object.");
      }
  }
  catch(ExceptionBase& e)
    { e.trace("in DataWriterBase::chooseInputVariable(String&)"); }
}


void DataWriterBase::chooseInputVariables(const String& xVarName, 
				     const String& yVarName)
{
  try{
    if (xVarName == yVarName) 
      ExceptionBase::raise("DataWriterBase::chooseInputVariables("
			   "Sring&, String&): you have selected\nthe "
			   "same input name for both inputs.  You "
			   "must select two unique names.");
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    if( sampler.nInputs() == 0)
      ExceptionBase::raise("DataWriterBase:chooseInputVariable() : there are "
			   "\nno input variables for the specified "
			   "DDace object.");
    else
      {
	Array<String> varNames(sampler.nInputs());
	ddaceObj_.getVarNames(varNames);
	int n = 0;
    
	for(int i = 0; i < sampler.nInputs(); i++)
	  {
	    if((varNames[i] == xVarName) || (varNames[i] == yVarName))
	      {
		if(n == 0)
		  {
		    inputIndices_.resize(2);
		    inputNames_.resize(2);
		  }
		inputIndices_[n] = i;
		inputNames_[n] = varNames[i];
		n++;
		if(n == 2)
		  break;
	      }
	  }
	if (n < 2) 
	  ExceptionBase::raise("DataWriterBase::chooseInputVariables"
			       "(String&, String&): at least\none of "
			       "the input names provided does not "
			       "match any input\nnames in the "
			       "DDace Object.");
      }
  }
  catch(ExceptionBase& e)
    { e.trace("in DataWriterBase::chooseInputVariables(String&, String&)"); }
}

void DataWriterBase::chooseAllInputVariables()
{
  try{
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    if(sampler.nInputs() == 0) 
      ExceptionBase::raise("DataWriterBase::chooseAllInputVariables(): "
			   "\nno inputs registered in the DDace object.");
    else
      {
	inputNames_.resize(sampler.nInputs());
	inputIndices_.resize(sampler.nInputs());
	ddaceObj_.getVarNames(inputNames_);
	
	for(int i = 0; i < sampler.nInputs(); i++)
	  {
	    inputIndices_[i] = i;
	  }
      }
  }
  catch(ExceptionBase& e)
    { e.trace("in DataWriterBase::chooseAllInputVariables()"); }
}

void DataWriterBase::chooseVariables(Array<String>& inputNames, 
				     Array<String>& outputNames)
{
  try {
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    if(sampler.nInputs() == 0) 
      ExceptionBase::raise("DataWriterBase::chooseInputVariables"
			   "(Array<String>&, Array<String>&): "
			   "\nno inputs registered in the DDace object.");
    else
      {
	inputNames_.resize(inputNames.length());
	inputIndices_.resize(inputNames.length());

	Array<String> varNames(sampler.nInputs());
	ddaceObj_.getVarNames(varNames);
	
	for(int i = 0; i < inputNames_.length(); i++)
	  {
	    for(int j = 0; j < sampler.nInputs(); j++)
	      {
		if(varNames[j] == inputNames[i])
		  {
		    inputIndices_[i] = j;
		    inputNames_[i] = varNames[j];
		    break;
		  }
	      }
	    ExceptionBase::raise("DataWriterBase::chooseVariables() : user "
				 "selected input/nvariable " + inputNames[i] +
				 " not found.");
	  }
      }

    if( outputNames.length() == 0)
      ExceptionBase::raise("DataWriterBase:chooseVariables() : there are "
			   "\nno output variables for the specified "
			   "DDace object.");
    else 
      {
	outputNames_.resize(outputNames.length());
	outputIndices_.resize(outputNames.length());
	
	Array<String> tmpOutNames(0);
	if(ddaceObj_.getNumOutputs() == 0)
	  ExceptionBase::raise("DataWriterBase:chooseVariables() : there are "
			       "\nno output variables for the specified "
			       "DDace object.");
	
	ddaceObj_.getOutputNames(tmpOutNames);
	
	for(int k = 0; k < outputNames.length(); k++)
	  {
	    for(int l = 0; l < tmpOutNames.length(); l++) 
	      {
		if(tmpOutNames[l] == outputNames[k])
		  {
		    outputIndices_[k] = l;
		    outputNames_[k] = tmpOutNames[l];
		    break;
		  }
	      }
	    ExceptionBase::raise("DataWriterBase::chooseVariables() : user "
				 "selected output/nvariable " + outputNames[k] +
				 " not found.");
	  }
      }
  }
  catch(ExceptionBase& e)
    { e.trace("in DataWriterBase::chooseVariables()"); }
}

void DataWriterBase::chooseAllVariables()
{
  try{
    DDaceSampler sampler;
    ddaceObj_.getSampler(sampler);
    if( sampler.nInputs() == 0)
      ExceptionBase::raise("DataWriterBase:chooseAllVariables() : there are "
			   "\nno input variables for the specified "
			   "DDace object.");
    inputNames_.resize(sampler.nInputs());
    ddaceObj_.getVarNames(inputNames_);
    if(ddaceObj_.getNumOutputs() == 0)
      ExceptionBase::raise("DataWriterBase:chooseAllVariables() : there are "
			   "\nno output variables for the specified "
			   "DDace object.");
    ddaceObj_.getOutputNames(outputNames_);
    
    inputIndices_.resize(inputNames_.length());
    outputIndices_.resize(outputNames_.length());
    
    for(int i = 0; i < sampler.nInputs(); i++) inputIndices_[i] = i;
    for(int j = 0; j < outputIndices_.length(); j++) outputIndices_[j] = j;
  }
  catch(ExceptionBase& e)
    { e.trace("in DataWriterBase::chooseAllVariables()"); }
}


void  DataWriterBase::chooseOutputVariable(const String& outputName)
{
  try {
    Array<String> tmpOutNames(0);
    if(ddaceObj_.getNumOutputs() == 0)
      ExceptionBase::raise("DataWriterBase:chooseOutputVariable() : there are "
			   "\nno output variables for the specified "
			   "DDace object.");
    ddaceObj_.getOutputNames(tmpOutNames);

    for (int i = 0; i <  tmpOutNames.length(); i++)
      {
	if(tmpOutNames[i] == outputName)
	  {
	    outputIndices_.resize(1);
	    outputNames_.resize(1);
	    
	    outputIndices_[0] = i;
	    outputNames_[0] = tmpOutNames[i];
	    break;
	  }
      }
    if (outputIndices_.length() == 0) 
      ExceptionBase::raise("DataWriterBase::chooseOutputVariable(String&): "
			   "output name\nprovided does not match any output "
			   "names in the DDace Object.");
  }
  catch(ExceptionBase& e)
    { e.trace("in DataWriterBase::chooseOutputVariable(String&)"); }
}






