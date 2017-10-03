#include "Factor.h"

using std::vector;

namespace DDaceMainEffects {
	
Factor::Factor(){
		std::vector<int> emptyFactors;
		this->factors_ = emptyFactors; 
		this->nLevels_ = 0;
		Response emptyResponse;
		this->response_ = emptyResponse;
		vector<Response> emptyResponses;
		this->levelResponses_ = emptyResponses;
		this->nObs_ = 0;
}

Factor::Factor(std::vector<int> factors, int nLevels, Response response) 
	: response_(response), 
	  factors_(factors), 
	  nObs_(factors.size()), 
	  nLevels_(nLevels),
          //BUG:  assumes input data will have a value for every level           
	  //levelIndices_(nLevels, vector<int>()), 
	  		//BUG:  levelIndices_ & levelReponses_ are NOT parallel
          //levelIndices_(0),
          //BUG:  assumes input data will have a value for every level
	  //levelResponses_(nLevels)  
          levelResponses_(0)


{	
	// If nLevels is negative, then throw a runtime_error 
	if(nLevels_ <= 0) 
	  throw std::runtime_error
	  ("Error in Factor ctor: nLevels cannot be nonpositive.");
	  
	/*
	 * Break our (input,output) data pairs into groups of distinct
	 * input data values.  For example, if we have the following data set:
	 * {(5,1), (15,2), (10,3), (15,4), (5,5)} 
	 * Then we want to create the following groups 
	 * {all (input,output) data pairs where the input value is 5},
	 * {all (input,output) data pairs where the input value is 10},
	 * {all (input,output) data pairs where the input value is 15}.
	 * The output values of each group pair is saved in levelResponses_.
	 * levelResponses_ = {{1,5}, {3}, {2,4}}.
	 */
	for(int i=0; i< nLevels_; i++)
	{
		extractAllOutputDataValuesThatHaveThisInputDataValue(i);
	}

        /* BUG FIX:  Bad things happen if the input data does not */
        this->nLevels_ = levelResponses_.size();
}
	
	// Get the factors
	std::vector<int> Factor::getFactors(){
		return(this->factors_);
	}
	
    //Get the response
    Response Factor::getResponse() {
    	return(this->response_);
    }

void Factor::extractAllOutputDataValuesThatHaveThisInputDataValue
	     (int inputValue) {

	if(inputValue > nLevels_)
	{
		throw std::runtime_error("The specified level does not exist.");
	}


	vector<double> 	levelResp;
	vector<int>	currLevelIndices;

	for(int i = 0; i < nObs_; i++)
	{
	   // Our (input,output) data points are actually two parallel
           // arrays.  Factor contains the input data.        
           // Response contains the matching output data 
           // We are looking for a particular input data value.    
           // Whenever we find one, we save the corresponding output 
           // in levelResp
	   if(factors_[i]==inputValue)
	   {
	      // Save the input data point and output data point 
	      levelResp.push_back(response_[i]);
	      //currLevelIndices.push_back(i);
	   }
	}

        // levelIndices[i] contains a vector of the input data points 
        // whose value is i.  

	//BUG:  levelIndices_ & levelReponses_ are NOT parallel	
	//levelIndices_[lvl]=currLevelIndices;
          //if (currLevelIndices.size()==0) return;
          //levelIndices_.push_back(currLevelIndices);     

	/* BUG:  This code blows up if the vector levelResp is empty */
	//Response levelResponse(levelResp);
        if (levelResp.size()==0) return;
	Response levelResponse(levelResp);

	
	/* levelResponse[i] contains a vector of all the output data points */
	/* that were acquired when the input data point had a value of i.   */

        /* BUG:  This code blows up if any of the levelResp are empty. */
	//levelResponses_[lvl] = levelResponse;
        levelResponses_.push_back(levelResponse);
}

double Factor::sumOfSquaresBetweenGroups()
{
	double within = sumOfSquaresWithinGroups();
	return response_.getSumOfSquaresPop() - within;
}

double Factor::sumOfSquaresWithinGroups()
{	
	
	double within = 0;
	for(int i = 0; i < nLevels_; i++)
	{
		within += levelResponses_[i].getSumOfSquaresPop();
	}
	return within;
}

double Factor::varianceBetweenGroups()
{
	return sumOfSquaresBetweenGroups() / doFBetween();
}

double Factor::varianceWithinGroups()
{
	return sumOfSquaresWithinGroups() / doFWithin();
}

double Factor::Fdata()
{
	return varianceBetweenGroups() / varianceWithinGroups();
}

int Factor::getNumberOfLevels() const
{
	return nLevels_;
}

int Factor::getNumberOfObservations() const
{
	return nObs_;
}

int Factor::doFWithin() const
{
	return getNumberOfObservations() - getNumberOfLevels();
}

int Factor::doFBetween() const
{
	return getNumberOfLevels() - 1;
}

double Factor::getLevelVariance(int x)
{
	return levelResponses_[x].getVariancePop();
}

double Factor::getLevelSumOfSquares(int x)
{
	return levelResponses_[x].getSumOfSquaresPop();
}

double Factor::getLevelSum(int x)
{
	return levelResponses_[x].getSumPop();
}

double Factor::getLevelAverage(int x)
{
	return levelResponses_[x].getAveragePop();
}

vector<double> Factor::getAllLevelAverages()
{
	vector<double> rtn;
	for(int i = 0; i < getNumberOfLevels(); i++)
	{
		rtn.push_back(levelResponses_[i].getAveragePop());
	}
	return rtn;
}

}//namespace
