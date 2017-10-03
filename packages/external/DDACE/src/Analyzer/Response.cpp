#include "Response.h"

namespace DDaceMainEffects {
	

Response::Response(std::vector<double> resp) 
	: responses_(resp)
{
	
	// Check to make sure that the response vector that was passed in is
	// not empty.
	if(resp.empty()) 
	    throw std::runtime_error("Error in Response ctor: An empty vector was passed");
}

Response::Response(Response& other) : responses_(other.responses_)
{
}

Response::Response(const Response& other) : responses_(other.responses_)
{
}

double Response::getSumPop() const {	
	return Statistics::sum(responses_);
}

double Response::getAveragePop() const
{
	return Statistics::average(responses_);
}

double Response::getSumOfSquaresPop() const
{
	double mean = Statistics::average(responses_);
	return Statistics::sumOfSquares(responses_,mean);
}

double Response::getVariancePop() const
{
	return Statistics::variance(responses_);
}

double Response::operator[] ( unsigned int index ) const // note "const" here
{
    return responses_[index];
}
double & Response::operator[] ( unsigned int index )
{
    return responses_[index];
}

}//namespace
