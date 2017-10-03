#ifndef RESPONSE_H
#define RESPONSE_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include "Statistics.h"

namespace DDaceMainEffects {

class Response
{
  public:
	/**
	 *  Default ctor for Response
	**/
	Response() {};

	/**
	 *  ctor for Response
	 *  The internal data member response_ is set by the resp param.
	 *  @param vector<double> resp: this is the data for the Responses
	 *  @throws runtime_error if resp is empty
	**/
	Response(std::vector<double> resp);

	/**
	 *	ctor for Response
	 *	This is the non-constant copy constructor. It takes the value of
	 * 	other's response and assigns it for itself.
	 *	@param Response& other: the response from which this will obtain
	 *	its data.
	**/
	Response(Response& other);

	/**
	 *	copy ctor for Response
	 *	This is the constant copy constructor. 
	 *	@param const Response other: the response from which 
         *      this will obtain its data
	**/
	Response(const Response& other);
	
	
	/**
	 *	double getSumPop() const
	 *	returns the sum of the responses. Uses the Statistics class
	 *	for all of the mathematics.
	**/	
	double getSumPop() const;
	
	/**
	 *	double getAveragePop() const
	 *	returns the average of the responses. Uses the Statistics class
	 *	for all of the mathematics.
	**/
	double getAveragePop() const;

	/**
	 *	double getSumOfSquaresPop() const
	 *	returns the sum of squares of the population. 
	**/
	double getSumOfSquaresPop() const;

	/**
	 *	double getVariancePop() const
	 *	returns the variance of the population. 
	**/
	double getVariancePop() const;
	
	/**
	 *	int getNumOfObservation() const
	 *	returns the number of observations in the responses.
	**/
	int getNumOfObservations() const { return responses_.size(); };	

	/**
	 *	double operator[](unsigned int index) const
	 *	this is the bracket operator. it allows for array like access to the
	 *	responses. This is the constant operator, used for reading from the
	 *	vector.
	 *	@param unsigned int index: index that you wish to access.
	**/
	double operator[](unsigned int index) const;

	/**
	 *	double & operator[](unsigned int index)
	 *	this is the bracket operator used to manipulate the internal response_
	 *	vector in an array-like manner.
	 *	Note: THIS ALLOWS FOR THE MANIPULATION OF THE response_ VECTOR THIS
	 *	CAN BE VERY DANGEROUS.
	 *	@param unsigned int index: index that you wish to change.
	**/
	double & operator[](unsigned int index);
  //private:
	
	// internal vector of doubles to hold the response data
	std::vector<double>	responses_;
};

}//namespace

#endif
