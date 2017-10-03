#ifndef ONEWAYANOVA_H
#define ONEWAYANOVA_H

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
#include "ddace_config.h"
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <iostream>
#include <stdexcept>
#include "Factor.h"
#include "Response.h"
#include <sstream>


namespace DDaceMainEffects {

class OneWayANOVA
{
  public:
	/** C'tor for OneWayANOVA
	 * 	the constructor assigns factors and resp to internal variables
	 * 	after checking to make sure that the factors and the response are
	 *	both well formed.
	 *  @param: std::vector<Factor> factors, this represents the matrix of
	 *  factors. note that all of the factors must be coded as integers,
	 *  starting at zero.
	 */
	OneWayANOVA(std::vector<Factor> factors);
	
	

	/**
	 *	The getFactor(int i) function returns a Factor that is stored in
	 *  vector of factors that is stored as a private data member in this
	 *  class.
	 *  Note that this function is const to avoid accidental manipulation
	 *  of the factors_ data member.
	 *	@param: int i, this corresponds to the desired factor's location
	 *	in the vector of Factors (factors_). 
	**/ 
	const Factor getFactor(int i) const { return factors_[i]; };
	

	/**
	 *	Get an ANOVA table for one factor
	 *  @param int factor: this specifies which factor is used for the
	 *  ANOVA calculations.
	 *  @return ANOVA table
	**/
	std::string getANOVATable(int factor);	

	/**
	 *	Gets an ANOVA Table for each of the factors in order.
	 *  @param ANOVA Table of all the factors
	**/
	std::string getANOVATables();
	
	/**
	 *	The printANOVATable(int factor) function uses cout to output an
	 * 	ANOVA table for factor
	 *  @param int factor: this specifies which factor is used for the
	 *  ANOVA calculations.
	**/
	void printANOVATable(int factor);	

	/**
	 *	printANOVATable()
	 *	Prints out an ANOVA Table for each of the factors in order.
	**/
	void printANOVATables();	
	

  private:
	// This is the vector that holds all of the factors. 
	std::vector<Factor> factors_;
};

}//end namespace

#endif
