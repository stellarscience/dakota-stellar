#ifndef MainEffectsConverter_H_
#define MainEffectsConverter_H_

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
#include "ddace_config.h"
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include "OneWayANOVA.h"
#include "Factor.h"
#include "Response.h"
#include "ValueAndRowIndexAndColumnIndex.h"
#include "VectorCountingNumbersAndCount.h"


class MainEffectsConverter
{
public:
	MainEffectsConverter();
	virtual ~MainEffectsConverter();
	
     /**
      * Convert input data and output data into a format that
      * MainEffects can use.
      * <p>
      * vectorInputDataPoints[indexOfRun][indexOfVariable] 
      * contains a table of input data values.  Each input variable
      * is assigned a column.  Each row contains the data values from
      * one run.
      * <p>
      * EXAMPLE: if vectorInputDataPoints contains: <br>
      * <PRE> 
      *       variable1  variable2
      * run1: 11         12 
      * run2: 21         22
      * run3: 31         32
      * </PRE>
      * then <br>
      * vectorInputDataPoints[0][0] = 11 <br>
      * vectorInputDataPoints[0][1] = 12 <br>
      * vectorInputDataPoints[1][0] = 21 <br>
      * vectorInputDataPoints[1][1] = 22 <br>
      * vectorInputDataPoints[2][0] = 31 <br>
      * vectorInputDataPoints[2]10] = 32 <br>
      * <p>
      * vectorOutputDataPoints[indexOfRun][indexOfVariable] 
      * contains a table of output data values.  Each input variable
      * is assigned a column.  Each row contains the data values from
      * one run.
      * <p
      * EXAMPLE: if vectorOutputDataPoints contains: <br>
      * <PRE> 
      *       variable1  variable2
      * run1: 111        112 
      * run2: 121        122
      * run3: 131        132
      * </PRE>
      * then <br>
      * vectorInputDataPoints[0][0] = 111 <br>
      * vectorInputDataPoints[0][1] = 112 <br>
      * vectorInputDataPoints[1][0] = 121 <br>
      * vectorInputDataPoints[1][1] = 122 <br>
      * vectorInputDataPoints[2][0] = 131 <br>
      * vectorInputDataPoints[2]10] = 132 <br>
      * <p>
      * Each input data point is assigned an index number.
      * For example, the lowest number can be assigned the index number 0,
      * the next number can be assigned the index number 1, etc.
      * <p>
      * If we convert our example input data points into index values, 
      * then we would see:
      * <p>
      * <PRE> 
      *       variable1  variable2
      * run1: 0          3
      * run2: 1          4
      * run3: 2          5
      * </PRE>
      * <p>
      * This method returns vector of MainEffects Factors.
      * Each factor is a table: <br>
      * 1)The 1st column contains a column from the table of indices <br>
      * 2)The 2nd column contains the number of different indices <br>
      * 3)The 3rd column contains one of the columns from the output table <br> 
      * <p>
      * The number of factors is equal to the number of input data variables
      * times the number of output data variables.  In this example, the 
      * return value would contain these 4 tables: <br>
      * table 1: <br>
      * 0 6 111     0 6 112     3 6 111     3 6 112 <br>
      * 1 6 121     1 6 122     4 6 112     4 6 122 <br>
      * 2 6 131     2 6 132     5 6 113     5 6 132 <br>
      * @param vectorInputDataPoints
      * Contains the entire set of input data values.
      * Each element contains one run.  Each run is a vector.
      * Each element in the run contains one value for every input variable.
      * The vector must contain at least one run.
      * @param vectorOutputDataPoints
      * Contains the entire set of output data values.
      * Each element contains one run.  Each run is a vector.
      * Each element in the run contains one value for every output variable.
      * The vector must contain at least one run.
      * @return A vector of MainEffect's Factor.  Each factor contains
      * a) column of indices which point to one of the input columns, 
      * b) the number of values which are availble to the input variable, 
      * c) column of output values.
      */ 	
	std::vector<DDaceMainEffects::Factor> convert
	  (std::vector<std::vector<double> >& vectorInputDataPoints,
	   std::vector<std::vector<double> >& vectorOutputDataPoints);
                
                
    /**
     * Convert a 2-D vector of doubles into a 1-D array of doubles 
     * @param vectorDoubles The 2-D vector doubles that are to be 
     * converted into an array
     * @return The input vector of doubles converted into an array.
     */            
    ValueAndRowIndexAndColumnIndex *convertTableOfDoublesToArray
            (std::vector<std::vector<double> >& vectorDoubles);
                
    
	
    /** 
      * Replace each and every input data point        
      * with a counting number.  
      * <p>
      * Algorithm: <br>
      * Map the lowest value to the counting number 0. <br>
      * Map all numbers which are very close to the    
      * lowest value to the counting number 0. <br>
      * Map the lowest value that remains to the counting number 1. <br>
      * Map all numbers which are very close to the 
      * lowest value to the counting number 1. <br>
      * Map the lowest value that remains to the counting number 2. <br>
      * Map all numbers which are very close to the 
      * lowest value to the counting number 2. <br>
      * etc.
      * <p>
      * Example:
      * If the vector of doubles looks like this: <br>
      * &nbsp;&nbsp;&nbsp;&nbsp; 11  12   13 <br>
      * &nbsp;&nbsp;&nbsp;&nbsp; 21  22   23 <br>
      * &nbsp;&nbsp;&nbsp;&nbsp; 31  32   33 <br>
      * then we would return: <br>
      * &nbsp;&nbsp;&nbsp;&nbsp; 0   1   2 <br>
      * &nbsp;&nbsp;&nbsp;&nbsp; 3   4   5 <br>
      * &nbsp;&nbsp;&nbsp;&nbsp; 6   7   8 <br>
      * @param vectorDoubles Contains the entire set of input data values.
      * Each element contains one run.  Each run is a vector.
      * Each element in the run contains one value for every variable.
      * @return Replica of the input vector where all of the doubles
      * have been replaced with counting numbers.
      */  	
      VectorCountingNumbersAndCount convertAllDoublesToCountingNumbers
              (std::vector<std::vector<double> >& vectorDoubles);
                          
                          
    /**
      * Given the data from all the input variables for all runs
      * and given the data from all the output variables for all runs
      * and given the index value of one of the input variables
      * and given the index value of one of the output variables
      * then slice out the indexed variables (e.g. slice out
      * one column from the input table and slice out one column
      * from the output table).
      * <p>
      * Example:  If we have 3 input variables 
      * and if we have 5 runs
      * and if the input data values are: <br>
      * 1   2  3 <br>
      * 4   5  6 <br>
      * 7   8  9 <br>
      * 10 11 12 <br>
      * 13 14 15 <br>
      * and if the selected index is 1 
      * then we would slice out the 2nd column (i.e. {2, 5, 8 11, 14}).
      * @param vectorInputDataPoints The data from all input variables
      * for all runs.  Each element contains one run.  Each run 
      * contains one vector.  Each run element contains one value for
      * each input variable.  The vector must contains one or more runs.
      * @param vectorOutputDataPoints The data from all output variables
      * for all runs.  Each element contains one run.  Each run 
      * contains one vector.  Each run element contains one value for
      * each output variable.  The vector must contain one or more runs.
      * @param indexOfInputVariable Selects one of the input variables.
      * The first input variable that is inside vectorInputDataPoints has
      * an index value of 0.
      * @param indexOfOutputVariable Selects one of the output variables.
      * The first input variable that is inside vectorOutputDataPoints has
      * an index value of 0. 
      * @param numberOfValuesAvailableForOneInputVar The number of values
      * that are avaiable for the selected input variable.
      * @return The data that belongs to the selected input variable is
      * sliced out of vectorInputDataPoints.  The data that belongs to the
      * selected output variable is sliced out of vectorOutputDataPoints.
      */
    DDaceMainEffects::Factor 
         sliceOutOneInputVarAndOneOutputVar
           (std::vector< std::vector<int> >& vectorInputDataPoints,   
            std::vector< std::vector<double> >& vectorOutputDataPoints,  
            int indexOfInputVariable,              
            int indexOfOutputVariable,
            int numberOfValuesAvailableForOneInputVar);	    
            
            
    /**
      * Retrieve every data value from every input variable.  How?
      * Copy each and every data value that is in a 2-D table into a 1-D
      * table.
      * @param vectorInputDataPoints 
      * Contains the entire set of input data values.  The data is arranged
      * in a 2-D table.
      * @return Every value from every input variable.  The data is arranged
      * in a 1-D table.
      */  
    std::vector<double> extractAllInputDataValues
                   (std::vector<std::vector<double> >& vectorInputDataPoints);
                   
                   
	
};	

#endif 
