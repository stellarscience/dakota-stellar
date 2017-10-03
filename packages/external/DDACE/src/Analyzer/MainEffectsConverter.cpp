#include "MainEffectsConverter.h"

    MainEffectsConverter::
        MainEffectsConverter() {}
    
    MainEffectsConverter::
        ~MainEffectsConverter() {}
    
    ValueAndRowIndexAndColumnIndex*
        MainEffectsConverter::convertTableOfDoublesToArray
            (std::vector<std::vector<double> >& vectorDoubles) {

        /* error check */
   	    if (vectorDoubles.size()==0) {
   	    	ValueAndRowIndexAndColumnIndex *values = 
   	            new ValueAndRowIndexAndColumnIndex[0];
   	        return(values);       	    	
   	    }
   	    
   	    /* how many doubles do we have? */
   	    int numberOfRows = vectorDoubles.size();
   	    int numberOfColumns = vectorDoubles[0].size();
   	    int numberOfDoubles = numberOfRows * numberOfColumns;


        /* create an empty array of ValueAndRowIndexAndColumnIndex */
   	    ValueAndRowIndexAndColumnIndex *values = 
   	        new ValueAndRowIndexAndColumnIndex[numberOfDoubles];
   	        
   	    /* copy the doubles from the input vector into our array */    
   	    int indexArray = 0;
   	    for (int indexRow=0; indexRow<numberOfRows; indexRow++) {
   	    	for (int indexCol=0; indexCol<numberOfColumns; indexCol++) {
   	    		values[indexArray++] = 
   	    		    ValueAndRowIndexAndColumnIndex
   	    		        (vectorDoubles[indexRow][indexCol],
   	    		         indexRow,
   	    		         indexCol);
   	    	}
   	    }
   	    
   	    /* return our array of ValueAndRowIndexAndColumnIndex */
   	    return(values);
    }

                          	
    VectorCountingNumbersAndCount MainEffectsConverter::
           convertAllDoublesToCountingNumbers
              (std::vector<std::vector<double> >& vectorDoubles) {
	
        /* error check */
   	    if (vectorDoubles.size()==0) {
   	    	VectorCountingNumbersAndCount x;
   	    	return(x);
   	    }
   	    
   	    /* how many doubles do we have? */
   	    int numberOfRows = vectorDoubles.size();
   	    int numberOfColumns = vectorDoubles[0].size();
   	    int numberOfDoubles = numberOfRows * numberOfColumns;
   	    
   	    
   	    /* 
   	     * copy all of the doubles 
   	     * into an array of ValueAndRowIndexAndColumnIndex 
   	     */
   	    ValueAndRowIndexAndColumnIndex *values = 
   	        convertTableOfDoublesToArray(vectorDoubles);
   	    
   	    
        /* sort the array */
        qsort
            (values,             //array to be sorted
             numberOfDoubles,    //size of array
             sizeof(ValueAndRowIndexAndColumnIndex), //size of one element
             ValueAndRowIndexAndColumnIndex::compare);  //compare function
                      
              
       
       /* We have to be very careful of round off errors.     */
       /* i.e. 9.999999, 10.000000, 10.0000001                */
       /* are all probably the same number.                   */
       /* We will say that 2 numbers are the same if          */
       /* they are very very very close to each other.        */
       double maxMinusMin = 
           values[numberOfDoubles-1].value - values[0].value;

       double averageDistance =
              maxMinusMin / (double)numberOfDoubles;
       double epsilon = averageDistance/ 100.0;   
              
       /* We are going to create a PARALLEL array of counting numbers. */
       /* That is, we are going to map each and every double in the    */
       /* input array to a counting number.                            */
       std::vector<std::vector<int> > vectorCountingNumbers;
       for (int i=0; i<numberOfRows; i++) {
          	vectorCountingNumbers.push_back(std::vector<int>(numberOfColumns));
       }
           
       
       
       /* 
        * Go ahead and map the first number in our array 
        * to the counting number 0.  How?  The first number in the array
        * was actually pulled from a row and column in the input vector.
        * Find the corresponding row and column in the output vector;
        * set that cell to the counting number 0.
        */
       int countingNumber = 0;
       double value = values[0].value;
       vectorCountingNumbers[values[0].indexRow][values[0].indexColumn] 
           = countingNumber;
                       

       /* process the remaing values in our array */ 
       for (int i=1; i<numberOfDoubles; i++) {
       	     	

          /* is it time to jump to the next counting number????? */
   	      if (fabs(values[i].value-value) > epsilon){ 
     	      value = values[i].value;
    	      countingNumber++; 	  
   	      }
   	         	      

   	      /* Map this value to a counting number. */
   	      /* How?  This value was pulled from a row and column */
   	      /* in the input vector.  Find the corresponding row and */
   	      /* column in the output vector; insert the value of the */
   	      /* count into that cell.                                */
          vectorCountingNumbers[values[i].indexRow][values[i].indexColumn] 
               = countingNumber;         
       }//for
   
   
       /* no memory leaks */
       delete[] values;
   
       /* We have successfully mapped all of the doubles that were */
       /* in the input vector into counting numbers.  The counting */
       /* numbers are in the output vector; this vector PARALLELS  */
       /* the input vector.  That is, the counting number that     */
       /* corresponds to the double value in the first row and     */
       /* first column can be found in the first row and first     */
       /* column of the output vector.                             */
       VectorCountingNumbersAndCount vectorCountingNumbersAndCount;
       vectorCountingNumbersAndCount.vectorCountingNumbers =
           vectorCountingNumbers;
       vectorCountingNumbersAndCount.count = countingNumber + 1;
       return(vectorCountingNumbersAndCount);
   
    }//method


    DDaceMainEffects::Factor 
         MainEffectsConverter::
             sliceOutOneInputVarAndOneOutputVar
	           (std::vector<std::vector<int> >& vectorInputDataPoints,   
	            std::vector<std::vector<double> >& vectorOutputDataPoints,  
	            int inputVarIndex,              
	            int outputVarIndex,
	            int numberOfValuesAvailableForOneInputVar){
	            	

        /* error check */
       DDaceMainEffects::Factor emptyFactor;
       if (vectorInputDataPoints.size()==0) return(emptyFactor);
       if (vectorOutputDataPoints.size()==0) return(emptyFactor);
       if (inputVarIndex<0) return(emptyFactor);
       if (outputVarIndex<0) return(emptyFactor);
       //does vector have at least inputVarIndex # of columns? */
       if (inputVarIndex>=vectorInputDataPoints[0].size()) {
           return(emptyFactor);
       }
       //does vector have at least outputVarIndex # of columns? */
       if (outputVarIndex>=vectorOutputDataPoints[0].size()) {
           return(emptyFactor);
       }        
       //does input var have at least one value to select from? */
       if (numberOfValuesAvailableForOneInputVar<=0)
           return(emptyFactor);
       //do the 2 vectors have the same number of rows (or samples or runs) */
       if (vectorInputDataPoints.size()!=vectorOutputDataPoints.size()) {
       	return(emptyFactor);
       }
       

        /* I need a vector of ints to hold the */
        /* data points from the selected input variable */
        std::vector<int> vectorDataFromSelectedInputVar;


        /* I need a vector of doubles to hold the */
        /* data points form the selected output variable */
        std::vector<double> vectorDataFromSelectedOutputVar;


        /* populate the two vectors 
         * i.e. copy the targeted input data column into 
         * vectorDataFromSelectedInputVar AND copy the   
         * tagetd output data column into
         * vectorDataFromSelectedOutputVar
         */
        int numberOfRuns = vectorInputDataPoints.size();
        for (int indexRun=0; indexRun<numberOfRuns; indexRun++) {
            vectorDataFromSelectedInputVar.push_back
                 (vectorInputDataPoints[indexRun][inputVarIndex]);
            double value = 
                vectorOutputDataPoints[indexRun][outputVarIndex];
            vectorDataFromSelectedOutputVar.push_back(value);                  
        }
        

        /* Create a DDACE MainEffects Response */
        DDaceMainEffects::Response ddaceMainEffectsResponse
                   (vectorDataFromSelectedOutputVar);



        /* Create a DDACE MainEffects Factor */
        DDaceMainEffects::Factor ddaceMainEffectsFactor
               (vectorDataFromSelectedInputVar, 
                numberOfValuesAvailableForOneInputVar, 
                ddaceMainEffectsResponse);



        /* return our DDACE MainEffects Factor */
        return(ddaceMainEffectsFactor);
        
    } 
   
	std::vector<DDaceMainEffects::Factor> 
	    MainEffectsConverter::convert
               (std::vector<std::vector<double> >& vectorInputDataPoints,
                std::vector<std::vector<double> >& vectorOutputDataPoints) {
                	
         /* error check */
         std::vector<DDaceMainEffects::Factor> emptyVector;
         if (vectorInputDataPoints.size()==0) return(emptyVector);
         if (vectorOutputDataPoints.size()==0) return(emptyVector);
               	
                	 
               
         /* Replace every INPUT data value with a counting number */

         VectorCountingNumbersAndCount vectorCountingNumbersAndCount =
              convertAllDoublesToCountingNumbers(vectorInputDataPoints);  
         std::vector<std::vector<int> > vectorInputIndicies =              
                vectorCountingNumbersAndCount.vectorCountingNumbers;
         int numberOfCountingNumbers = vectorCountingNumbersAndCount.count;

              
         /* Create an empty vector.  We will store our Factors here. */
         std::vector<DDaceMainEffects::Factor> vectorOfFactors;
         
        /* What do we have right now??????                           */
        /* vectorInputDataPoints contains the equivalent counting    */
        /* number values of every INPUT var for every run.           */
        /* vectorOutputDataPoints contains the data                  */
        /* values of every OUTPUT var for every run.                 */

        /* Let's try to merge the two vectors into something that    */
        /* DDACE's MainEffects can use.  i.e. create a vector of     */
        /* Factors.  Each element is a table.                        */
        /* Table 1 contains index values from column 1 of the input  */
        /*                  data values from column 1 of the output  */
        /* Table 2 contains index values from column 1 of the input  */
        /*                  data values from column 2 of the output  */        
        /* etc.                                                      */
        
        
        /* How many columns are in the input table? */
        int numberOfInputVariables = vectorInputDataPoints[0].size();
        
        /* How man columns are in the output table? */
        int numberOfOutputVariables = vectorOutputDataPoints[0].size();  
        
        
        /* pair input column 1 with output column 1 */
        /* pair input column 1 with output column 2 */
        /* pair input column 1 with output column 3 */
        /* etc.                                     */
        /* pair input column 2 with output column 1 */
        /* pair input column 2 with output column 2 */
        /* pair input column 2 with output column 3 */
        /* etc.                                     */        
        for (int indexInput=0; 
             indexInput<numberOfInputVariables; 
             indexInput++) {
        for (int indexOutput=0;
             indexOutput<numberOfOutputVariables; 
             indexOutput++) {


             /* slice out the selected input var & selected output var */
             DDaceMainEffects::Factor factor = 
                 sliceOutOneInputVarAndOneOutputVar
                      (vectorInputIndicies,    //data from all input vars
                       vectorOutputDataPoints, //data from all output vars
                       indexInput,             //slice out this input var
                       indexOutput,            //slice out this output var 
                       numberOfCountingNumbers);  //# of different input values

             /* append the slice to our vectorOfFactors */
             vectorOfFactors.push_back(factor);

        }//for indexOutput
        }//for indexInput
         
        /* return the vector of factors */
        return(vectorOfFactors);                                                       	
                	
    }
