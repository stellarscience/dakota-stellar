#ifndef MAINEFFECTSEXCELOUTPUT_H_
#define MAINEFFECTSEXCELOUTPUT_H_

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
#include "ddace_config.h"
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <sstream>
#include "MainEffectsConverter.h"
#include "Factor.h"


class MainEffectsExcelOutput {
	public:
	
	    MainEffectsExcelOutput();
	    
	    ~MainEffectsExcelOutput();
	    
	    std::string computeExcelOutput
	        (std::vector<std::vector<double> > vectorInputData,
             std::vector<std::vector<double> > vectorOutputData);
             
    protected:
        std::string outputColumnHeaders
                    (int numInputs, int numOutputs);
                     
        std::string outputMainEffects
         (int inputVarIndex, int numInputs,
          int outputVarIndex, int numOutputs,
          DDaceMainEffects::Factor factor);
          
        std::string outputMainEffects
         (int inputVarIndex, int numInputs,
          int outputVarIndex, int numOutputs,
          DDaceMainEffects::Factor factor,
          int indexOfInputValue);          
};

#endif 
