#ifndef DDACEANALYZER_H
#define DDACEANALYZER_H

#include "DDaceAnalyzerBase.h"

/**
 * DDaceAnalyzer:
 *
 */

class DDaceAnalyzer
{
 public:
  
  DDaceAnalyzer() : ptr_(0){;}
  
  DDaceAnalyzer(const DDaceAnalyzerBase& base) : ptr_(base.clone()){;}
  
  void getOutputName(String& outputName) const 
    { ptr_->getOutputName(outputName); }
  
  void getOutputNames(Array<String>& outputNames) const
    { ptr_->getOutputNames(outputName); }
  
  void getOutputData(Array<double>& outputData) const
    { ptr_->getOutputData(outputData); }

 protected:
  
  SmartPtr<DDaceAnalyzerBase> ptr_;
}
