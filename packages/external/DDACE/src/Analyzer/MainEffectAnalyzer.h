// *******************************************************************
// *******************************************************************
// Definition for the class ImportanceAnalyzer
// Author       : Charles Tong
// Date         : December, 1998
// *******************************************************************
// *******************************************************************

#ifndef _MAINEFFECTANALYZERH_
#define _MAINEFFECTANALYZERH_

#include <iostream>
#include <cmath>
#include "Array.h"
#include "String.h"

class MainEffectAnalyzer {

public:
   MainEffectAnalyzer()  {;}
   ~MainEffectAnalyzer() {;}

   void generateGroupMeans(int Nsamples, int Nintervals, 
                     const Array<double> x, const Array<double> y,
                     Array<double>& xout, Array<double>& ymean);

private:
   void computeMeanVariance(int Ninputs, int Noutputs, int Nsamples, 
                     const Array<double> y, Array<double>& amean,
                     Array<double>& avariance);

   void sortDouble(int leng, Array<double>& array, Array<double>& array2);

};

#endif 

