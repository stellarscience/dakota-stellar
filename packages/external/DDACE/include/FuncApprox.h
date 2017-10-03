#ifndef FUNCAPPROX_H
#define FUNCAPPROX_H

#include <iostream>
#include "FuncApproxBase.h"
#include "SmartPtr.h"


/** 
 * FuncApprox is a handle class to FuncApproxBase
 * Author       : Charles Tong
 * Date         : September, 1997
 *
 * Modified by Kevin Long 4/22/1999
 */

class FuncApprox 
{
 public:
	
	FuncApprox(const FuncApproxBase& base);
	
	void setBounds(const std::vector<double>& lowerBounds, 
			 const std::vector<double>& upperBounds);

	void setNPtsPerDim(int nPtsPerDim) ;

	int getNPtsPerDim() const ;

	void generateGridData(const std::vector<DDaceSamplePoint>& samplePts,
				const std::vector<double>& sampleResults,
				std::vector<std::vector<double> >& gridPts,
				std::vector<double>& gridResults) ;

	void generate2DGridData(const std::vector<DDaceSamplePoint>& samplePts,
				const std::vector<double>& sampleResults,
				int var1, 
				int var2,
				const std::vector<double>& settings,
				std::vector<std::vector<double> >& gridPts,
				std::vector<double>& gridResults) ;
	
	 void write2DGridData(const std::string& filename,
			const std::vector<DDaceSamplePoint>& samplePts,
			const std::vector<double>& sampleResults,
			int var1, 
			int var2,
			const std::vector<double>& settings) ;
	 
	 double evaluatePoint(const std::vector<double> & x) const ;
	 void evaluatePoints(const std::vector<std::vector<double> >& pts, 
				 std::vector<double>& y) const ;
 protected:
	 SmartPtr<FuncApproxBase> ptr_;
};


#endif // FUNCAPRROX_H

