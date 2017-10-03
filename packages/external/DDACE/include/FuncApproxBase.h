/**
 * FuncApprox Base class
 * Author       : Charles Tong
 * Date         : September, 1997
 * Modified by Kevin Long 4/22/1999
 *
 */

#ifndef FUNCAPPROXBASE_H
#define FUNCAPPROXBASE_H

#include <iostream>
#include <stdexcept>
#include "DDaceSamplePoint.h"



class FuncApproxBase
{
 public:
  FuncApproxBase();
  FuncApproxBase(int nInputs, int nSamples);
  virtual ~FuncApproxBase() {;}
  
  virtual void setBounds(const std::vector<double>& lowerBounds, 
			 const std::vector<double>& upperBounds);
  
  virtual void setNPtsPerDim(int nPtsPerDim) ;
  
  virtual int getNPtsPerDim() const {return nPtsPerDim_;}
  
  virtual void generateGridData(const std::vector<DDaceSamplePoint>& samplePts,
				const std::vector<double>& sampleResults,
				std::vector<std::vector<double> >& gridPts,
				std::vector<double>& gridResults) = 0 ;
  
  virtual void generate2DGridData(const std::vector<DDaceSamplePoint>& samplePts,
				  const std::vector<double>& sampleResults,
				  int var1, 
				  int var2,
				  const std::vector<double>& settings,
				  std::vector<std::vector<double> >& gridPts,
				  std::vector<double>& gridResults) = 0 ;
  
  virtual void write2DGridData(const std::string& filename,
			       const std::vector<DDaceSamplePoint>& samplePts,
			       const std::vector<double>& sampleResults,
			       int var1, 
			       int var2,
			       const std::vector<double>& settings) = 0 ;
  
  virtual double evaluatePoint(const std::vector<double> & x) const = 0 ;
  virtual void evaluatePoints(const std::vector<std::vector<double> >& pts, 
			      std::vector<double>& y) const = 0 ;
  
  virtual FuncApproxBase* clone() const = 0;
 protected:

  /**
   * generate a grid on a 2D plane defined by two of our variables.
   * The other vars are set to values given in the "settings" array.
   */
  void generate2DGridPoints(std::vector<std::vector<double> >& gridPts,
			    int var1, 
			    int var2,
			    const std::vector<double>& settings);
  /// generate a grid in the Ninputs-dimensional sample space
  void generateGridPoints(std::vector<std::vector<double> >& gridPts);
  
  int nSamples_;
  int nInputs_;
  int nPtsPerDim_;
  std::vector<double> lowerBounds_;
  std::vector<double> upperBounds_;
  
  static int defaultNPtsPerDim_;
};



inline void FuncApproxBase::setBounds(const std::vector<double>& lowerBounds, 
				      const std::vector<double>& upperBounds)
{ 
  lowerBounds_ = lowerBounds;
  upperBounds_ = upperBounds;
}

inline void FuncApproxBase::setNPtsPerDim(int nPtsPerDim)   
{
  nPtsPerDim_ = nPtsPerDim;
}

#endif // FUNCAPRROXBASE_H


