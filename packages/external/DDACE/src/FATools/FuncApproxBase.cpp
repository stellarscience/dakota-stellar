#include "FuncApproxBase.h"

int FuncApproxBase::defaultNPtsPerDim_ = 10;

FuncApproxBase::FuncApproxBase()
  : nSamples_(1), 
    nInputs_(1), 
    nPtsPerDim_(defaultNPtsPerDim_),
    lowerBounds_(1, 0.0),
    upperBounds_(1, 0.0)
{;}


FuncApproxBase::FuncApproxBase(int nInputs, int nSamples)
	: nSamples_(nSamples), 
		nInputs_(nInputs), 
		nPtsPerDim_(defaultNPtsPerDim_),
		lowerBounds_(nInputs, 0.0),
		upperBounds_(nInputs, 0.0)
{
	;
}


// generate a 2D grid.

void FuncApproxBase::generate2DGridPoints(std::vector<std::vector<double> >& gridPts,
					  int var1, 
					  int var2,
					  const std::vector<double>& settings)
{
  
  int i, j; // stupid SGI doesn't meet standards for loop scoping.
  
  // grid size is N x N.
  gridPts.resize(nPtsPerDim_ * nPtsPerDim_);
  
  // set all coordinates to preset values
  for (i=0; i< (int) gridPts.size(); i++)
    {
      gridPts[i] = settings;
    }
  
  // now vary the two active coordinates
  
  double dx1 = upperBounds_[var1] - lowerBounds_[var1];
  double dx2 = upperBounds_[var2] - lowerBounds_[var2];
  
  for (i=0; i<nPtsPerDim_; i++)
    {
      double x1 = lowerBounds_[var1] 
	+ ((double) i)/((double) (nPtsPerDim_-1)) * dx1;
      for (j=0; j<nPtsPerDim_; j++)
	{
	  double x2 = lowerBounds_[var2] 
	    + ((double) j)/((double) (nPtsPerDim_-1)) * dx2;
	  gridPts[i*nPtsPerDim_ + j][var1] = x1;
	  gridPts[i*nPtsPerDim_ + j][var2] = x2;
	}
    }
}


// generate an N-dimensional grid
void FuncApproxBase::generateGridPoints(std::vector<std::vector<double> >& gridPts)
{
  throw std::runtime_error("FuncApproxBase::generateGridPoints not implemented yet");
}
