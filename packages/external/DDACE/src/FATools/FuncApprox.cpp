#include "FuncApprox.h"

FuncApprox::FuncApprox(const FuncApproxBase& base)
	: ptr_(base.clone())
{
	;
}

void FuncApprox::setBounds(const std::vector<double>& lowerBounds, 
			 const std::vector<double>& upperBounds)
{
	ptr_->setBounds(lowerBounds, upperBounds);
}

void FuncApprox::setNPtsPerDim(int nPtsPerDim)
{
	ptr_->setNPtsPerDim(nPtsPerDim);
}

int FuncApprox::getNPtsPerDim() const 
{
	return ptr_->getNPtsPerDim();
}

void FuncApprox::generateGridData(const std::vector<DDaceSamplePoint>& samplePts,
				const std::vector<double>& sampleResults,
				std::vector<std::vector<double> >& gridPts,
				std::vector<double>& gridResults)
{
	ptr_->generateGridData(samplePts, sampleResults, gridPts, gridResults);
}

void FuncApprox::generate2DGridData(const std::vector<DDaceSamplePoint>& samplePts,
				const std::vector<double>& sampleResults,
				int var1, 
				int var2,
				const std::vector<double>& settings,
				std::vector<std::vector<double> >& gridPts,
				std::vector<double>& gridResults)
{
	ptr_->generate2DGridData(samplePts, sampleResults, var1, var2, 
				 settings, gridPts, gridResults);
}


void FuncApprox::write2DGridData(const std::string& filename,
			 const std::vector<DDaceSamplePoint>& samplePts,
			 const std::vector<double>& sampleResults,
			 int var1, 
			 int var2,
			 const std::vector<double>& settings) 
{
	ptr_->write2DGridData(filename, samplePts, sampleResults, var1, var2, 
				settings);
}


double FuncApprox::evaluatePoint(const std::vector<double> & x) const
{
  return ptr_->evaluatePoint(x);
}

void FuncApprox::evaluatePoints(const std::vector<std::vector<double> >& pts, 
				std::vector<double>& y) const 
{
	ptr_->evaluatePoints(pts, y);
}

