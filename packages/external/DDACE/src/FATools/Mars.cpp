#include "Mars.h"

using namespace std;

int Mars::defaultNBasis_ = 15;
int Mars::defaultMaxVarPerBasis_ = 2;

Mars::Mars()
  :
  FuncApproxBase(),
  weights_(0),
  nBasis_(defaultNBasis_),
  maxVarPerBasis_(defaultMaxVarPerBasis_),
  varFlags_(0),
  fm_(0),
  im_(0)
{;}
 
Mars::Mars(int nInputs, int nSamples)
	: FuncApproxBase(nInputs, nSamples),
		weights_(nSamples, 1.0), 
		nBasis_(defaultNBasis_),
		maxVarPerBasis_(defaultMaxVarPerBasis_),
		varFlags_(nInputs, 1),
		fm_(0),
		im_(0)
{

	// set some magic numbers. See the Mars documentation for explanations
	int fmLength = 3 + nBasis_ * (5 * maxVarPerBasis_ + nSamples_ + 6) 
		+ 2 * nInputs_ + nSamples_;
	fm_.resize(fmLength);

	int imLength = 21 + nBasis_ * ( 3 * maxVarPerBasis_ + 8);
	im_.resize(imLength);
}



Mars::Mars(int nInputs, int nSamples, int nBasis, int maxVarPerBasis, 
					 const std::vector<int>& varFlags)
	: FuncApproxBase(nInputs, nSamples),
		weights_(nSamples_, 1.0), 
		nBasis_(nBasis),
		maxVarPerBasis_(maxVarPerBasis),
		varFlags_(varFlags),
		fm_(0),
		im_(0)
{

	// set some magic numbers. See the Mars documentation for explanations
	int fmLength = 3 + nBasis_ * (5 * maxVarPerBasis_ + nSamples_ + 6) 
		+ 2 * nInputs_ + nSamples_;
	fm_.resize(fmLength);

	int imLength = 21 + nBasis_ * ( 3 * maxVarPerBasis_ + 8);
	im_.resize(imLength);
}


void Mars::generateGridData(const std::vector<DDaceSamplePoint>& samplePts,
			const std::vector<double>& sampleResults,
			std::vector<std::vector<double> >& gridPts,
			std::vector<double>& gridResults) 
{
	if ((int) samplePts.size() != nSamples_)
		{
			throw runtime_error("mismatched sample point array size");
		}

	// compute the regression coefficients
	preProcess(samplePts, sampleResults);
	
	// fill in the array of grid points 
	generateGridPoints(gridPts);
	
	// evaluate the approximate function at the grid points
	evaluatePoints(gridPts, gridResults);
}

void Mars::generate2DGridData(const std::vector<DDaceSamplePoint>& samplePts,
			      const std::vector<double>& sampleResults,
			      int var1, 
			      int var2,
			      const std::vector<double>& settings,
			      std::vector<std::vector<double> >& gridPts,
			      std::vector<double>& gridResults) 
{
  if ((int) samplePts.size() != nSamples_)
    {
      throw runtime_error("mismatched sample point array size");
    }

  // compute the regression coefficients
  preProcess(samplePts, sampleResults);
	
	// fill in the array of grid points, varying only in the var1 and var2 
	// directions.
  generate2DGridPoints(gridPts, var1, var2, settings);
	
	// evaluate the approximate function at the grid points
	evaluatePoints(gridPts, gridResults);
}

void Mars::write2DGridData(const std::string& filename,
			 const std::vector<DDaceSamplePoint>& samplePts,
			 const std::vector<double>& sampleResults,
			 int var1, 
			 int var2,
			 const std::vector<double>& settings)
{
	int i;
	if ((int) samplePts.size() != nSamples_)
		{
			throw runtime_error("mismatched sample point array size");
		}

	// compute the regression coefficients
	//cerr << "preprocessing..." << endl;
	preProcess(samplePts, sampleResults);
	
	// fill in the array of grid points, varying only in the var1 and var2 
	// directions.
	std::vector<std::vector<double> > gridPts;
	//cerr << "generating point..." << endl;
	generate2DGridPoints(gridPts, var1, var2, settings);
	
	// evaluate the approximate function at the grid points
	std::vector<double> gridResults;
	//cerr << "evaluating point..." << endl;
	evaluatePoints(gridPts, gridResults);

	// write the grid points to file 

	ofstream os(filename.c_str());
	if (!os) throw runtime_error("Mars::write2DGridData could not open file");
		
	os << nPtsPerDim_ << endl;

	// print out the var1 coordinates
	double dx1 = upperBounds_[var1] - lowerBounds_[var1];
	for (i=0; i< nPtsPerDim_; i++)
		{
			os << lowerBounds_[var1] + ((double) i)/((double) (nPtsPerDim_-1))*dx1
				 << endl;
		}
	
	// print out the var2 coordinates
	double dx2 = upperBounds_[var2] - lowerBounds_[var2];
	for (i=0; i< nPtsPerDim_; i++)
		{
			os << lowerBounds_[var2] + ((double) i)/((double) (nPtsPerDim_-1))*dx2
				 << endl;
		}

	// print out the function values
	for (i=0; i< (int) gridResults.size(); i++)
		{
			os << gridResults[i] << endl;
		}
}






double Mars::evaluatePoint(const std::vector<double>& pt) const 
{
	std::vector<std::vector<double> > x(1, pt);
	std::vector<double> y;

	evaluatePoints(x, y);

	return y[0];
}


FuncApproxBase* Mars::clone() const 
{
	FuncApproxBase* rtn = new Mars(*this);
	if (rtn == 0) throw std::bad_alloc();
	return rtn;
}


//-------------------------------------------------------------------------
//
//              Interface to the MARS Fortran code
//
//-------------------------------------------------------------------------	

void Mars::preProcess(const std::vector<DDaceSamplePoint>& samplePts,
		      const std::vector<double>& sampleResults)
{
  int i, j;
  
  // This code just calls the fortran subroutine "mars." 
  // The only complication is the need to pack everything into 
  // stupid C arrays so we can talk to fortran.
  
  // fws, dws, and iws are workspace arrays for the call to fortran.
  // Aren't ya glad we don't write code in Fortran anymore?
  int length = nSamples_ * (nBasis_ + 4) + 6 * nSamples_ + 2 * nInputs_ 
    + 4 * nBasis_; 
  float* fws = new float[length];

  length = nSamples_ * nBasis_ + 10 * nBasis_;
  double* dws = new double[length];
  
  length = nSamples_ * nInputs_ + 2 * nSamples_;
  int* iws  = new int[length];
  
	
  // pack the sample points and sample results into C arrays.
  float* xlocal = new float[nSamples_ * nInputs_];
  float* ylocal = new float[nSamples_];
  
  for (i=0; i<nSamples_; i++) 
    {
      ylocal[i] = sampleResults[i];
      for (j=0; j<nInputs_; j++) 
	{
	  if (samplePts[i].length() != nInputs_)
	    {
	      throw runtime_error("Mars::preProcess bad sample point size");
	    }
	  xlocal[j*nSamples_ + i] = samplePts[i][j];
	}
    }

  // copy weights into C arrays 
  float* cWeights = new float[weights_.size()];
  if (cWeights==0) throw std::bad_alloc(); 
  for (i=0; i< (int) weights_.size(); i++) cWeights[i] = weights_[i];
  
	// copy var flags into C arrays 
  int* cVarFlags = new int[varFlags_.size()];
  if (cVarFlags==0) throw std::bad_alloc();
  for (i=0; i< (int) varFlags_.size(); i++) cVarFlags[i] = varFlags_[i];

  // allocate a C array for fm
  float* cFm = new float[fm_.size()];
  if (cFm==0) throw std::bad_alloc();

  // allocate a C array for im
  int* cIm = new int[im_.size()];
  if (cIm==0) throw std::bad_alloc();

  // Call the fortran code. As always with fortran everything
  // gets passed by address.
  MARS_F77(nSamples_, nInputs_, xlocal, ylocal, cWeights, nBasis_,
 	   maxVarPerBasis_, cVarFlags, cFm, cIm, fws, dws, iws);

  // now get the fm_ and im_ array data out of the stupid C form.
  
  for (i=0; i< (int) fm_.size(); i++) fm_[i] = cFm[i];
  for (i=0; i< (int) im_.size(); i++) im_[i] = cIm[i];
  
  // free the stupid C arrays
  
  delete [] fws;
  delete [] iws;
  delete [] dws;
  delete [] xlocal;
  delete [] ylocal;
  delete [] cWeights;
  delete [] cVarFlags;
  delete [] cFm;
  delete [] cIm;
}


void Mars::evaluatePoints(const std::vector<std::vector<double> >& pts,
			 std::vector<double>& y)
	const
{
	int i, j;
	int modelFlag = 2; // magic number specifying piecewise-cubic mars model

	// allocate C arrays

	int n = pts.size();
	y.resize(n);

	float* sp = new float[n * 4];
	if (sp==0) throw std::bad_alloc(); 
	
	float* xlocal = new float[n * nInputs_];
	if (xlocal==0) throw std::bad_alloc();
	
	float* ylocal = new float[n];
	if (ylocal==0) throw std::bad_alloc();

	// copy the input point coordinates into a C array
	for (i=0; i<n; i++)
		{
			for (j=0; j<nInputs_; j++)
				{
					xlocal[i+j*n] = pts[i][j];
					//					cerr << "point: " << i << " " << j << " " << xlocal[i+j*n] << endl;
				}
		}

	// allocate a C array for fm
	float* cFm = new float[fm_.size()];
	if (cFm==0) throw std::bad_alloc();
	for (i=0; i< (int) fm_.size(); i++) cFm[i] = fm_[i];	

	// allocate a C array for im
	int* cIm = new int[im_.size()];
	if (cIm==0) throw std::bad_alloc();
	for (i=0; i< (int) im_.size(); i++) cIm[i] = im_[i];
	
	//omit: cerr << "calling fmod..." << endl;
	FMODM_F77(modelFlag, n, xlocal, cFm, cIm, ylocal, sp);

	// omit: cerr << "copying points..." << endl;
	for (i=0; i<n; i++) y[i] = ylocal[i];

	delete [] sp;
	delete [] xlocal;
	delete [] ylocal;
	delete [] cFm;
	delete [] cIm;
}
