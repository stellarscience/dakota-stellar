//
// Tony Giunta
// 30 April 2001
// 01 May   2002, modified constructor arguments
//
//
// This is the Central Composite Design (CCD) sampling technique for 
// n-variables. CCD is intended for "square" or "rectangular" design
// spaces that have limits defined by upper and lower bounds.  Info
// on the CCD sampling method can be found in any classical design
// of experiments book (see reference below).  
//
// In most CCD descriptions, the n-dimensional design space is 
// normalized so that each variable (dimension) takes on values on
// (-1,+1).  So, a value of zero places a sample at the mid-point
// value of the variable, -1 is the lower bound, and +1 is the upper bound.
// In this implementation, this scaling is not used and the lower and
// upper bounds correspond to the user settings of the design parameters.
//
// CCD sampling selects the points in the following manner:
// (a) at the center of the n-dimensional square (i.e., hypercube), 
// (b) all vertices of the hypercube (which is a 2^n full-factorial 
//     pattern), and
// (c) the "axial" points that lie along the coordinate axes. 
//
// Note that the "axial" points are created by setting n-1 variables
// to zero, and then setting the remaining variable values to be
// +/- alpha, where "alpha" is typically given values between 1.0
// and Sqrt(n).  If alpha=1.0, the sampling method is called a 
// Face-Centered Central Composite Design. 
//
// *** This initial CCD implementation will be hardwired with a
// *** value of alpha=1.0.  Eventually, it would be good to read-in
// *** a user-supplied value of alpha (with 1.0 as default) for CCD 
// *** sampling.
//
// The total number of samples for CCD is:  1 + 2^n + 2*n
// where the 2*n quantity comes from the "axial" points.
//
//
// Reference: 
// R.H. Myers and D.C. Montgomery, "Response Surface Methodology: 
// Process and Product Improvement Using Designed Experiments," 
// Wiley, New York, 1995, pp. 297-318.


#include "DDaceCentralCompositeSampler.h"

std::string DDaceCentralCompositeSampler::typeName_ = 
  "DDaceCentralCompositeSampler";

DDaceCentralCompositeSampler::DDaceCentralCompositeSampler(int nSamples, 
  int nInputs, const std::vector<Distribution>& dist): 
  DDaceSamplerBase(nSamples, nInputs, false, dist)
{ 
  /* Check that nInputs equals dist.size(), this must be true or else
   * there is the possibility of a seg fault in the method getSamples().
   */
  if (nInputs != (int) dist.size())
    throw std::runtime_error("DDaceCentralCompositeSampler: nInputs not equal to dist.length()");

}

std::vector<DDaceSamplePoint>& DDaceCentralCompositeSampler::
getSamples(std::vector<DDaceSamplePoint>& samplePoints) const
{
  int i,j,k,currentIndex,s=0;
  samplePoints.resize(nSamples_);

  // extract the lower and upper bounds for the variables from "dist"
  std::vector<double> lower_bounds(nInputs_);
  std::vector<double> upper_bounds(nInputs_);
  for(i=0;i<nInputs_;i++){
    lower_bounds[i] = dist_[i].lowerBound();
    upper_bounds[i] = dist_[i].upperBound();
  }

  // the first sample is the "center point" of the n-dimensional hypercube
  std::vector<double> x(nInputs_);
  for(k=0; k<(nInputs_); k++) 
    x[k] = 0.5*(lower_bounds[k]+upper_bounds[k]);
  samplePoints[s] = DDaceSamplePoint(s, x);
  s++;

  // Next, create the "axial" sample points along the coordinate axes.
  for(i=0; i<(nInputs_); i++) {
    for(j=0; j<(nInputs_); j++)
      x[j] = 0.5*(lower_bounds[j]+upper_bounds[j]);
    x[i] = lower_bounds[i];
    samplePoints[s] = DDaceSamplePoint(s, x);
    s++;
    x[i] = upper_bounds[i];
    samplePoints[s] = DDaceSamplePoint(s, x);
    s++;
  }

  // Last, create the two-level full-factorial samples
  std::vector<double> stepSize(nInputs_);
  for(j=0; j<nInputs_; j++) {
    x[j]         = 0.0; 
    stepSize[j]  = upper_bounds[j]-lower_bounds[j];
  }
  currentIndex = nInputs_ - 1;
  //
  // fullFactorialPoints() function uses recursion to step through
  // all combinations of upper and lower values, for an arbitrary
  // number of variables.
  //
  fullFactorialPoints(currentIndex, nInputs_, lower_bounds, upper_bounds, 
		      stepSize, x, samplePoints, s);

   return samplePoints;
}

void DDaceCentralCompositeSampler::fullFactorialPoints(int currentIndex, 
  int numVars, const std::vector<double>& lower_bounds, 
  const std::vector<double>& upper_bounds, const std::vector<double>& stepSize, 
  std::vector<double>& x, std::vector<DDaceSamplePoint>& samplePoints, int& sIndex) const
{
  int nextIndex;
  int step;

  // Note: setting the "step" size in this for-loop limits this to a two-level
  // full factorial sampling method. For an n-level factorial, change the

  // values in the stepSize array and the step indices as needed.
  for( step=0; step<=1; step++ ) {
    x[currentIndex] = lower_bounds[currentIndex] + 
      static_cast<double>(step)*stepSize[currentIndex];
    if(currentIndex != 0) {
      nextIndex = currentIndex - 1; 
      fullFactorialPoints(nextIndex, numVars, lower_bounds, upper_bounds, 
		 stepSize, x, samplePoints, sIndex);
    }
    else {
      // If currentIndex == 0, then add the current array "x"
      // to the collection of DDace sample points.
      //
      samplePoints[sIndex] = DDaceSamplePoint(sIndex, x);
      sIndex++;
    }
  }
}

DDaceSamplerBase* DDaceCentralCompositeSampler::clone() const
{
  DDaceSamplerBase* rtn = new DDaceCentralCompositeSampler(*this);
  if (rtn==0)
    throw std::bad_alloc();
  return rtn;
}

void DDaceCentralCompositeSampler::print(std::ostream& os) const
{
  os << "METHOD Central Composite Design" << std::endl;
  os << "SAMPLES " << nSamples_ << std::endl;
}
std::vector<std::vector<int> > DDaceCentralCompositeSampler::getP() const 
{
        throw std::runtime_error("DDaceCentralCompositeSampler::getP not defined for base class");
        std::vector<std::vector<int> > tmp;
        return tmp;
}

