//
// This is a quasi-Box-Behnken sampling technique for n-variables. It
// creates a set of samples that encompass all two-variable combinations
// of n-variables, where each two-variable combination results in a 
// two-level factorial sample set.  Thus, the total number of samples
// is 4*(n \choose 2), where (n \choose 2) = n!/( 2! (n-2)! ) = 0.5*n*(n-1).
//
// The "quasi" part of the quasi-BB sampling method is due to the fact
// that for n > 5, BB sampling uses combination of three variables. 
// This three-variable combination is difficult to code up in a general
// manner, although someone more clever than I could do this at a later
// date....
//
// Reference: G.E.P.Box and D.W.Behnken, "Some new three level designs
// for the study of quantitative variables," Technometrics, Vol. 2, No. 4,
// 1960, pp. 455-475.
//
// Tony Giunta
// 19 March 2001
// 01 May   2002, modified constructor arguments
//


#include "DDaceBoxBehnkenSampler.h"

std::string DDaceBoxBehnkenSampler::typeName_ = "DDaceBoxBehnkenSampler";

DDaceBoxBehnkenSampler::DDaceBoxBehnkenSampler(int nSamples, int nInputs, 
  const std::vector<Distribution>& dist) :
  DDaceSamplerBase(nSamples, nInputs, false, dist) 
{
  /**
   * Check that nInputs equals dist.size(), this must be true or else
   * there is the possibility of a seg fault in the method getSamples().
   */
  if (nInputs != (int) dist.size())
    {
      throw std::runtime_error("DDaceBoxBehnkenSampler: nInputs not equal to dist.length()");
    }
}

std::vector<DDaceSamplePoint>& DDaceBoxBehnkenSampler::getSamples(std::vector<DDaceSamplePoint>& samplePoints) const
{
  int i,j,k,level,s=0;
  samplePoints.resize(nSamples_);
  
  // extract the lower and upper bounds for the variables from "dist"
  std::vector<double> lowerBnds(nInputs_);
  std::vector<double> upperBnds(nInputs_);
  for(i=0;i<nInputs_;i++){
    lowerBnds[i] = dist_[i].lowerBound();
    upperBnds[i] = dist_[i].upperBound();
  }

  // the first sample is the "center point" of the n-dimensional hypercube
  std::vector<double> x(nInputs_);
  for(k=0; k<(nInputs_); k++) 
    x[k] = 0.5*(lowerBnds[k]+upperBnds[k]);
  samplePoints[s] = DDaceSamplePoint(s, x);
  s++;

  // this creates all of the two-variable factorial samples
  for(i=0; i<=(nInputs_-2); i++) {
    for(j=i+1; j<=(nInputs_-1); j++) {
      for(level=1; level<=4; level++) {
	std::vector<double> x(nInputs_); // is a new 'x' needed each pass thru loop?
	for(k=0; k<(nInputs_); k++) 
	  x[k] = 0.5*(lowerBnds[k]+upperBnds[k]);
	if(level==1){
	  x[i] = upperBnds[i];
	  x[j] = upperBnds[j];
	}
	if(level==2){
	  x[i] = upperBnds[i];
	  x[j] = lowerBnds[j];
	}
	if(level==3){
	  x[i] = lowerBnds[i];
	  x[j] = upperBnds[j];
	}
	if(level==4){
	  x[i] = lowerBnds[i];
	  x[j] = lowerBnds[j];
	}
	samplePoints[s] = DDaceSamplePoint(s, x);
	s++;
      } // end For(level) loop
    }   // end For(j) loop
  }     // end For(i) loop
  
	return samplePoints;

}

DDaceSamplerBase* DDaceBoxBehnkenSampler::clone() const
{
	DDaceSamplerBase* rtn = new DDaceBoxBehnkenSampler(*this);
	if (rtn==0) throw std::bad_alloc();
	return rtn;
}

void DDaceBoxBehnkenSampler::print(std::ostream& os) const
{
        os << "METHOD BoxBehnken" << std::endl;
	os << "SAMPLES " << nSamples_ << std::endl;
}
	
std::vector<std::vector<int> > DDaceBoxBehnkenSampler::getP() const 
{
        throw std::runtime_error("DDaceBoxBehnkenSampler::getP not defined for base class");
        std::vector<std::vector<int> > tmp;
        return tmp;
}
	
