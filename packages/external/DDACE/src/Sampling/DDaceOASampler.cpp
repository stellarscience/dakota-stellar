#include "DDaceOASampler.h"

using namespace std;

std::string DDaceOASampler::typeName_ = "DDaceOASampler";

extern "C" {
  int bose_link( int n, int ninputs, int str, int ***AA );
  int OA_strength( int q,int nrow,int ncol,int** A,int *str,int verbose );
}


DDaceOASampler::DDaceOASampler(int nSamples, bool noise,
				 const std::vector<Distribution>& dist)
	: DDaceSamplerBase(nSamples, dist.size(), noise, dist),
		symbolMap_(0),
		nSymbols_(0)
{
	// The following code is adapted from Charles Tong's original
	// implementation

	// compute number of symbols
	double tmp = (double) nSamples_;
	tmp = pow(tmp, 0.5000001);
	nSymbols_ = (int) tmp;
	int nSamples1 = nSymbols_ * nSymbols_;
	if (nSamples1 < nSamples_)
		{
			int nSamples2 = (nSymbols_+1)*(nSymbols_+1);
			if ( (nSamples_ - nSamples1) < (nSamples2 - nSamples_))
				{
					nSamples_ = nSamples1;
				}
			else
				{
					nSamples_ = nSamples2;
					nSymbols_++;
				}
		}

	// set up sample pattern
	initPattern();
}

DDaceOASampler::DDaceOASampler(int nSamples, int nInputs, bool noise)
	: DDaceSamplerBase(nSamples, nInputs, noise),
		symbolMap_(0),
		nSymbols_(0)
{
	// The following code is adapted from Charles Tong's original
	// implementation

	// compute number of symbols
	double tmp = (double) nSamples_;
	tmp = pow(tmp, 0.5000001);
	nSymbols_ = (int) tmp;
	int nSamples1 = nSymbols_ * nSymbols_;
	if (nSamples1 < nSamples_)
		{
			int nSamples2 = (nSymbols_+1)*(nSymbols_+1);
			if ( (nSamples_ - nSamples1) < (nSamples2 - nSamples_))
				{
					nSamples_ = nSamples1;
				}
			else
				{
					nSamples_ = nSamples2;
					nSymbols_++;
				}
		}

	// set up sample pattern
	initPattern();

}

void DDaceOASampler::initPattern()
{
	int i,j ;
	// we need to use C arrays here to connect to the C code. 
	int** pTmp;
	int status = bose_link(nSamples_, nInputs_, 2, &pTmp);
	if (pTmp == 0) std::bad_alloc();

	if (status >= 0)
		{
			if (status != nSamples_)
				{
					cerr << "DDaceOASampler: num samples adjusted to " << status << endl;
					nSamples_ = status;
				}
		}
	else
		{
			throw runtime_error("DDaceOASampler::initPattern: bose cannot generate points");
		}
	
	
	// randomize the pattern


        std::vector<int> ivec(nSymbols_);
	for (i=0; i<nInputs_; i++)
		{
			ivec = randomIVector(nSymbols_);
			for (j=0; j<nSamples_; j++)
				{
					pTmp[j][i] = ivec[pTmp[j][i]];
				}
		}
	
	int strength;
	OA_strength( nSymbols_, nSamples_, nInputs_, 
                pTmp, &strength, 0);
	
	if (strength != 2)
		{
			throw runtime_error("Orthogonal Array Sampling : failed strength 2 test");
		}
	
	// copy into safe array and throw out C array
	symbolMap_.resize(nSamples_);
	for (i=0; i<nSamples_; i++)
		{
			symbolMap_[i].resize(nInputs_);
			for (j=0; j<nInputs_; j++)
				{
					symbolMap_[i][j] = pTmp[i][j];
				}
			delete[] pTmp[i];
		}
	delete[] pTmp;
}

/**
 * If a constructor which doesn't set the distribution is used to create the
 * DDaceOASampler object the getSamples() method will throw an
 * a runtime_error unless the setDist() method is used first to
 * pass in an array of distributions of length nInputs_.
 */
vector<DDaceSamplePoint>& DDaceOASampler::getSamples(std::vector<DDaceSamplePoint>& samplePoints) const
{
	double denom = (double) nSymbols_;
	samplePoints.resize(nSamples_);
	std::vector<double> x(nInputs_);
	double p;

  /**
   * Check that nInputs equals dist.size(), this must be true or else
   * there is the possibility of a seg fault in the method getSamples().
   */
        if (nInputs_ != (int) dist_.size())
            throw runtime_error("DDaceOASampler::getSamples: nInputs not equal to dist.length()");

	for (int s=0; s<nSamples_; s++)
		{
			for (int i=0; i<nInputs_; i++)
				{
					p = (double) symbolMap_[s][i];
					if (noise_)
						{
							p = (p + DistributionBase::uniformUnitDeviate())/denom;
						}
					else
						{
							p = (p + 0.5)/denom;
						}
					x[i] = dist_[i].getDeviate(p);
				}
			samplePoints[s] = DDaceSamplePoint(s, x);
		}
	return samplePoints;
}
	
DDaceSamplerBase* DDaceOASampler::clone() const
{
	DDaceSamplerBase* rtn = new DDaceOASampler(*this);
	if (rtn==0) throw std::bad_alloc();
	return rtn;
}

void DDaceOASampler::print(ostream& os) const
{
	os << "<OrthogonalArray "
		 << "samples=\"" << nSamples_ << "\" ";
	os << "perturb=\"";
	if (noise_)
		{
			os << "true\" ";
		}
	else
		{
			os << "false\" ";
		}
	os << "seed=\"" << DistributionBase::seed() << "\"/>";
}
	
